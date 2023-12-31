library(biomaRt)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(EnsDb.Hsapiens.v86)
library(reticulate)
library(scanMiR)
library(foreach)
library(doParallel)
library(doSNOW)

get_transcripts_ensembl_id <- function(){
  
  # Be careful with dbplyr and BiocFileCache packages!!
  
  print('Getting ensembl_transcript_id of protein coding transcripts...')
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  genes_of_interest_annot <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id",
                                                "uniprotswissprot", "hgnc_symbol", "chromosome_name",
                                                "gene_biotype"), 
                                   filters="transcript_biotype", 
                                   values=list("protein_coding"),
                                   mart=mart)
  # To compare with miRDB results, we use: genes_of_interest_annot <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "uniprotswissprot", "hgnc_symbol", "chromosome_name", "gene_biotype", "refseq_mrna"), filters="transcript_biotype", values=list("protein_coding"), mart=mart)
  
  genes_of_interest_annot <- subset(genes_of_interest_annot, chromosome_name %in% c(1:22, 'X', 'Y', 'MT'))

  # Get transcripts annotated in UniProt:
  genes_of_interest_annot <- genes_of_interest_annot[genes_of_interest_annot$uniprotswissprot != "",]
  # To compare with miRDB, we use: genes_of_interest_annot <- genes_of_interest_annot[genes_of_interest_annot$uniprotswissprot != "" & genes_of_interest_annot$refseq_mrna != "",]
  # And: genes_of_interest_annot <- genes_of_interest_annot[!duplicated(genes_of_interest_annot$refseq_mrna),]
  
  return(genes_of_interest_annot)
}

get_transcripts_seqs <- function(is_threeUTR, genes_of_interest_annot){
  
  edb <- EnsDb.Hsapiens.v86
  bsg <- BSgenome.Hsapiens.NCBI.GRCh38

  # 3'UTR sequences of all protein-coding transcripts (U are replaced by T):
  all_sign_transcripts_threeUTR <- extractTranscriptSeqs(bsg, threeUTRsByTranscript(edb, filter = AnnotationFilter::TxIdFilter(
                                                                                    genes_of_interest_annot$ensembl_transcript_id)))
  if(is_threeUTR == TRUE){
    print('Getting 3UTR sequence of protein coding transcripts...')
    return(all_sign_transcripts_threeUTR)
  }else{
    print('Getting complete sequence of protein coding transcripts...')
    # Complete sequences of all protein-coding transcripts (U are replaced by T):
    all_sign_transcripts <- extractTranscriptSeqs(bsg, exonsBy(edb, 'tx', filter = AnnotationFilter::TxIdFilter(
                                                                          genes_of_interest_annot$ensembl_transcript_id)))
    
    # Select transcripts included in all_sign_transcripts_threeUTR:
    all_sign_transcripts <- all_sign_transcripts[names(all_sign_transcripts) %in% names(all_sign_transcripts_threeUTR)]
    
    return(all_sign_transcripts)
  }
}

init_python <- function(python_path){
  use_python(python_path) # 'C:/Users/marta/anaconda3/python.exe'
  py_config()
}

# Calculate length of 5'UTR+CDS (named ORF).
length_ORF <- function(dna_string_set_transcripts, dna_string_set_threeUTR){
  
  print('Calculating length of ORF...')
  
  length_ORF <- foreach (i=1:length(dna_string_set_threeUTR), .combine=c, .packages='Biostrings') %dopar% {
    
    width(dna_string_set_transcripts[names(dna_string_set_transcripts)==
                                                       names(dna_string_set_threeUTR[i])])-width(dna_string_set_threeUTR[i])
  }
}

# Return KdModelList using CNN from McGeary et al. (2019) (https://github.com/kslin/miRNA_models/)
generate_KdModelList <- function(mirna_seq){
  
  foreach(i=1:length(mirna_seq[,1]), .packages=c('reticulate', 'scanMiR')) %dopar% {
    
    system(sprintf('python generate_12mer_kds.py --name %s --mirseq %s --mirlen 10 
                   --load_model model-100 --outfile kd_%s.txt', 
                   mirna_seq[i,1], mirna_seq[i,2], mirna_seq[i,1]), wait=TRUE)
    
    obj_KdModel_i <- getKdModel(kd=sprintf('kd_%s.txt', mirna_seq[i,1]), mirseq=mirna_seq[i,2], name=mirna_seq[i,1])
    
    file.remove(sprintf('kd_%s.txt', mirna_seq[i,1]))
    
    return(obj_KdModel_i)
  }
}

search_aggr_select_objectives <- function(obj_KdModelList, all_sign_transcripts, mirna_seq, rows_number){
  
  print('Searching, aggregating and selecting...')

  # Progress bar settings:
  pb <- txtProgressBar(max=length(obj_KdModelList), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  result_scanmir <- foreach(i=1:length(obj_KdModelList), .packages='scanMiR', .options.snow=opts, .verbose=T, .combine=rbind) %dopar% {

    search_i <- findSeedMatches(seqs=all_sign_transcripts, seeds=obj_KdModelList[[i]], verbose=TRUE, ret='data.frame')
    aggr_i <- aggregateMatches(search_i)
    select_i <- aggr_i[order(aggr_i$repression)[1:rows_number],]
    select_i$mirna <- rep(mirna_seq[i,1], rows_number)
    
    select_i
  }
  
  close(pb)
  return(result_scanmir)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  args_length <- length(args)
  if(args_length > 4){stop('error: too many arguments.')} 
  if(args_length < 3){stop('error: missing arguments.')}
  if(is.numeric(as.numeric(args[3])) == FALSE){stop('error: the amount of targets must be a number')}
  
  # Initializing working directory and Python
  setwd(args[1])
  if(args_length == 4){init_python(args[4])}
  
  # Parallel processing configuration:
  n_cores <- round(detectCores()/2)
  cl <- makeCluster(n_cores, type="SOCK")
  registerDoSNOW(cl)
  
  genes_of_interest_annot <- get_transcripts_ensembl_id()
  all_sign_transcripts <- get_transcripts_seqs(FALSE, genes_of_interest_annot)
  
  # Add a column named ORF.length in all_sign_transcripts object (important to findSeedMatches()):
  mcols(all_sign_transcripts)$ORF.length <- length_ORF(all_sign_transcripts, get_transcripts_seqs(TRUE, genes_of_interest_annot))

  mirna_seq <- read.csv(args[2], header=TRUE, sep=',')
  
  KdModelList_interest <- generate_KdModelList(mirna_seq)
  
  result_scanmir <- search_aggr_select_objectives(KdModelList_interest, all_sign_transcripts, mirna_seq, args[3])
  
  write.csv(result_scanmir, file="result_scanmir_azoos.csv", row.names=FALSE)
  
  stopCluster(cl)
}


main()


