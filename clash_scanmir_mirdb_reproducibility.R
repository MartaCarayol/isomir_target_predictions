library(biomaRt)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(reticulate)
library(scanMiR)
library(dplyr)

## First specific objective: 

# Access dataset hsapiens_gene_ensembl of the database ensembl:
mart <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Obtain specific information only for genes coding for proteins:
genes_of_interest_annot <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id",
                                              "uniprotswissprot", "hgnc_symbol", "chromosome_name", 
                                              "start_position", "end_position", "strand", "description",
                                              "gene_biotype", "refseq_mrna"), 
                                 filters="transcript_biotype", 
                                 values=list("protein_coding"),
                                 mart=mart) 

# Keep only those transcripts that are found on chromosomes 1-22 or on the X, Y or MT chromosomes:
genes_of_interest_annot <- subset(genes_of_interest_annot, chromosome_name %in% 
                                    c(1:22,'X', 'Y', 'MT'))

genes_of_interest_annot$chromosome_name <- paste0("chr", 
                                                  genes_of_interest_annot$chromosome_name)

# Keep transcripts that are annotated in UniProt and contain a RefSeq ID:
genes_of_interest_annot <- genes_of_interest_annot[genes_of_interest_annot$
                                                     uniprotswissprot != "" & 
                                                     genes_of_interest_annot$refseq_mrna != "",]

# Remove entries where refseq_mrna is associated with different ENSEMBL transcript ID: 
genes_of_interest_annot <- genes_of_interest_annot[!duplicated(genes_of_interest_annot$refseq_mrna),]

edb <- EnsDb.Hsapiens.v86
bsg <- BSgenome.Hsapiens.NCBI.GRCh38

# Get complete sequence of the transcripts:
all_sign_transcripts <- extractTranscriptSeqs(bsg, exonsBy(edb, 'tx', 
                                                           filter = AnnotationFilter::TxIdFilter(
                                                             genes_of_interest_annot$ensembl_transcript_id)))

# Check if the complete sequence of all protein-coding transcripts exists:
dim(as.data.frame(all_sign_transcripts))
dim(genes_of_interest_annot)

# Get 3'UTR sequence of the transcripts:
all_sign_transcripts_threeUTR <- extractTranscriptSeqs(bsg, threeUTRsByTranscript(edb, 
                                                                                  filter = AnnotationFilter::TxIdFilter
                                                                                  (genes_of_interest_annot$ensembl_transcript_id)))

# Check if the 3'UTR sequence of all complete sequences exists:
dim(as.data.frame(all_sign_transcripts_threeUTR))
dim(as.data.frame(all_sign_transcripts))

# Keep only transcripts for those 3'UTR sequence exists:
all_sign_transcripts <- all_sign_transcripts[names(all_sign_transcripts) %in% 
                                               names(all_sign_transcripts_threeUTR)]


# Get 10.1016/j.cell.2013.03.043 data:
clash <- read.table(file.path("C:/Users/marta/Desktop/masterBioinformaticaBioestadistica",
                              "/TFM/mmc1.txt"), header=TRUE)[,c(2,5,6,11,15,21,22,24)]

# Keep miRNA-mRNA interactions that have been classified into classes I-III:
clash_seed <- subset(clash, clash$folding_class %in% c('I','II','III'))

# Split the columns microRNA_name and mRNA_name:
clash_seed <- separate(clash_seed, col=microRNA_name, into=c("mirbase_ID", "X1", 
                                                             "mirna", "X2"), sep="_", remove=TRUE)[,c(-2, -4)]
clash_seed <- separate(clash_seed, col=mRNA_name, into=c("ENSG_ID", "ENST_ID", 
                                                         "hgnc_symbol", "X3"), sep="_", remove=TRUE)[,c(-7)]

# Keep miRNA-mRNA interactions where mRNA is in all_sign_transcripts:
clash_seed <- subset(clash_seed, clash_seed$ENST_ID %in% 
                       subset(genes_of_interest_annot$ensembl_transcript_id, 
                              genes_of_interest_annot$ensembl_transcript_id %in%
                                names(all_sign_transcripts)))

# Keep miRNA-mRNA interactions where log2_target_enrichment is not NA:
clash_seed_log2 <- subset(clash_seed, !is.na(clash_seed$log2_target_enrichment))

# Keep names of the 15 miRNA with the highest number of targets:
mirna_max_obj <- tail(table(clash_seed_log2$mirna)[order(table(clash_seed_log2$mirna))],15)

# Calculate the average log2 fold change of the 15 miRNAs with the highest number of targets: 
means_log2 <- c()

for (i in names(mirna_max_obj)){
  means_log2 <- c(means_log2, mean(clash_seed_log2$
                                     log2_target_enrichment[clash_seed_log2$mirna == i]))
}

# Keep miRNA-mRNA interactions for the 5 miRNA with the highest average log2 fold change:
clash_seed_sub <- subset(clash_seed_log2, clash_seed_log2$mirna %in% 
                           names(mirna_max_obj[tail(order(means_log2),5)]))

unique(clash_seed_sub$mirna)

# Create a .txt file with the miRNA in FASTA format for miRDB use: 
for (mirna in unique(clash_seed_sub$mirna)){
  write(c(paste0(">", mirna), unique(clash_seed_sub$miRNA_seq[clash_seed_sub$mirna==mirna])),
        "C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM/data_clash_objetivo1.txt",
        append=TRUE)
}

# Create a .csv file for scanMiR use:
write.csv(data.frame("mirna"=unique(clash_seed_sub$mirna), "miRNA_seq"=unique(clash_seed_sub$miRNA_seq)),
          "C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM/data_clash_objetivo1.csv",
          row.names=FALSE)

# Set directory where the scripts are:
setwd('C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM')

# Run scanMiR script:
system('Rscript script_scanmir.R C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM data_clash_objetivo1.csv 200 result_scanmir C:/Users/marta/anaconda3/python.exe')

# Select Python version: 
use_python("C:\\Users\\marta\\anaconda3\\python.exe")
py_config()

# Run Python script:
system('python automated_mirdb_mirna.py data_clash_objetivo1.txt result_mirdb Human', wait=TRUE)

## Second specific objective: 

result_mirdb <- read.table(file.path("C:/Users/marta/Desktop",
                                     "/masterBioinformaticaBioestadistica/TFM/result_mirdb.csv"), 
                           header=TRUE, sep=',')

result_scanmir <- read.table(file.path("C:/Users/marta/Desktop",
                                       "/masterBioinformaticaBioestadistica/TFM/result_scanmir.csv"), 
                             header=TRUE, sep=',')

# Keep miRNA-mRNA interactions where repression < -1.5:
result_scanmir <- as.data.frame(result_scanmir %>%
                                  filter(repression < -1.5))

# Keep miRNA-mRNA interactions where mRNA is in all_sign_transcripts:
result_scanmir <- subset(result_scanmir, result_scanmir$transcript %in%
                           subset(genes_of_interest_annot$ensembl_transcript_id, 
                                  genes_of_interest_annot$ensembl_transcript_id %in%
                                    names(all_sign_transcripts)))

# Keep miRNA-mRNA interactions where mRNA is in all_sign_transcripts:
result_mirdb <- subset(result_mirdb, result_mirdb$GenBankAccession %in% 
                         subset(genes_of_interest_annot$refseq_mrna,
                                genes_of_interest_annot$ensembl_transcript_id %in%
                                  names(all_sign_transcripts)))

# Add ensembl_transcript_id column:
ensembl_transcript_id <- c()
for(i in result_mirdb$GenBankAccession){
  ensembl_transcript_id <- c(ensembl_transcript_id, 
                             subset(genes_of_interest_annot$ensembl_transcript_id, 
                                    genes_of_interest_annot$refseq_mrna == i))
}

result_mirdb$ensembl_transcript_id <- ensembl_transcript_id

# Keep targets obtained by scanMiR that have at least a target site of type 8mer or 7mer in the 3'UTR region:
result_scanmir_3UTR <- subset(result_scanmir, result_scanmir[,3]!=0 | result_scanmir[,4]!=0)

# Keep targets experimentally obtained that have at least a target site of type 8mer, 7mer or 9mer in the 3'UTR region:
clash_seed_sub_3UTR <- subset(clash_seed_sub, !is.na(clash_seed_sub$X3.UTR) & 
                             (clash_seed_sub$seed_type == "8-mer" | 
                              clash_seed_sub$seed_type == "7-mer" | 
                              clash_seed_sub$seed_type == "9-mer"))

# Calculate percentage of the results obtained with different approaches.
# args:
#   - df_x_mirna_enst_id: first argument for semi_join().
#   - df_y_mirna_enst_id: second argument for semi_join().
# return: percentage of the results obtained with df_x_mirna_enst_id that have been
# obtained in df_y_mirna_enst_id.
percentage_results <- function(df_x_mirna_enst_id, df_y_mirna_enst_id){
  
  names(df_x_mirna_enst_id)[1] <- "mirna"
  names(df_y_mirna_enst_id)[1] <- "mirna"
  names(df_x_mirna_enst_id)[2] <- "enst_id"
  names(df_y_mirna_enst_id)[2] <- "enst_id"
  
  result_inner_join <- semi_join(df_x_mirna_enst_id, df_y_mirna_enst_id, 
                                 by=c("mirna", "enst_id"))
  
  return((dim(result_inner_join)[1]/dim(df_x_mirna_enst_id)[1])*100)
}

# Percentage of the results obtained with scanMiR that have been obtained experimentally:
percentage_results(result_scanmir[,c(9,1)], clash_seed_sub[,c(2,5)])

# Percentage of the results experimentally obtained that have been obtained with scanMiR:
percentage_results(clash_seed_sub[,c(2,5)], result_scanmir[,c(9,1)])

# Percentage of the results obtained with scanMiR that have been obtained with miRDB:
percentage_results(result_scanmir_3UTR[,c(9,1)], result_mirdb[,c(2,7)])

# Percentage of the results obtained with miRDB that have been obtained with scanMiR:
percentage_results(result_mirdb[,c(2,7)], result_scanmir_3UTR[,c(9,1)])

# Percentage of the results obtained with miRDB that have been obtained experimentally:
percentage_results(result_mirdb[,c(2,7)], clash_seed_sub_3UTR[,c(2,5)])

# Percentage of the results obtained experimentally that have been obtained with miRDB:
percentage_results(clash_seed_sub_3UTR[,c(2,5)], result_mirdb[,c(2,7)])
