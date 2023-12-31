library(biomaRt)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(reticulate)


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

read.csv("C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM/data_clash_objetivo1.csv", header=TRUE, sep=',')

# Set directory where the scripts are:
setwd('C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM')

# Run scanMiR script:
system('Rscript script_scanmir.R C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM data_clash_objetivo1.csv 200 C:/Users/marta/anaconda3/python.exe')

# 
use_python("C:\\Users\\marta\\anaconda3\\python.exe")
py_config()

# Ejecutamos el script de Python:
system('python automated_mirdb_mirna.py data_clash_objetivo1.txt result_mirdb Human',
       wait=TRUE)

# REPRODUCIBILIDAD

# Cargamos los resultados obtenidos con miRDB:
result_mirdb <- read.table(file.path("C:/Users/marta/Desktop",
                                     "/masterBioinformaticaBioestadistica/TFM/result_mirdb.csv"), 
                           header=TRUE, sep=',')

# Tenemos 825 objetivos predichos para los 5 miRNA:
dim(result_mirdb)

# Nos quedamos solo con los transcritos codificantes de proteínas que se encuentran
# en all_sign_transcripts:
result_mirdb <- subset(result_mirdb, result_mirdb$GenBankAccession %in% 
                         subset(genes_of_interest_annot$refseq_mrna,
                                genes_of_interest_annot$ensembl_transcript_id %in%
                                  names(all_sign_transcripts)))

# Estos 5 miRNA regulan a 753 transcritos codificantes de proteínas:
dim(result_mirdb)

# Añadimos la columna ensembl_transcript_id a los resultados obtenidos utilizando
# miRDB:
ensembl_transcript_id <- c()
for(i in result_mirdb$GenBankAccession){
  ensembl_transcript_id <- c(ensembl_transcript_id, 
                             subset(genes_of_interest_annot$ensembl_transcript_id, 
                                    genes_of_interest_annot$refseq_mrna == i))
}

result_mirdb$ensembl_transcript_id <- ensembl_transcript_id

# Queremos comparar los resultados obtenidos experimentalmente, con scanMiR y con miRDB a nivel de transcrito. Vamos a trabajar con el *ENSEMBL transcript ID*.

# Primero, cambiamos los nombres de los miRNA en los resultados de scanMiR (ya que
# están en el formato de miRBase, mientras que en los resultados experimentales y
# de miRDB están diferentes):
levels(agg_matches_ORF_repr$miRNA) <- c("miR-106b", "miR-17", "miR-20a", "miR-27b", 
                                        "miR-30c")

# Filtramos los transcritos objetivo obtenidos por scanMiR que presentan mínimo 
# un sitio objetivo del tipo 8mer o 7mer en la región 3'UTR:
agg_matches_ORF_repr_3UTR <- subset(agg_matches_ORF_repr, agg_matches_ORF_repr[,4]!=0 |
                                      agg_matches_ORF_repr[,5]!=0)

# Filtramos los transcritos objetivo obtenidos experimentalmente donde se ha 
# encontrado un sitio objetivo del tipo 8mer, 7mer o 9mer en la región 3'UTR (vamos
# a incluir los sitios objetivo del tipo 9mer):
clash_seed_sub_3UTR <- subset(clash_seed_sub, !is.na(clash_seed_sub$X3.UTR) & 
                               (clash_seed_sub$seed_type == "8-mer" | 
                                  clash_seed_sub$seed_type == "7-mer" | 
                                  clash_seed_sub$seed_type == "9-mer"))

# Utilizamos la función semi_join(x, y, by=c()) que nos permite seleccionar las filas
# en x que tienen una o más coincidencias en y. Así, obtendremos el número de 
# objetivos de x que también han sido obtenidos en y.

# Un 0.80% de los resultados obtenidos con scanMiR han sido obtenidos experimentalmente:
scanmir_clash <- semi_join(agg_matches_ORF_repr, clash_seed_sub, by=c("miRNA"="mirna", 
                                                                    "transcript"="ENST_ID"))
dim(scanmir_clash)[1]/dim(agg_matches_ORF_repr)[1]*100

# Un 3.05% de los resultados obtenidos experimentalmente han sido obtenidos con 
# scanMiR:
clash_scanmir <- semi_join(clash_seed_sub, agg_matches_ORF_repr, 
                          by=c("mirna"="miRNA", "ENST_ID"="transcript"))
dim(clash_scanmir)[1]/dim(clash_seed_sub)[1]*100

# Un 15.15% de los resultados obtenidos con scanMiR (resultados filtrados para 
# quedarnos con aquellos transcritos que presentan mínimo un sitio objetivo del 
# tipo 8mer o 7mer en la región 3'UTR) han sido obtenidos con miRDB:
scanmir_mirdb <- semi_join(agg_matches_ORF_repr_3UTR, result_mirdb,
                           by=c("miRNA"="mirna", "transcript"="ensembl_transcript_id"))
dim(scanmir_mirdb)[1]/dim(agg_matches_ORF_repr_3UTR)[1]*100

# Un 20.05% de los resultados obtenidos con miRDB han sido obtenidos con scanMiR 
# (resultados filtrados para quedarnos con aquellos que presentan mínimo un sitio
# objetivo del tipo 8mer o 7mer en la región 3'UTR):
mirdb_scanmir <- semi_join(result_mirdb, agg_matches_ORF_repr_3UTR,
                           by=c("mirna"="miRNA", "ensembl_transcript_id"="transcript"))
dim(mirdb_scanmir)[1]/dim(result_mirdb)[1]*100

# Un 0.66% de los resultados obtenidos con miRDB han sido obtenidos experimentalmente
# (resultados filtrados para quedarnos con aquellos transcritos donde se ha 
# encontrado un sitio objetivo del tipo 8mer, 7mer o 9mer en la región 3'UTR):
mirdb_clash <- semi_join(result_mirdb, clash_seed_sub_3UTR, by=c("mirna",
                                                               "ensembl_transcript_id"="ENST_ID"))
dim(mirdb_clash)[1]/dim(result_mirdb)[1]*100

# Un 8.93% de los resultados obtenidos experimentalmente (resultados filtrados 
# para quedarnos con aquellos transcritos donde se ha encontrado un sitio objetivo 
# del tipo 8mer, 7mer o 9mer en la región 3'UTR) han sido obtenidos con miRDB:
clash_mirdb <- semi_join(clash_seed_sub_3UTR, result_mirdb, 
                        by=c("mirna", "ENST_ID"="ensembl_transcript_id"))
dim(clash_mirdb)[1]/dim(clash_seed_sub_3UTR)[1]*100


# Creamos la función percentage_results() para obtener los porcentajes mostrados
# anteriormente. Los argumentos de la función son dos data.frames, x e y, que
# corresponden con la x e y de la función semi_join(), los cuales únicamente 
# tienen 2 columnas: la primera corresponde con los mirna y la segunda con los 
# ENSEMBL transcript ID:
percentage_results <- function(df_x_mirna_enst_id, df_y_mirna_enst_id){
  
  names(df_x_mirna_enst_id)[1] <- "mirna"
  names(df_y_mirna_enst_id)[1] <- "mirna"
  names(df_x_mirna_enst_id)[2] <- "enst_id"
  names(df_y_mirna_enst_id)[2] <- "enst_id"
  
  result_inner_join <- semi_join(df_x_mirna_enst_id, df_y_mirna_enst_id, 
                                 by=c("mirna", "enst_id"))
  
  return((dim(result_inner_join)[1]/dim(df_x_mirna_enst_id)[1])*100)
}

# Por ejemplo:
percentage_results(clash_seed_sub_3UTR[,c(2,5)], result_mirdb[,c(2,7)])

