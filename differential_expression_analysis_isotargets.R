library(dplyr)
library(digest)
library(knitr)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(tibble)
library(stringr)

## Third specific objective: 

# Load scanmiR_isomirs.RData file:
scanmiR_isomirs <- load(file.path("C:/Users/marta/Desktop",
                                  "/masterBioinformaticaBioestadistica/TFM/scanmiR_isomirs.RData"))

# Show the first rows of isomiR.df.all:
head(isomiR.df.all)

# Check range of sequence lengths for scanMiR:
range(apply(isomiR.df.all["Read"], 2, nchar))

# Set empty values to 0 and transform to numeric vector the columns 6:22:
isomiR.df.all[isomiR.df.all == ""] <- 0

for(i in 6:22){
  isomiR.df.all[,i] <- as.numeric(isomiR.df.all[,i])
}

# Keep rows with at least 10 reads in 4 samples:
isomir_subset <- isomiR.df.all[rowSums(isomiR.df.all[6:22] >= 10) >= 4,]

# Set directory where the script is:
setwd('C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM')

# Create a .csv file for scanMiR use:
write.csv(isomir_subset[,c("UID.x", "Read")], file="azoos_isomir_subset.csv",
          row.names=FALSE)

# Run scanMiR script:
system('Rscript script_scanmir.R
       C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM
       azoos_isomir_subset.csv 200 result_scanmir_azoos C:/Users/marta/anaconda3/python.exe')

result_scanmir_azoos <- read.csv(file="result_scanmir_azoos.csv", header=TRUE, sep=',')

# Check if transcript and mirna are character vectors:
class(result_scanmir_azoos$transcript)
class(result_scanmir_azoos$mirna)

# Keep miRNA-mRNA interactions with repression < -1. Group rows by mirna column and
# assign a hash code to each isomiR (same hash code to isomiR with the same 200 targets).
df_filtered <- as.data.frame(result_scanmir_azoos %>%
                               filter(repression < -1))
dim(df_filtered)
length(unique(df_filtered$mirna))

df_filtered_mirna_transcript_code <- as.data.frame(df_filtered %>%
                                                     group_by(mirna) %>%
                                                     summarize(transcript_code = digest::digest(sort(transcript))) %>%
                                                     ungroup())

# Number of unique hash code:
length(unique(df_filtered_mirna_transcript_code$transcript_code))

# Repeat with repression < -1.5:
df_filtered_rep1_5 <- as.data.frame(result_scanmir_azoos %>%
                                      filter(repression < -1.5))
dim(df_filtered_rep1_5)
length(unique(df_filtered_rep1_5$mirna))

df_filtered_mirna_transcript_code_rep1_5 <- as.data.frame(df_filtered_rep1_5 %>%
                                                            group_by(mirna) %>%
                                                            summarize(transcript_code = digest::digest(sort(transcript))) %>%
                                                            ungroup())

# Number of unique hash code:
length(unique(df_filtered_mirna_transcript_code_rep1_5$transcript_code))

# Repeat with repression < -2:
df_filtered_rep2 <- as.data.frame(result_scanmir_azoos %>%
                                    filter(repression < -2))
dim(df_filtered_rep2)
length(unique(df_filtered_rep2$mirna))

df_filtered_mirna_transcript_code_rep2 <- as.data.frame(df_filtered_rep2 %>%
                                                          group_by(mirna) %>%
                                                          summarize(transcript_code = digest::digest(sort(transcript))) %>%
                                                          ungroup())

# Number of unique hash code:
length(unique(df_filtered_mirna_transcript_code_rep2$transcript_code))

# Create a table with repression value, number of selected isomiR and number of isoTargets generated:
kable(data.frame("Repression"=c(-1, -1.5, -2), "miRNA_isomiR_number"=c(2871, 2705, 2080),
                 "isoTargets_number"=c(2309, 1877, 927)))

# Keep rows in isomir_subset when isomiR is in df_filtered_mirna_transcript_code_rep1_5. 
# The column with hash code is added to the new data.frame:
isomir_subset_rep1_5 <- inner_join(isomir_subset, 
                                   df_filtered_mirna_transcript_code_rep1_5,
                                   by=c("UID.x" = "mirna"))


## Fourth specific objective: 

# Group the isomiR counts belonging to the same isoTargets group in each sample:
isotargets_count <- as.data.frame(isomir_subset_rep1_5 %>%
                                    group_by(transcript_code) %>%
                                    summarize_at(vars(6:23), sum) %>%
                                    ungroup())

# Create a data.frame with information for the samples (condition and subcondition):
coldata <- data.frame(condition = c(rep("AO", 4), rep("AS", 9), rep("DS", 4)),
                      subcondition = c(rep("AO.N", 4), rep("AS.REQneg", 4), 
                                       rep("AS.REQpos", 5), rep("DS", 4)))

# Transform condition and subcondition in factor:
coldata$condition <- factor(coldata$condition)
coldata$subcondition <- factor(coldata$subcondition)

# Each row corresponds to a sample:
row.names(coldata) <- c("AO.N_1", "AO.N_2", "AO.N_3", "AO.N_4", "AS.REQneg_1", 
                        "AS.REQneg_2", "AS.REQneg_3", "AS.REQneg_4", "AS.REQpos_1",
                        "AS.REQpos_2", "AS.REQpos_3", "AS.REQpos_4", "AS.REQpos_5",
                        "DS_1", "DS_2", "DS_3", "DS_4")

# isoTargets (transcript_code) must be the row names:
row.names(isotargets_count) <- isotargets_count$transcript_code

# Remove rowCountsSum column:
isotargets_count <- isotargets_count[,2:18]

# Show the first rows of isotargets_count:
head(isotargets_count)

# Create the DESeqDataSet objects:
dds_c <- DESeqDataSetFromMatrix(countData=isotargets_count, colData=coldata,
                                design=~condition)

dds_sc <- DESeqDataSetFromMatrix(countData=isotargets_count, colData=coldata,
                                 design=~subcondition)

# Do differential expression analysis at isoTargets level:
dds_c <- DESeq(dds_c)
dds_sc <- DESeq(dds_sc)

# Show results when comparing AS and AO groups:
res_AS_AO <- results(dds_c, contrast=c("condition", "AS", "AO"), lfcThreshold=0.26,
                     alpha=0.05, altHypothesis="greaterAbs")
summary(res_AS_AO)

# Show results when comparing AS.REQneg and AS.REQpos groups:
res_ASREQneg_pos <- results(dds_sc, contrast=c("subcondition", "AS.REQneg", 
                                               "AS.REQpos"), lfcThreshold=0.26, alpha=0.05, 
                            altHypothesis="greaterAbs")
summary(res_ASREQneg_pos)

# Volcano plot for AS and AO conditions:
EnhancedVolcano(res_AS_AO, lab="", x="log2FoldChange", y="padj", pCutoff=0.05, 
                FCcutoff=0.26, title="isoTargets_condition_AS_vs_AO", ylim=c(0,7), 
                xlim=c(-2,7), legendLabSize=12, legendIconSize=5, 
                legendPosition="right")

# Volcano plot for AS.REQneg and AS.REQpos subconditions:
EnhancedVolcano(res_ASREQneg_pos, lab="", x="log2FoldChange", y="padj", 
                pCutoff=0.05, FCcutoff=0.26, title="isoTargets_subcondition_ASREQneg_pos", 
                ylim=c(0,5), xlim=c(-5, 6), legendLabSize=12, legendIconSize=5, 
                legendPosition = "right")

# Keep statistically significant |LFC|>0.26 isoTargets:
isotargets_sig_c <- subset(res_AS_AO, padj<0.05)
isotargets_sig_sc <- subset(res_ASREQneg_pos, padj<0.05)

# Number of canonical miRNA and other isomiR variants in differentially expressed isoTargets (AS and AO conditions):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Number of canonical miRNA and other isomiR variants in differentially expressed isoTargets (AS.REQneg and AS.REQpos subconditions):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Calculate shrunken log2 fold changes (method apeglm):
res_AS_AO_LFC <- lfcShrink(dds_c, coef="condition_AS_vs_AO", res=res_AS_AO, 
                           type="apeglm")

# Show results when comparing AS and AO groups with shrunken log2 fold changes 
# (|LFC|>0.26):
summary(res_AS_AO_LFC)

resultsNames(dds_sc)

# Because subcondition_AS.REQneg_vs_AS.REQpos is not in resultsNames(dds_sc):
resultsNames(dds_sc)

# Select AS.REQpos as reference level. Use DESeq() after changing the reference level:
dds_sc$subcondition <- relevel(dds_sc$subcondition, ref = "AS.REQpos")
dds_sc <- DESeq(dds_sc)
resultsNames(dds_sc)
res_ASREQneg_pos_LFC <- lfcShrink(dds_sc, coef="subcondition_AS.REQneg_vs_AS.REQpos", 
                                  res=res_ASREQneg_pos, type="apeglm")
summary(res_ASREQneg_pos_LFC)

# Check if p-values and adjusted p-values change when working with shrunken log2 fold changes:
head(subset(res_ASREQneg_pos, padj<0.05 & abs(log2FoldChange)>0.26),3)
head(subset(res_ASREQneg_pos_LFC, padj<0.05 & abs(log2FoldChange)>0.26),3)

# MA-plots with log2 fold changes or shrunken log2 fold changes on the y-axis (AS and AO conditions): 
par(mfrow=c(1,2))

plotMA(res_AS_AO, ylab="log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

plotMA(res_AS_AO_LFC, ylab="shrunken log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

# Add a title shared by both plots:
mtext(text="condition_AS_vs_AO", side=3, line=-2, outer=TRUE)

# MA-plots with log2 fold changes or shrunken log2 fold changes on the y-axis (AS.REQneg and AS.REQpos subconditions): 
par(mfrow=c(1,2))

plotMA(res_ASREQneg_pos, ylab="log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

plotMA(res_ASREQneg_pos_LFC, ylab="shrunken log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

# Add a title shared by both plots:
mtext(text="subcondition_AS.REQneg_vs_AS.REQpos", side=3, line=-2, outer=TRUE)

# Do regularized logarithm transformation:
vsd_c <- vst(dds_c, blind=FALSE)
vsd_sc <- vst(dds_sc, blind=FALSE)

# Plot standard deviation of log2 scaled counts for each isoTargets against the mean:
meanSd_log2_scaled_count <- vsn::meanSdPlot(assay(normTransform(dds_c)), ranks=FALSE, plot=FALSE)
meanSd_log2_scaled_count$gg + ggtitle("log2 scaled counts")

# Plot standard deviation of regularized log transformation counts for each isoTargets against the mean:
meanSd_reg_log_transf_counts <- vsn::meanSdPlot(assay(vsd_c), ranks=FALSE, plot=FALSE)
meanSd_reg_log_transf_counts$gg + ggtitle("Regularized log transformation counts")

# PCA plot using isoTargets with a statistically significant |LFC|>0.26 (AO and AS conditions):
PCA_sig_c <- plotPCA(vsd_c[rownames(isotargets_sig_c), c(1:13)], intgroup="condition")

PCA_sig_c + geom_text(aes(label = rownames(coldata)[1:13]), size=3, 
                      angle=15, color="black") + coord_fixed(ylim=c(-15, 17), xlim=c(-31, 27))

# PCA plot using isoTargets with a statistically significant |LFC|>0.26 (AS.REQneg and AS.REQpos subconditions):
PCA_sig_sc <- plotPCA(vsd_sc[rownames(isotargets_sig_sc), c(5:13)], 
                      intgroup="subcondition")

PCA_sig_sc + geom_text(aes(label = rownames(coldata)[5:13]), size=3, 
                       angle=70, color="black") + coord_fixed(ylim=c(-8, 7), xlim=c(-11, 11))

# Create two data.frame with names of samples to be compared:
df_c <- as.data.frame(colData(dds_c)[1:13,][1])
df_sc <- as.data.frame(colData(dds_sc)[5:13,][2])

# Heatmap for the AS and AO conditions: 
pheatmap(assay(vsd_c)[rownames(isotargets_sig_c), c(1:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_c, scale="row", border_color=NA)

# Heatmap for the AS.REQneg and AS.REQpos subconditions: 
pheatmap(assay(vsd_sc)[rownames(isotargets_sig_sc), c(5:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_sc, scale="row", border_color=NA)

# Similar procedure for differential expression analysis at isomiR level:
isomir_count <- isomir_subset[,6:22]

# Create the DESeqDataSet objects:
dds_isomir_c <- DESeqDataSetFromMatrix(countData=isomir_count, colData=coldata,
                                       design=~condition)

dds_isomir_sc <- DESeqDataSetFromMatrix(countData=isomir_count, colData=coldata,
                                        design=~subcondition)

# Do differential expression analysis at isomiR level:
dds_isomir_c <- DESeq(dds_isomir_c)
dds_isomir_sc <- DESeq(dds_isomir_sc)

# Show results when comparing AS and AO groups:
res_isomir_AS_AO <- results(dds_isomir_c, contrast=c("condition", "AS", "AO"), 
                            lfcThreshold=0.26, alpha=0.05, altHypothesis="greaterAbs")
summary(res_isomir_AS_AO)

# Show results when comparing AS.REQneg and AS.REQpos groups:
res_isomir_ASREQneg_pos <- results(dds_isomir_sc, contrast=c("subcondition", "AS.REQneg", 
                                                             "AS.REQpos"), lfcThreshold=0.26, alpha=0.05, 
                                   altHypothesis="greaterAbs")
summary(res_isomir_ASREQneg_pos)

# Volcano plot for AS and AO conditions:
EnhancedVolcano(res_isomir_AS_AO, lab="", x="log2FoldChange", y="padj", pCutoff=0.05, 
                FCcutoff=0.26, title="isomiR_condition_AS_vs_AO", ylim=c(0,7), 
                xlim=c(-2,7), legendLabSize=12, legendIconSize=5, 
                legendPosition="right")

# Volcano plot for AS.REQneg and AS.REQpos subconditions:
EnhancedVolcano(res_isomir_ASREQneg_pos, lab="", x="log2FoldChange", y="padj", 
                pCutoff=0.05, FCcutoff=0.26, title="isomiR_subcondition_ASREQneg_pos", 
                ylim=c(0,5), xlim=c(-5, 6), legendLabSize=12, legendIconSize=5, 
                legendPosition = "right")

# Keep statistically significant |LFC|>0.26 isomiR:
isomir_sig_c <- subset(res_isomir_AS_AO, padj<0.05)
isomir_sig_sc <- subset(res_isomir_ASREQneg_pos, padj<0.05)

# Number of canonical miRNA and other isomiR variants differentially expressed (AS and AO conditions):
as.data.frame(isomir_sig_c) %>%
  mutate(mirna_canonico = row.names(isomir_sig_c) %in% subset(isomir_subset$UID.x, isomir_subset$variant == "NA")) %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Number of canonical miRNA and other isomiR variants differentially expressed (AS.REQneg and AS.REQpos subconditions):
as.data.frame(isomir_sig_sc) %>%
  mutate(mirna_canonico = row.names(isomir_sig_sc) %in% subset(isomir_subset$UID.x, isomir_subset$variant == "NA")) %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Do regularized logarithm transformation:
vsd_isomir_c <- vst(dds_isomir_c, blind=FALSE)
vsd_isomir_sc <- vst(dds_isomir_sc, blind=FALSE)

# PCA plot using isomiR with a statistically significant |LFC|>0.26 (AO and AS conditions):
PCA_sig_isomir_c <- plotPCA(vsd_isomir_c[rownames(isomir_sig_c), c(1:13)], 
                      intgroup="condition")

PCA_sig_isomir_c + geom_text(aes(label = rownames(coldata)[1:13]), size=3, 
                       angle=15, color="black") + coord_fixed(ylim=c(-18, 17), xlim=c(-42, 34))

# PCA plot using isomiR with a statistically significant |LFC|>0.26 (AS.REQneg and AS.REQpos subconditions):
PCA_sig_isomir_sc <- plotPCA(vsd_isomir_sc[rownames(isomir_sig_sc), c(5:13)], 
                      intgroup="subcondition")

PCA_sig_isomir_sc + geom_text(aes(label = rownames(coldata)[5:13]), size=3, 
                       angle=70, color="black") + coord_fixed(ylim=c(-8, 7), xlim=c(-11, 11))

# Create two data.frame with names of samples to be compared:
df_isomir_c <- as.data.frame(colData(dds_isomir_c)[1:13,][1])
df_isomir_sc <- as.data.frame(colData(dds_isomir_sc)[5:13,][2])

# Heatmap for the AS and AO conditions: 
pheatmap(assay(vsd_isomir_c)[rownames(isomir_sig_c), c(1:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_isomir_c, scale="row")

# Heatmap for the AS.REQneg and AS.REQpos subconditions: 
pheatmap(assay(vsd_isomir_sc)[rownames(isomir_sig_sc), c(5:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_isomir_sc, scale="row", border_color=NA)

# Similar procedure for differential expression analysis at miRNA level:
mirna_count <- as.data.frame(isomir_subset %>%
                               group_by(parent) %>%
                               summarize_at(vars(5:21), sum) %>%
                               ungroup() %>%
                               column_to_rownames("parent"))

# Create the DESeqDataSet objects:
dds_mirna_c <- DESeqDataSetFromMatrix(countData=mirna_count, colData=coldata,
                                      design=~condition)

dds_mirna_sc <- DESeqDataSetFromMatrix(countData=mirna_count, colData=coldata,
                                       design=~subcondition)

# Do differential expression analysis at miRNA level:
dds_mirna_c <- DESeq(dds_mirna_c)
dds_mirna_sc <- DESeq(dds_mirna_sc)

# Show results when comparing AS and AO groups:
res_mirna_AS_AO <- results(dds_mirna_c, contrast=c("condition", "AS", "AO"), 
                           lfcThreshold=0.26, alpha=0.05, altHypothesis="greaterAbs")
summary(res_mirna_AS_AO)

# Show results when comparing AS.REQneg and AS.REQpos groups:
res_mirna_ASREQneg_pos <- results(dds_mirna_sc, contrast=c("subcondition", "AS.REQneg", 
                                                           "AS.REQpos"), lfcThreshold=0.26, alpha=0.05, 
                                  altHypothesis="greaterAbs")
summary(res_mirna_ASREQneg_pos)

# Volcano plot for AS and AO conditions:
EnhancedVolcano(res_mirna_AS_AO, lab="", x="log2FoldChange", y="padj", pCutoff=0.05, 
                FCcutoff=0.26, title="miRNA_condition_AS_vs_AO", ylim=c(0,7), 
                xlim=c(-2,7), legendLabSize=12, legendIconSize=5, 
                legendPosition="right")

# Volcano plot for AS.REQneg and AS.REQpos subconditions:
EnhancedVolcano(res_mirna_ASREQneg_pos, lab="", x="log2FoldChange", y="padj", 
                pCutoff=0.05, FCcutoff=0.26, title="miRNA_subcondition_ASREQneg_pos", 
                ylim=c(0,5), xlim=c(-5, 6), legendLabSize=12, legendIconSize=5, 
                legendPosition = "right")

# Keep statistically significant |LFC|>0.26 miRNA:
mirna_sig_c <- subset(res_mirna_AS_AO, padj<0.05)
mirna_sig_sc <- subset(res_mirna_ASREQneg_pos, padj<0.05)

# Number of canonical miRNA and other isomiR variants in differentially expressed miRNA (AS and AO conditions):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(parent %in% row.names(mirna_sig_c)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Number of canonical miRNA and other isomiR variants in differentially expressed miRNA (AS.REQneg and AS.REQpos subconditions):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(parent %in% row.names(mirna_sig_sc)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Do regularized logarithm transformation:
vsd_mirna_c <- varianceStabilizingTransformation(dds_mirna_c, blind=FALSE)
vsd_mirna_sc <- varianceStabilizingTransformation(dds_mirna_sc, blind=FALSE)

# PCA plot using miRNA with a statistically significant |LFC|>0.26 (AO and AS conditions):
PCA_sig_mirna_c <- plotPCA(vsd_mirna_c[rownames(mirna_sig_c), c(1:13)], 
                            intgroup="condition")

PCA_sig_mirna_c + geom_text(aes(label = rownames(coldata)[1:13]), size=3, 
                             angle=15, color="black") + coord_fixed(ylim=c(-5, 9), xlim=c(-19, 15))

# PCA plot using miRNA with a statistically significant |LFC|>0.26 (AS.REQneg and AS.REQpos subconditions):
PCA_sig_mirna_sc <- plotPCA(vsd_mirna_sc[rownames(mirna_sig_sc), c(5:13)], 
                             intgroup="subcondition")

PCA_sig_mirna_sc + geom_text(aes(label = rownames(coldata)[5:13]), size=3, 
                              angle=70, color="black") + coord_fixed(ylim=c(-7, 6), xlim=c(-8, 7))

# Create two data.frame with names of samples to be compared:
df_mirna_c <- as.data.frame(colData(dds_mirna_c)[1:13,][1])
df_mirna_sc <- as.data.frame(colData(dds_mirna_sc)[5:13,][2])

# Heatmap for the AS and AO conditions: 
pheatmap(assay(vsd_mirna_c)[rownames(mirna_sig_c), c(1:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_mirna_c, scale="row", border_color=NA)

# Heatmap for the AS.REQneg and AS.REQpos subconditions: 
pheatmap(assay(vsd_mirna_sc)[rownames(mirna_sig_sc), c(5:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_mirna_sc, scale="row", border_color=NA)

# Study of isoTargets. 

# Select isomiR of isoTargets with more than 1 isomiR:
isomir_arrange <- as.data.frame(isomir_subset_rep1_5[,c(2:4,24)] %>%
                                  arrange(transcript_code))

isomir_isotargets_multiple <- as.data.frame(isomir_arrange %>%
                                              group_by(transcript_code) %>%
                                              mutate(isotargets_multiple=n()) %>%
                                              filter(isotargets_multiple>1))

# Number of isoTargets with more than 1 isomiR:
dim(as.data.frame(isomir_isotargets_multiple %>%
                    group_by(transcript_code) %>%
                    summarize(count_distinct_parent=n_distinct(parent))))[1]

# Number of canonical miRNA and other isomiR variants of isoTargets with more than one isomiR:
dim(isomir_isotargets_multiple)[1]

# Percentage of 5p isomiR in isoTargets with more than 1 isomiR:
(dim(as.data.frame(isomir_isotargets_multiple %>% filter(str_detect(variant, "5p"))))[1]/
dim(as.data.frame(isomir_arrange %>% filter(str_detect(variant, "5p"))))[1])*100

# Check if the isoTargets of these 5p isomiR also contain the canonical miRNA. 
# Canonical_isotargets is the sum of all variants in the transcript_code-parent group when 
# canonical miRNA is in the transcript_code-parent group. All 5p variants and their 
# canonical miRNA are in distinct isoTargets:
semi_join(isomir_arrange, (as.data.frame(isomir_isotargets_multiple %>%
                                           filter(str_detect(variant, "5p")))), by="transcript_code") %>%
  group_by(transcript_code, parent) %>%
  summarize(canonical_isotargets = sum(str_detect(variant, "5p") & 
                                         any(variant == "NA"))) %>%
  ungroup() %>%
  summarize(sum(canonical_isotargets))

# Number of canonical miRNA corresponding to 5p variants in isomir_arrange:
isomir_arrange %>%
  filter(str_detect(variant, "5p")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Percentage of 3p isomiR in isoTargets with more than 1 isomiR:
(dim(as.data.frame(isomir_isotargets_multiple %>% filter(str_detect(variant, "3p"))))[1]/
dim(as.data.frame(isomir_arrange %>% filter(str_detect(variant, "3p"))))[1])*100

# Check if the isoTargets of these 3p isomiR also contain the canonical miRNA. 
# Canonical_isotargets is the sum of all variants in the transcript_code-parent group when 
# canonical miRNA is in the transcript_code-parent group:
semi_join(isomir_arrange, (as.data.frame(isomir_isotargets_multiple %>%
                                           filter(str_detect(variant, "3p")))), by="transcript_code") %>%
  group_by(transcript_code, parent) %>%
  summarize(canonical_isotargets = sum(str_detect(variant, "3p") & 
                                         any(variant == "NA"))) %>%  
  ungroup() %>%
  summarize(sum(canonical_isotargets))

# Number of canonical miRNA corresponding to 3p variants in isomir_arrange:
isomir_arrange %>%
  filter(str_detect(variant, "3p")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Percentage of iso_snv_seed in isoTargets with more than 1 isomiR:
(dim(as.data.frame(isomir_isotargets_multiple %>% filter(str_detect(variant, "snv_seed"))))[1]/
dim(as.data.frame(isomir_arrange %>% filter(str_detect(variant, "snv_seed"))))[1])*100

# Number of canonical miRNA corresponding to iso_snv_seed in isomir_arrange:
isomir_arrange %>%
  filter(str_detect(variant, "iso_snv_seed")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Percentage of iso_snv (without iso_snv_seed) in isoTargets with more than 1 isomiR:
(dim(as.data.frame(isomir_isotargets_multiple %>%
                    filter(str_detect(variant, "iso_snv")) %>%
                    filter(!str_detect(variant, "iso_snv_seed"))))[1]/
dim(as.data.frame(isomir_arrange %>%
                    filter(str_detect(variant, "iso_snv")) %>%
                    filter(!str_detect(variant, "iso_snv_seed"))))[1])*100

# Check if the isoTargets of these iso_snv variants also contain the canonical miRNA. 
# Canonical_isotargets is the sum of all variants in the transcript_code-parent group when 
# canonical miRNA is in the transcript_code-parent group:
semi_join(isomir_arrange, (as.data.frame(isomir_arrange %>%
                                           filter(str_detect(variant, "iso_snv")) %>%
                                           filter(!str_detect(variant, "iso_snv_seed")))), by="transcript_code") %>%
  group_by(transcript_code, parent) %>%
  summarize(canonical_isotargets = sum(str_detect(variant, "iso_snv") & 
                                         any(variant == "NA"))) %>%
  ungroup() %>%
  summarize(sum(canonical_isotargets))

# Number of canonical miRNA corresponding to iso_snv (without iso_snv_seed) in isomir_arrange:
isomir_arrange %>%
  filter(str_detect(variant, "iso_snv")) %>%
  filter(!str_detect(variant, "iso_snv_seed")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Plot the results in pie charts:
par(mfrow=c(2,2), mar=c(2.5,2.5,2.5,2.5))

colors <- c("#99FF99", "#FFD699")

five <- 100
cat_five <- "diff_isoTargets"

pie(five, labels=cat_five, col=colors, cex=0.75, main="5p")

three <- c((1651-518)/1651*100, 100-((1651-518)/1651*100))
cat_three <- c("diff_isoTargets", "same_isoTargets")

pie(three, labels=cat_three, col=colors, cex=0.75, main="3p")

snv_seed <- 100
cat_snv_seed <- "diff_isoTargets"

pie(snv_seed, labels=cat_snv_seed, col=colors, cex=0.75, main="snv_seed")

snv_noseed <- c((303-12)/303*100, 100-((303-12)/303*100))
cat_snv_noseed <- c("diff_isoTargets", "same_isoTargets")

pie(snv_noseed, labels=cat_snv_noseed, col=colors, cex=0.75, main="snv_noseed")

# Study of the results obtained with the 3 different approaches. 

# Percentage of isomiR of isoTargets differentially expressed in AS and AO conditions
# that match those obtained using the isomiR level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
  summarize(count(UID.x %in% row.names(isomir_sig_c))))[,1]/
dim(isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_c)))[1])*100

# Percentage of isomiR differentially expressed at isomiR level in AS and AO conditions
# that match those obtained using the isoTargets level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
    summarize(count(UID.x %in% row.names(isomir_sig_c))))[,1]/
dim(isomir_sig_c)[1])*100

# Percentage of isomiR of isoTargets differentially expressed in AS.REQneg and
# AS.REQpos subconditions that match those obtained using the isomiR level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
  summarize(count(UID.x %in% row.names(isomir_sig_sc))))[1]/
dim(isomir_subset_rep1_5[,c(2:4,24)] %>% 
    filter(transcript_code %in% row.names(isotargets_sig_sc)))[1])*100

# Percentage of isomiR differentially expressed at isomiR level in AS.REQneg
# and AS.REQpos subconditions that match those obtained using the isoTargets level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
    summarize(count(UID.x %in% row.names(isomir_sig_sc))))[,1]/
dim(isomir_sig_sc)[1])*100

# Percentage of isomiR of isoTargets differentially expressed in AS and AO 
# conditions that match those obtained using the miRNA level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
  summarize(count(parent %in% row.names(mirna_sig_c))))[,1]/
dim(isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_c)))[1])*100

# Percentage of isomiR differentially expressed at miRNA level in AS and AO conditions
# that match those obtained using the isoTargets level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
    summarize(count(parent %in% row.names(mirna_sig_c))))[,1]/
dim(isomir_subset %>% filter(parent %in% row.names(mirna_sig_c)))[1])*100

# Percentage of isomiR of isoTargets differentially expressed in AS.REQneg and 
# AS.REQpos subconditions that match those obtained using the miRNA level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
    summarize(count(parent %in% row.names(mirna_sig_sc))))[,1]/
dim(isomir_subset_rep1_5[,c(2:4,24)] %>% 
    filter(transcript_code %in% row.names(isotargets_sig_sc)))[1])*100

# Percentage of isomiR differentially expressed at miRNA level in AS.REQneg and
# AS.REQpos subconditions that match those obtained using the isoTargets level approach:
((isomir_subset_rep1_5[,c(2:4,24)] %>%
    filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
    summarize(count(parent %in% row.names(mirna_sig_sc))))[,1]/
dim(isomir_subset %>% filter(parent %in% row.names(mirna_sig_sc)))[1])*100

# Compare the p-value and adjusted p-value for each differentially expressed isoTargets
# with the mean p-value and adjusted p-value obtained for each isomiR.

# Keep differentially expressed isoTargets composed of more than 1 isomiR 
# and add the UID.x column (AS and AO conditions):
isotargets_sig_c_UID <- as.data.frame(isotargets_sig_c) %>%
  filter(row.names(isotargets_sig_c) %in% isomir_isotargets_multiple$transcript_code) %>%
  merge(isomir_isotargets_multiple %>% 
          select(transcript_code, UID.x),
        by.x = 0, by.y = "transcript_code", all.x = TRUE) %>%
  mutate(transcript_code=Row.names) %>%
  select(-Row.names)

# Create a data.frame with the columns transcript_code, pvalue and padj for differentially 
# expressed isoTargets composed of more than 1 isomiR (AS and AO conditions):
isotargets_multiple_sig_c <- as.data.frame(isotargets_sig_c) %>%
  filter(row.names(isotargets_sig_c) %in% isomir_isotargets_multiple$transcript_code) %>%
  rownames_to_column("transcript_code") %>%
  select(c(transcript_code, pvalue, padj))

# Create a data.frame with the columns transcript_code, pvalue_mean and padj_mean 
# with the mean p-values and adjusted p-values of those isomiR (isomiR level) 
# that are part of differentially expressed isoTargets composed of more than 1 isomiR
# (AS and AO conditions):
isomir_multiple_c <- as.data.frame(as.data.frame(res_isomir_AS_AO) %>%
                                     rownames_to_column("UID.x") %>%
                                     inner_join(isotargets_sig_c_UID, by = c("UID.x"="UID.x")) %>%
                                     group_by(transcript_code) %>%
                                     summarize(pvalue_mean = mean(pvalue.x, na.rm=TRUE), padj_mean=mean(padj.x, 
                                                                                                        na.rm=TRUE)))

# Compare the p-values with the mean p-values, and the adjusted p-values with the
# mean adjusted p-values for each isoTargets (AS and AO conditions):
isotargets_multiple_sig_c %>%
  inner_join(isomir_multiple_c, by="transcript_code") %>%
  mutate(pvalue_isotargets_low=pvalue<pvalue_mean, padj_isotargets_low=padj<padj_mean)

# Keep differentially expressed isoTargets composed of more than 1 isomiR 
# and add the UID.x column (AS.REQneg and AS.REQpos subconditions):
isotargets_sig_sc_UID <- as.data.frame(isotargets_sig_sc) %>%
  filter(row.names(isotargets_sig_sc) %in% isomir_isotargets_multiple$transcript_code) %>%
  merge(isomir_isotargets_multiple %>% 
          select(transcript_code, UID.x),
        by.x = 0, by.y = "transcript_code", all.x = TRUE) %>%
  mutate(transcript_code=Row.names) %>%
  select(-Row.names)

# Create a data.frame with the columns transcript_code, pvalue and padj for differentially 
# expressed isoTargets composed of more than 1 isomiR (AS.REQneg and AS.REQpos subconditions):
isotargets_multiple_sig_sc <- as.data.frame(isotargets_sig_sc) %>%
  filter(row.names(isotargets_sig_sc) %in% isomir_isotargets_multiple$transcript_code) %>%
  rownames_to_column("transcript_code") %>%
  select(c(transcript_code, pvalue, padj))

# Create a data.frame with the columns transcript_code, pvalue_mean and padj_mean 
# with the mean p-values and adjusted p-values of those isomiR (isomiR level) 
# that are part of differentially expressed isoTargets composed of more than 1 isomiR
# (AS.REQneg and AS.REQpos subconditions):
isomir_multiple_sc <- as.data.frame(as.data.frame(res_isomir_ASREQneg_pos) %>%
                                     rownames_to_column("UID.x") %>%
                                     inner_join(isotargets_sig_sc_UID, by = c("UID.x"="UID.x")) %>%
                                     group_by(transcript_code) %>%
                                     summarize(pvalue_mean = mean(pvalue.x, na.rm=TRUE), padj_mean=mean(padj.x, 
                                                                                                        na.rm=TRUE)))

# Compare the p-values with the mean p-values, and the adjusted p-values with the
# mean adjusted p-values for each isoTargets (AS.REQneg and AS.REQpos subconditions):
isotargets_multiple_sig_sc %>%
  inner_join(isomir_multiple_sc, by="transcript_code") %>%
  mutate(pvalue_isotargets_low=pvalue<pvalue_mean, padj_isotargets_low=padj<padj_mean)

