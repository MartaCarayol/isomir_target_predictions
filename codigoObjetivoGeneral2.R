# OBJETIVO GENERAL 2

# Cargamos el archivo scanmiR_isomirs.RData, que contiene varios objetos:

scanmiR_isomirs <- load(file.path("C:/Users/marta/Desktop",
                                  "/masterBioinformaticaBioestadistica/TFM/scanmiR_isomirs.RData"))

# Mostramos los objetos que contiene:
scanmiR_isomirs

# Mostramos los primeros elementos de isomiR.df.all:
head(isomiR.df.all)

# Importante para aplicar scanMiR: rango de longitudes de las secuencias:
range(apply(isomiR.df.all["Read"], 2, nchar))

# Las columnas 6:22 corresponden con vectores de tipo character. Los transformamos
# en vectores numéricos para filtrar:
isomiR.df.all[isomiR.df.all == ""] <- 0

for(i in 6:22){
  isomiR.df.all[,i] <- as.numeric(isomiR.df.all[,i])
}

# De esta forma nos quedamos con 2900 miRNA/isomiR:
isomir_subset <- isomiR.df.all[rowSums(isomiR.df.all[6:22] >= 10) >= 4,]
dim(isomir_subset)

# Establecemos el directorio de trabajo:
setwd('C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM')

# Creamos el archivo .csv para aplicar el script de scanMiR:
write.csv(isomir_subset[,c("UID.x", "Read")], file="azoos_isomir_subset.csv", 
          row.names=FALSE)

# Aplicamos el script de scanMiR:
system('Rscript script_scanmir.R 
       C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM 
       azoos_isomir_subset.csv 200 C:/Users/marta/anaconda3/python.exe')

# Una vez finalizado, cargamos el archivo .csv generado con los resultados:
result_scanmir_azoos <- read.csv(file="result_scanmir.csv", header=TRUE, sep=',')

# Mostramos los primeros elementos:
head(result_scanmir_azoos)

# Tenemos 200 transcritos objetivo para cada uno de los 2900 miRNA/isomiR (en total,
# 580000 filas):
dim(result_scanmir_azoos)

library(dplyr)
library(digest)

# Primero, comprobamos que transcript y mirna corresponden con vectores de tipo 
# character:
class(result_scanmir_azoos$transcript)
class(result_scanmir_azoos$mirna)

# Filtramos por repression = -1:
df_filtered <- as.data.frame(result_scanmir_azoos %>%
                               filter(repression < -1))

# Nos quedamos con 531077 de 580000 transcritos objetivo:
dim(df_filtered)

# Nos quedamos con 2871 de 2900 miRNA/isomiR
length(unique(df_filtered$mirna))

# Utilizamos el código anterior para asignar un hash a cada miRNA/isomiR 
# (transcript_code):
df_filtered_mirna_transcript_code <- as.data.frame(df_filtered %>%
                                                     group_by(mirna) %>%
                                                     summarize(transcript_code = digest::digest(sort(transcript))) %>%
                                                     ungroup())

# Hemos obtenido 2309 hash únicos:
length(unique(df_filtered_mirna_transcript_code$transcript_code))

# Repetimos el mismo procedimiento utilizando un valor de repression = -1.5:
df_filtered_rep1_5 <- as.data.frame(result_scanmir_azoos %>%
                                      filter(repression < -1.5))

# Nos quedamos con 337404 de 580000 transcritos objetivo:
dim(df_filtered_rep1_5)

# Nos quedamos con 2705 de 2900 miRNA/isomiR
length(unique(df_filtered_rep1_5$mirna))

# Utilizamos el código anterior para asignar un hash a cada miRNA/isomiR 
# (transcript_code):
df_filtered_mirna_transcript_code_rep1_5 <- as.data.frame(df_filtered_rep1_5 %>%
                                                            group_by(mirna) %>%
                                                            summarize(transcript_code = digest::digest(sort(transcript))) %>%
                                                            ungroup())

# Hemos obtenido 1877 hash únicos (isoTargets):
length(unique(df_filtered_mirna_transcript_code_rep1_5$transcript_code))

# Repetimos el mismo procedimiento utilizando un valor de repression = -2:
df_filtered_rep2 <- as.data.frame(result_scanmir_azoos %>%
                                    filter(repression < -2))

# Nos quedamos con 101983 de 580000 transcritos objetivo:
dim(df_filtered_rep2)

# Nos quedamos con 2080 de 2900 miRNA/isomiR
length(unique(df_filtered_rep2$mirna))

# Utilizamos el código anterior para asignar un hash a cada miRNA/isomiR 
# (transcript_code):
df_filtered_mirna_transcript_code_rep2 <- as.data.frame(df_filtered_rep2 %>%
                                                          group_by(mirna) %>%
                                                          summarize(transcript_code = digest::digest(sort(transcript))) %>%
                                                          ungroup())

# Hemos obtenido 927 hash únicos:
length(unique(df_filtered_mirna_transcript_code_rep2$transcript_code))

# Mostramos los primeros elementos:
head(df_filtered_mirna_transcript_code_rep2, 20)

# Para conocer el número de objetivos de cada miRNA/isomiR podemos utilizar:
df_count <- as.data.frame(result_scanmir_azoos %>%
                            filter(repression < -2) %>%
                            dplyr::count(mirna)) 

library(knitr)

# Representamos en una tabla los distintos filtros: valor de repression, número 
# de miRNA/isomiR seleccionados y número de isoTargets generados:
kable(data.frame("Repression"=c(-1, -1.5, -2), "miRNA_isomiR_number"=c(2871, 2705, 2080),
                 "isoTargets_number"=c(2309, 1877, 927)))

# Utilizamos inner_join() del paquete dplyr para seleccionar aquellas filas de 
# isomir_subset que contienen los contajes de miRNA/isomiR presentes en 
# df_filtered_mirna_transcript_code_rep1_5. Se añade al nuevo data.frame la columna
# transcript_code de df_filtered_mirna_transcript_code_rep1_5:
isomir_subset_rep1_5 <- inner_join(isomir_subset, 
                                   df_filtered_mirna_transcript_code_rep1_5,
                                   by=c("UID.x" = "mirna"))

# Para agrupar los contajes de los miRNA/isomiR que pertenecen al mismo grupo 
# isoTargets en cada muestra, agrupamos las filas utilizando la columna 
# transcript_code (group_by(transcript_code)) y utilizamos la función summarize_at()
# para aplicar la función sum() a las variables de la 6:23, que corresponden con 
# las distintas muestras y la columna rowCountsSum. Finalmente, utilizamos ungroup()
# para obtener el tibble final sin agrupar (utilizamos as.data.frame() para obtener
# un data.frame):
isotargets_count <- as.data.frame(isomir_subset_rep1_5 %>%
                                    group_by(transcript_code) %>%
                                    summarize_at(vars(6:23), sum) %>%
                                    ungroup())

# Mostramos los primeros elementos:
head(isotargets_count)

# Construimos el data.frame con información de las muestras (variable condition, 
# con 3 niveles: AO, AS o DS, y variable subcondition, con 4 niveles: AO, AS.REQneg,
# AS.REQpos y DS) (veremos por qué utilizamos estas dos variables a continuación):
coldata <- data.frame(condition = c(rep("AO", 4), rep("AS", 9), rep("DS", 4)),
                      subcondition = c(rep("AO.N", 4), rep("AS.REQneg", 4), 
                                       rep("AS.REQpos", 5), rep("DS", 4)))

# Transformamos las variables condition y subcondition en variables de tipo factor:
coldata$condition <- factor(coldata$condition)
coldata$subcondition <- factor(coldata$subcondition)

# Cada una de las filas corresponde con una muestra. Las muestras deben aparecer 
# en el orden en el que están en las columnas de isotargets_count:
row.names(coldata) <- c("AO.N_1", "AO.N_2", "AO.N_3", "AO.N_4", "AS.REQneg_1", 
                        "AS.REQneg_2", "AS.REQneg_3", "AS.REQneg_4", "AS.REQpos_1",
                        "AS.REQpos_2", "AS.REQpos_3", "AS.REQpos_4", "AS.REQpos_5",
                        "DS_1", "DS_2", "DS_3", "DS_4")

# Mostramos los primeros elementos:
head(coldata)

# En isotargets_count los isoTargets no pueden aparecer como una columna más
# (transcript_code), sino que deben corresponder con los nombres de las filas 
# (las columnas solo deben ser las distintas muestras, también eliminamos 
# rowCountsSum):
row.names(isotargets_count) <- isotargets_count$transcript_code

isotargets_count <- isotargets_count[,2:18]

# Mostramos los primeros elementos:
head(isotargets_count)

library(DESeq2)

# Construimos los objetos DESeqDataSet:
dds_c <- DESeqDataSetFromMatrix(countData=isotargets_count, colData=coldata,
                                design=~condition)

dds_sc <- DESeqDataSetFromMatrix(countData=isotargets_count, colData=coldata,
                                 design=~subcondition)

# Llevamos a cabo el análisis de expresión diferencial:
dds_c <- DESeq(dds_c)
dds_sc <- DESeq(dds_sc)

# Mostramos los resultados al comparar los grupos AS y AO:
res_AS_AO <- results(dds_c, contrast=c("condition", "AS", "AO"), lfcThreshold=0.26,
                     alpha=0.05, altHypothesis="greaterAbs")
summary(res_AS_AO)

# Mostramos los resultados al comparar los grupos AS.REQneg y AS.REQpos:
res_ASREQneg_pos <- results(dds_sc, contrast=c("subcondition", "AS.REQneg", 
                                               "AS.REQpos"), lfcThreshold=0.26, alpha=0.05, 
                            altHypothesis="greaterAbs")
summary(res_ASREQneg_pos)

library(EnhancedVolcano)

# Utilizamos la función EnhancedVolcano() con distintos argumentos para obtener
# el volcano plot de las condiciones o grupos AS y AO:
EnhancedVolcano(res_AS_AO, lab="", x="log2FoldChange", y="padj", pCutoff=0.05, 
                FCcutoff=0.26, title="isoTargets_condition_AS_vs_AO", ylim=c(0,7), 
                xlim=c(-2,7), legendLabSize=12, legendIconSize=5, 
                legendPosition="right")

# Repetimos el procedimiento anterior para obtener el volcano plot de las 
# subcondiciones AS.REQneg y AS.REQpos:
EnhancedVolcano(res_ASREQneg_pos, lab="", x="log2FoldChange", y="padj", 
                pCutoff=0.05, FCcutoff=0.26, title="isoTargets_subcondition_ASREQneg_pos", 
                ylim=c(0,5), xlim=c(-5, 6), legendLabSize=12, legendIconSize=5, 
                legendPosition = "right")

# Seleccionamos aquellos isoTargets con un |LFC|>0.26 estadísticamente significativos:
isotargets_sig_c <- subset(res_AS_AO, padj<0.05)
isotargets_sig_sc <- subset(res_ASREQneg_pos, padj<0.05)

# Los 112 isoTargets diferencialmente expresados en las condiciones AS y AO 
# corresponden con 200 miRNA/isomiR. Vemos cuántos son isomiR y cuántos miRNA:
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Repetimos el mismo procedimiento en el caso de las subcondiciones AS.REQneg y 
# AS.REQpos. Los 34 isoTargets diferencialmente expresados en estas subcondiciones
# corresponden con 58 miRNA/isomiR. Vemos cuántos son isomiR y cuántos miRNA:
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Utilizamos la función lfcShrink() para trabajar con shrunken log2 fold changes:
res_AS_AO_LFC <- lfcShrink(dds_c, coef="condition_AS_vs_AO", res=res_AS_AO, 
                           type="apeglm")

# Aplicamos la función summary() al objeto DESeqResults generado tras aplicar la
# función lfcShrink(), donde hemos utilizado el argumento res para utilizar un 
# objeto DESeqResults donde hemos especificado los valores cutoff de alpha y log2
# fold changes. Aunque aparece LFC > 0 y LFC < 0 en vez de 0.26 y -0.26, los 
# resultados hacen referencia a estos últimos valores:
summary(res_AS_AO_LFC)

# Comprobamos que si cambiamos 0.26 por otro valor, los resultados cambian y sigue
# apareciendo LFC > 0 y LFC < 0:
res_AS_AO_Threshold <- results(dds_c, contrast=c("condition", "AS", "AO"), 
                               lfcThreshold=0.50, alpha=0.05)
res_AS_AO_LFC_Threshold <- lfcShrink(dds_c, coef="condition_AS_vs_AO", 
                                     res=res_AS_AO_Threshold, type="apeglm")
summary(res_AS_AO_LFC_Threshold)

# Aunque el parámetro coef en lfcShrink() es similar a contrast en results(), solo
# podemos llevar a cabo contrastes con los coeficientes que aparecen en resultsNames(). 
# Por defecto, el nivel de referencia de la variable subcondition es DS, por lo que,
# para poder comparar AS.REQneg y AS.REQpos, tenemos que seleccionar como nivel de
# referencia uno de estos niveles (para obtener los mismos resultados que anteriormente,
# seleccionamos AS.REQpos). Para ello, utilizamos la función relevel() y posteriormente
# utilizamos DESeq() para que se apliquen los cambios tras el cambio de 
# nivel de referencia.
dds_sc$subcondition <- relevel(dds_sc$subcondition, ref = "AS.REQpos")
dds_sc <- DESeq(dds_sc)
resultsNames(dds_sc)
res_ASREQneg_pos_LFC <- lfcShrink(dds_sc, coef="subcondition_AS.REQneg_vs_AS.REQpos", 
                                  res=res_ASREQneg_pos, type="apeglm")
summary(res_ASREQneg_pos_LFC)

# Ejemplo para ver que los p-valores y p-valores ajustados no cambian al trabajar 
# con shrunken log2 fold changes, pero los shrunken log2 fold changes son diferentes
# a log2 fold changes:
head(subset(res_ASREQneg_pos, padj<0.05 & abs(log2FoldChange)>0.26),3)
head(subset(res_ASREQneg_pos_LFC, padj<0.05 & abs(log2FoldChange)>0.26),3)

# Utilizamos la función par(mfrow=c(1,2)) para dividir la pantalla en 1 fila y 2
# columnas, y representar los dos gráficos uno al lado del otro:
par(mfrow=c(1,2))

# Utilizamos la función plotMA() para obtener los MA-plots con log2 fold changes
# o shrunken log2 fold changes en el eye y (comparando las condiciones AS y AO).
# Utilizamos abline() para añadir dos líneas rojas horizontales que marcan los 
# valores de log2 fold changes o shrunken log2 fold changes de 0.26 y -0.26:
plotMA(res_AS_AO, ylab="log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

plotMA(res_AS_AO_LFC, ylab="shrunken log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

# Utilizamos la función mtext() par añadir un título compartido por ambos gráficos:
mtext(text="condition_AS_vs_AO", side=3, line=-2, outer=TRUE)

# Repetimos el procedimiento anterior para obtener los MA-plots con log2 fold changes
# o shrunken log2 fold changes en el eye y (comparando las subcondiciones AS.REQneg y
# AS.REQpos):
par(mfrow=c(1,2))

plotMA(res_ASREQneg_pos, ylab="log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

plotMA(res_ASREQneg_pos_LFC, ylab="shrunken log2 fold changes")
abline(h=c(-0.26, 0.26), col="red", lwd=2)

mtext(text="subcondition_AS.REQneg_vs_AS.REQpos", side=3, line=-2, outer=TRUE)

library(ggplot2)

# Utilizamos la función vst() para llevar a cabo la transformación logarítmica 
# regularizada (utilizamos el argumento blind=FALSE ya que no requerimos que la 
# información de las muestras especificada mediante el argumento design no se 
# tenga en cuenta):
vsd_c <- vst(dds_c, blind=FALSE)
vsd_sc <- vst(dds_sc, blind=FALSE)
# Se trata de objetos DESeqTransform:
class(vsd_c)

# Utilizamos la función normTransform() para calcular log2 scaled counts y assay()
# para acceder a dichos valores. Utilizamos meanSdPlot() para representar la 
# desviación estándar de log2 scaled counts para cada isoTargets frente a la media
# (utilizamos el argumento ranks=FALSE para utilizar en el eje x la media y 
# plot=FALSE para no representar el gráfico y añadir posteriormente el título):
meanSd_log2_scaled_count <- vsn::meanSdPlot(assay(normTransform(dds_c)), ranks=FALSE, 
                                            plot=FALSE)

# Representamos el gráfico con el título:
meanSd_log2_scaled_count$gg + ggtitle("log2 scaled counts")

# Utilizamos meanSdPlot() para representar la desviación estándar de regularized 
# log transformation counts para cada isoTargets frente a la media (de nuevo, 
# accedemos a los valores utilizando assay()):
meanSd_reg_log_transf_counts <- vsn::meanSdPlot(assay(vsd_c), ranks=FALSE, plot=FALSE)

# Representamos el gráfico con el título:
meanSd_reg_log_transf_counts$gg + ggtitle("Regularized log transformation counts")

# Cargamos el paquete cowplot para utilizar la función plot_grid():
library(cowplot)

# Llevamos a cabo la representación de baja dimensión de los datos utilizando como
# conjunto de variables iniciales todos los isoTargets/aquellos con un |LFC|>0.26 
# estadísticamente significativos (condiciones o grupos AO y AS) (utilizamos los
# data.frame isotargets_sig_c e isotargets_sig_sc definidos anteriormente):
PCA_c <- plotPCA(vsd_c[,c(1:13)], intgroup="condition")
PCA_sig_c <- plotPCA(vsd_c[rownames(isotargets_sig_c), c(1:13)], intgroup="condition")

# Utilizamos geom_text() para añadir y modificar las etiquetas de las muestras y 
# coord_fixed() para cambiar los límites de los ejes:
PCA_c <- PCA_c + geom_text(aes(label = rownames(coldata)[1:13]), size=2.5, angle=15,
                           color="black") + coord_fixed(ylim=c(-22, 22), xlim=c(-35, 35))
PCA_sig_c <- PCA_sig_c + geom_text(aes(label = rownames(coldata)[1:13]), size=2.5, 
                                   angle=15, color="black") + coord_fixed(ylim=c(-15, 17), xlim=c(-31, 27))

# Utilizamos la función plot_grid() para dividir la pantalla en 2 filas y 1
# columna, y representar los dos gráficos uno encima del otro (no utilizamos la
# función par() porque no son base graphics):
plot_grid(PCA_c, PCA_sig_c, nrow=2)

# Repetimos el procedimiento anterior con las subcondiciones o grupos AS.REQneg y
# AS.REQpos:
PCA_sc <- plotPCA(vsd_sc[,c(5:13)], intgroup="subcondition")
PCA_sig_sc <- plotPCA(vsd_sc[rownames(isotargets_sig_sc), c(5:13)], 
                      intgroup="subcondition")

PCA_sc <- PCA_sc + geom_text(aes(label = rownames(coldata)[5:13]), size=2.5, angle=70,
                             color="black") + coord_fixed(ylim=c(-20, 20), xlim=c(-24, 25))
PCA_sig_sc <- PCA_sig_sc + geom_text(aes(label = rownames(coldata)[5:13]), size=2.5, 
                                     angle=70, color="black") + coord_fixed(ylim=c(-8, 7), xlim=c(-11, 11))

plot_grid(PCA_sc, PCA_sig_sc, nrow=2)

library(pheatmap)

# Creamos dos data.frame con los nombres de las muestras a comparar (AS y AO, y 
# AS.REQneg y AS.REQpos):
df_c <- as.data.frame(colData(dds_c)[1:13,][1])
df_sc <- as.data.frame(colData(dds_sc)[5:13,][2])

# Obtenemos el heatmap para las muestras AS y AO utilizando los argumentos 
# descritos anteriormente:
pheatmap(assay(vsd_c)[rownames(isotargets_sig_c), c(1:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_c, scale="row")

# Obtenemos el heatmap para las muestras AS.REQneg y AS.REQpos:
pheatmap(assay(vsd_sc)[rownames(isotargets_sig_sc), c(5:13)], cluster_rows=TRUE, 
         cluster_cols=TRUE, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", clustering_method="complete", 
         show_rownames=FALSE, annotation_col=df_sc, scale="row")

# Análisis de expresión diferencial a nivel de miRNA e isomiR
# En este caso countData corresponde con isomir_count (el nombre de las filas 
# corresponde con la columna UID.x):
isomir_count <- isomir_subset[,6:22]

# Mostramos los primeros elementos:
head(isomir_count)

# Construimos los objetos DESeqDataSet:
dds_isomir_c <- DESeqDataSetFromMatrix(countData=isomir_count, colData=coldata,
                                       design=~condition)

dds_isomir_sc <- DESeqDataSetFromMatrix(countData=isomir_count, colData=coldata,
                                        design=~subcondition)

# Llevamos a cabo el análisis de expresión diferencial:
dds_isomir_c <- DESeq(dds_isomir_c)
dds_isomir_sc <- DESeq(dds_isomir_sc)

# Mostramos los resultados al comparar los grupos AS y AO:
res_isomir_AS_AO <- results(dds_isomir_c, contrast=c("condition", "AS", "AO"), 
                            lfcThreshold=0.26, alpha=0.05, altHypothesis="greaterAbs")
summary(res_isomir_AS_AO)

# Mostramos los resultados al comparar los grupos AS.REQneg y AS.REQpos:
res_isomir_ASREQneg_pos <- results(dds_isomir_sc, contrast=c("subcondition", "AS.REQneg", 
                                                             "AS.REQpos"), lfcThreshold=0.26, alpha=0.05, 
                                   altHypothesis="greaterAbs")
summary(res_isomir_ASREQneg_pos)

# Utilizamos la función EnhancedVolcano() con distintos argumentos para obtener
# el volcano plot de las condiciones o grupos AS y AO:
EnhancedVolcano(res_isomir_AS_AO, lab="", x="log2FoldChange", y="padj", pCutoff=0.05, 
                FCcutoff=0.26, title="isomiR_condition_AS_vs_AO", ylim=c(0,7), 
                xlim=c(-2,7), legendLabSize=12, legendIconSize=5, 
                legendPosition="right")

# Repetimos el procedimiento anterior para obtener el volcano plot de las 
# subcondiciones AS.REQneg y AS.REQpos:
EnhancedVolcano(res_isomir_ASREQneg_pos, lab="", x="log2FoldChange", y="padj", 
                pCutoff=0.05, FCcutoff=0.26, title="isomiR_subcondition_ASREQneg_pos", 
                ylim=c(0,5), xlim=c(-5, 6), legendLabSize=12, legendIconSize=5, 
                legendPosition = "right")

# Seleccionamos aquellos miRNA/isomiR con un |LFC|>0.26 estadísticamente significativos:
isomir_sig_c <- subset(res_isomir_AS_AO, padj<0.05)
isomir_sig_sc <- subset(res_isomir_ASREQneg_pos, padj<0.05)

# De 211 y 70 miRNA/isomiR, vemos cuántos son isomiR y cuántos miRNA. Comenzamos por
# las condiciones AS y AO:
as.data.frame(isomir_sig_c) %>%
  mutate(mirna_canonico = row.names(isomir_sig_c) %in% subset(isomir_subset$UID.x,
                                                              isomir_subset$variant == "NA")) %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Repetimos el mismo procedimiento en el caso de las subcondiciones AS.REQneg y 
# AS.REQpos:
as.data.frame(isomir_sig_sc) %>%
  mutate(mirna_canonico = row.names(isomir_sig_sc) %in% subset(isomir_subset$UID.x,
                                                               isomir_subset$variant == "NA")) %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Análisis de expresión diferencial a nivel de miRNA
# Cargamos el paquete tibble para utilizar la función column_to_rownames():
library(tibble)

# En este caso countData corresponde con mirna_count (el nombre de las filas 
# corresponde con la columna parent). Agrupamos los contajes de cada miRNA:
mirna_count <- as.data.frame(isomir_subset %>%
                               group_by(parent) %>%
                               summarize_at(vars(5:21), sum) %>%
                               ungroup() %>%
                               column_to_rownames("parent"))

# Mostramos los primeros elementos:
head(mirna_count)

# Construimos los objetos DESeqDataSet:
dds_mirna_c <- DESeqDataSetFromMatrix(countData=mirna_count, colData=coldata,
                                      design=~condition)

dds_mirna_sc <- DESeqDataSetFromMatrix(countData=mirna_count, colData=coldata,
                                       design=~subcondition)

# Llevamos a cabo el análisis de expresión diferencial:
dds_mirna_c <- DESeq(dds_mirna_c)
dds_mirna_sc <- DESeq(dds_mirna_sc)

# Mostramos los resultados al comparar los grupos AS y AO:
res_mirna_AS_AO <- results(dds_mirna_c, contrast=c("condition", "AS", "AO"), 
                           lfcThreshold=0.26, alpha=0.05, altHypothesis="greaterAbs")
summary(res_mirna_AS_AO)

# Mostramos los resultados al comparar los grupos AS.REQneg y AS.REQpos:
res_mirna_ASREQneg_pos <- results(dds_mirna_sc, contrast=c("subcondition", "AS.REQneg", 
                                                           "AS.REQpos"), lfcThreshold=0.26, alpha=0.05, 
                                  altHypothesis="greaterAbs")
summary(res_mirna_ASREQneg_pos)

# Utilizamos la función EnhancedVolcano() con distintos argumentos para obtener
# el volcano plot de las condiciones o grupos AS y AO:
EnhancedVolcano(res_mirna_AS_AO, lab="", x="log2FoldChange", y="padj", pCutoff=0.05, 
                FCcutoff=0.26, title="miRNA_condition_AS_vs_AO", ylim=c(0,7), 
                xlim=c(-2,7), legendLabSize=12, legendIconSize=5, 
                legendPosition="right")

# Repetimos el procedimiento anterior para obtener el volcano plot de las 
# subcondiciones AS.REQneg y AS.REQpos:
EnhancedVolcano(res_mirna_ASREQneg_pos, lab="", x="log2FoldChange", y="padj", 
                pCutoff=0.05, FCcutoff=0.26, title="miRNA_subcondition_ASREQneg_pos", 
                ylim=c(0,5), xlim=c(-5, 6), legendLabSize=12, legendIconSize=5, 
                legendPosition = "right")

# Seleccionamos aquellos isoTargets con un |LFC|>0.26 estadísticamente significativos:
mirna_sig_c <- subset(res_mirna_AS_AO, padj<0.05)
mirna_sig_sc <- subset(res_mirna_ASREQneg_pos, padj<0.05)

# Los 39 miRNA diferencialmente expresados en las condiciones AS y AO corresponden
# con 30 miRNA y 202 isomiR:
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(parent %in% row.names(mirna_sig_c)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

# Repetimos el mismo procedimiento en el caso de las subcondiciones AS.REQneg y 
# AS.REQpos. Los 26 miRNA diferencialmente expresados en estas subcondiciones
# corresponden con 19 miRNA y 166 isomiR:
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(parent %in% row.names(mirna_sig_sc)) %>%
  mutate(mirna_canonico = variant == "NA") %>%
  summarize(count_mirna_canonico = sum(mirna_canonico == TRUE),
            count_isomir = sum(mirna_canonico == FALSE))

library(stringr)

# Seleccionamos los miRNA/isomiR que forman parte de isoTargets que engloban más
# de 1 miRNA/isomiR (hay isoTargets con un único miRNA/isomiR):
isomir_isotargets_multiple <- as.data.frame(isomir_arrange %>%
                                              group_by(transcript_code) %>%
                                              mutate(isotargets_multiple=n()) %>%
                                              filter(isotargets_multiple>1))

# Tenemos 344 isoTargets que engloban más de 1 miRNA/isomiR:
dim(as.data.frame(isomir_isotargets_multiple %>%
                    group_by(transcript_code) %>%
                    summarize(count_distinct_parent=n_distinct(parent))))[1]

# Estos 344 isoTargets incluyen 1172 miRNA/isomiR de 2705 totales:
dim(isomir_isotargets_multiple)[1]

# Solo 5 isomiR que corresponden con variantes en la región 5'UTR (de 179 totales)
# se han agrupado en un isoTargets que englobe más de 1 miRNA/isomiR (aproximadamente
# un 3 %):
dim(as.data.frame(isomir_isotargets_multiple %>%
                    filter(str_detect(variant, "5p"))))[1]

dim(as.data.frame(isomir_arrange %>%
                    filter(str_detect(variant, "5p"))))[1]

# Sin embargo, 961 isomiR que corresponden con variantes en la región 3'UTR (de 
# 1788 totales) se han agrupado en un isoTargets que englobe más de 1 miRNA/isomiR
# (53.75 %):
dim(as.data.frame(isomir_isotargets_multiple %>%
                    filter(str_detect(variant, "3p"))))[1]

dim(as.data.frame(isomir_arrange %>%
                    filter(str_detect(variant, "3p"))))[1]

# Comprobamos si el isoTargets al que pertenece cada uno de estos 5 isomiR también
# contiene al miRNA canónico. Primero, seleccionamos todos los miRNA/isomiR que 
# forman parte de los isoTargets a los que pertenecen los 5 isomiR que corresponden
# con variantes en la región 5'UTR que se han agrupado en un isoTargets que englobe 
# más de 1 miRNA/isomiR:
semi_join(isomir_arrange, (as.data.frame(isomir_isotargets_multiple %>%
                                           filter(str_detect(variant, "5p")))), by="transcript_code") %>%
  # Agrupamos por transcript_code y parent y creamos una nueva columna 
  # (canonical_isotargets) que corresponde con la suma de todas las variantes que
  # hay en el grupo transcript_code-parent si en este está el miRNA canónico:
  group_by(transcript_code, parent) %>%
  summarize(canonical_isotargets = sum(str_detect(variant, "5p") & 
                                         any(variant == "NA"))) %>%
  # Finalmente desagrupamos y sumamos todos los resultados obtenidos (0):
  ungroup() %>%
  summarize(sum(canonical_isotargets))

# De las 179 variantes en la región 5'UTR, el data.frame isomir_arrange contiene 
# el miRNA canónico de 158. Según los resultados obtenidos, los 158 isomiR que 
# corresponden con variantes en la región 5'UTR y sus respectivos miRNA canónicos
# están en isoTargets distintos (100%):
isomir_arrange %>%
  filter(str_detect(variant, "5p")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Repetimos el mismo procedimiento con los isomiR que corresponden con variantes 
# en la región 3'UTR. 
semi_join(isomir_arrange, (as.data.frame(isomir_isotargets_multiple %>%
                                           filter(str_detect(variant, "3p")))), by="transcript_code") %>%
  group_by(transcript_code, parent) %>%
  summarize(canonical_isotargets = sum(str_detect(variant, "3p") & 
                                         any(variant == "NA"))) %>%  
  # En este caso obtenemos 518 isomiR que son variantes en la región 3'UTR que han
  # sido agrupados con su miRNA canónico en un mismo isoTargets:
  ungroup() %>%
  summarize(sum(canonical_isotargets))

# De esas 1788 variantes en la región 3'UTR, el data.frame isomir_arrange contiene 
# el miRNA canónico de 1651. Según los resultados obtenidos, 1651-518 son los isomiR
# que corresponden con variantes en la región 3'UTR para los cuales sus respectivos
# miRNA canónicos están en isoTargets distintos (68.63 %):
isomir_arrange %>%
  filter(str_detect(variant, "3p")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Ningún isomiR que corresponde con una variante iso_snv_seed (de 108 totales) se
# ha agrupado en un isoTargets que englobe más de 1 miRNA/isomiR (0 %):
dim(as.data.frame(isomir_isotargets_multiple %>%
                    filter(str_detect(variant, "snv_seed"))))[1]

dim(as.data.frame(isomir_arrange %>%
                    filter(str_detect(variant, "snv_seed"))))[1]

# De esas 108 variantes iso_snv_seed, el data.frame isomir_arrange contiene el 
# miRNA canónico de 100. Según los resultados obtenidos, todos los isomiR que 
# corresponden con variantes iso_snv_seed están en isoTargets distintos con respecto
# al isoTargets de su miRNA canónico correspondiente (100 %):
isomir_arrange %>%
  filter(str_detect(variant, "iso_snv_seed")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Sin embargo, 31 isomiR que corresponden con variantes iso_snv, sin contar las
# variantes iso_snv_seed (de 308 totales) se han agrupado en un isoTargets que 
# englobe más de 1 miRNA/isomiR (10.06 %):
dim(as.data.frame(isomir_isotargets_multiple %>%
                    filter(str_detect(variant, "iso_snv")) %>%
                    filter(!str_detect(variant, "iso_snv_seed"))))[1]

dim(as.data.frame(isomir_arrange %>%
                    filter(str_detect(variant, "iso_snv")) %>%
                    filter(!str_detect(variant, "iso_snv_seed"))))[1]

# Comprobamos si el isoTargets al que pertenece cada uno de estos 31 isomiR iso_snv
# también contiene al miRNA canónico:
semi_join(isomir_arrange, (as.data.frame(isomir_arrange %>%
                                           filter(str_detect(variant, "iso_snv")) %>%
                                           filter(!str_detect(variant, "iso_snv_seed")))), by="transcript_code") %>%
  group_by(transcript_code, parent) %>%
  summarize(canonical_isotargets = sum(str_detect(variant, "iso_snv") & 
                                         any(variant == "NA"))) %>%
  # Finalmente desagrupamos y sumamos todos los resultados obtenidos (12):
  ungroup() %>%
  summarize(sum(canonical_isotargets))

# De esas 308 variantes iso_snv (sin incluir iso_snv_seed), el data.frame 
# isomir_arrange contiene el miRNA canónico de 303. Según los resultados obtenidos,
# 303-12 son los isomiR que corresponden con variantes iso_snv (sin incluir 
# iso_snv_seed) para los cuales sus respectivos miRNA canónicos están en isoTargets
# distintos (96.04 %):
isomir_arrange %>%
  filter(str_detect(variant, "iso_snv")) %>%
  filter(!str_detect(variant, "iso_snv_seed")) %>%
  summarize(canonical = sum(parent %in% subset(isomir_arrange$parent, 
                                               isomir_arrange$variant == "NA")))

# Utilizamos el argumento mar en la función par() para modificar los márgenes:
par(mfrow=c(2,2), mar=c(2.5,2.5,2.5,2.5))

# Creamos los gráficos de tarta con la función pie():
colors <- c("#99FF99", "#FFD699")

five <- 100
cat_five <- "diff_isoTargets"

pie(five, labels=cat_five, col=colors, cex=0.75, main="5p")

three <- c(68.63, 100-68.63)
cat_three <- c("diff_isoTargets", "same_isoTargets")

pie(three, labels=cat_three, col=colors, cex=0.75, main="3p")

snv_seed <- 100
cat_snv_seed <- "diff_isoTargets"

pie(snv_seed, labels=cat_snv_seed, col=colors, cex=0.75, main="snv_seed")

snv_noseed <- c(96.04, 100-96.04)
cat_snv_noseed <- c("diff_isoTargets", "same_isoTargets")

pie(snv_noseed, labels=cat_snv_noseed, col=colors, cex=0.75, main="snv_noseed")

# De los 200  miRNA/isomiR que forman parte de los 112 isoTargets diferencialmente
# expresados en las condiciones AS y AO, 182 coinciden con los obtenidos utilizando
# el enfoque a nivel de miRNA/isomiR (91 %):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
  summarize(count(UID.x %in% row.names(isomir_sig_c)))

# Es decir, estos 182 miRNA/isomiR representan el 86.26 % de los miRNA/isomiR
# diferencialmente expresados obtenidos con el enfoque a nivel de miRNA/isomiR:
182/211*100

# De los 58 miRNA/isomiR que forman parte de los 34 isoTargets diferencialmente
# expresados en las subcondiciones AS.REQneg y AS.REQpos, 39 coinciden con los 
# obtenidos utilizando el enfoque a nivel de miRNA/isomiR (67.24 %):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
  summarize(count(UID.x %in% row.names(isomir_sig_sc)))

# Es decir, estos 39 miRNA/isomiR representan el 55.71 % de los miRNA/isomiR
# diferencialmente expresados obtenidos con el enfoque a nivel de miRNA/isomiR:
39/70*100

# De los 200  miRNA/isomiR que forman parte de los 112 isoTargets diferencialmente
# expresados en las condiciones AS y AO, 190 coinciden con los obtenidos utilizando
# el enfoque a nivel de miRNA (95 %):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_c)) %>%
  summarize(count(parent %in% row.names(mirna_sig_c)))

# Es decir, estos 190 miRNA/isomiR representan el 81.90 % de los miRNA/isomiR
# diferencialmente expresados obtenidos con el enfoque a nivel de miRNA:
190/232*100

# De los 58 miRNA/isomiR que forman parte de los 34 isoTargets diferencialmente
# expresados en las subcondiciones AS.REQneg y AS.REQpos, 44 coinciden con los 
# obtenidos utilizando el enfoque a nivel de miRNA/isomiR (75.86 %):
isomir_subset_rep1_5[,c(2:4,24)] %>%
  filter(transcript_code %in% row.names(isotargets_sig_sc)) %>%
  summarize(count(parent %in% row.names(mirna_sig_sc)))

# Es decir, estos 44 miRNA/isomiR representan el 23.78 % de los miRNA/isomiR
# diferencialmente expresados obtenidos con el enfoque a nivel de miRNA:
44/185*100

# Comparamos el p-valor y el p-valor ajustado obtenidos para cada isoTargets 
# diferencialmente expresado con la media de p-valores y p-valores ajustados 
# obtenidos para cada miRNA/isomiR individual.

# Seleccionamos aquellos isoTargets diferencialmente expresados compuestos por más
# de 1 isomiR/miRNA y añadimos la columna UID.x (se duplicarán las filas para cada
# isoTargets para tener una fila por cada miRNA/isomiR individual que forma parte
# de dicho grupo):
isotargets_sig_c_UID <- as.data.frame(isotargets_sig_c) %>%
  filter(row.names(isotargets_sig_c) %in% isomir_isotargets_multiple$transcript_code) %>%
  merge(isomir_isotargets_multiple %>% 
          select(transcript_code, UID.x),
        by.x = 0, by.y = "transcript_code", all.x = TRUE) %>%
  mutate(transcript_code=Row.names) %>%
  select(-Row.names)

# Mostramos los primeros elementos:
head(isotargets_sig_c_UID, 5)

# Creamos un data.frame con las columnas transcript_code, pvalue y padj de aquellos
# isoTargets diferencialmente expresados compuestos por más de 1 isomiR/miRNA:
isotargets_multiple_sig_c <- as.data.frame(isotargets_sig_c) %>%
  filter(row.names(isotargets_sig_c) %in% isomir_isotargets_multiple$transcript_code) %>%
  rownames_to_column("transcript_code") %>%
  select(c(transcript_code, pvalue, padj))

# Creamos un data.frame con las columnas transcript_code, pvalue_mean y padj_mean 
# con los valores medios de los p-valores y p-valores ajustados de aquellos 
# miRNA/isomiR (enfoque a nivel de miRNA/isomiR) que forman parte de los isoTargets
# diferencialmente expresados compuestos por más de 1 isomiR/miRNA:
isomir_multiple_c <- as.data.frame(as.data.frame(res_isomir_AS_AO) %>%
                                     rownames_to_column("UID.x") %>%
                                     inner_join(isotargets_sig_c_UID, by = c("UID.x"="UID.x")) %>%
                                     group_by(transcript_code) %>%
                                     summarize(pvalue_mean = mean(pvalue.x, na.rm=TRUE), padj_mean=mean(padj.x, 
                                                                                                        na.rm=TRUE)))

# Comparamos los p-valores con los p-valores medios, y los p-valores ajustados con
# los p-valores ajustados medios para cada isoTargets:
isotargets_multiple_sig_c %>%
  inner_join(isomir_multiple_c, by="transcript_code") %>%
  mutate(pvalue_isotargets_low=pvalue<pvalue_mean, padj_isotargets_low=padj<padj_mean)

