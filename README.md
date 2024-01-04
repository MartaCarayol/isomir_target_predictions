# isomir_target_predictions

This repository contains two scripts for predicting mRNA targets from a list of customized miRNA/isomiR sequences using [scanMiR](https://github.com/ETHZ-INS/scanMiR/tree/master) and [miRDB](https://mirdb.org/custom.html). Additionally, it also presents a novel approach to perform differential expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (differential expression analysis at the isoTargets level).

## Requirements

- Trained model from [github.com/kslin/miRNA_models/tree/master/cnn](https://github.com/kslin/miRNA_models/tree/master/cnn): generate_12_mer_kds.py, model-100.data-00000-of-00001, model-100.index, model-100.meta and utils.py. Ensure that these files are located in the directory specified as the first argument in script_scanmir.

- Python packages for automated_mirdb_mirna.py: 

	- biopython
  	- selenium
 	- pandas
 	- beautifulsoup4
	- argparse
	- logging
	- sqlite3

- R packages for script_scanmir.R: 

	- biomaRt (Bioconductor package)
	- BSgenome.Hsapiens.NCBI.GRCh38 (Bioconductor package)
	- EnsDb.Hsapiens.v86 (Bioconductor package)
	- reticulate
	- scanMiR (Bioconductor package)
	- foreach
	- doParallel
	- doSNOW

- Python packages for script_scanmir.R:

	- optparse
	- itertools
	- numpy
	- pandas
	- tensorflow.compat.v1 (module for TensorFlow 1.x compatibility functions)
 	- requests
	- utils (custom module from [github.com/kslin/miRNA_models/blob/master/cnn/utils.py](https://github.com/kslin/miRNA_models/blob/master/cnn/utils.py))

- R packages for clash_scanmir_mirdb_reproducibility.R:

	- biomaRt (Bioconductor package)
	- BSgenome.Hsapiens.NCBI.GRCh38 (Bioconductor package)
	- EnsDb.Hsapiens.v86 (Bioconductor package)
	- tidyr
	- reticulate
	- scanMiR (Bioconductor package)
	- dplyr

- R packages for differential_expression_analysis_isotargets.R:

	- dplyr
	- digest
	- knitr
	- DESeq2 (Bioconductor package)
	- EnhancedVolcano
	- ggplot2
	- cowplot
	- pheatmap
	- tibble
	- stringr

For example, to install the R packages for differential_expression_analysis_isotargets.R, use ```install.packages(c("dplyr", "digest", "knitr", "EnhancedVolcano", "ggplot2", "cowplot", "pheatmap", "tibble", "stringr"))``` and for Bioconductor package DESeq2:

```
Install DESeq2 from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

To install a Python package using conda in the base Anaconda environment, use the ```conda install``` command followed by the package name. For example, ```conda install numpy```.

This project has been developed using Python 3.9.18 and R 4.3.0. 

## Modules

### script_scanmir.R

This module contains an R script developed to simplify the usage of [scanMiR](https://github.com/ETHZ-INS/scanMiR/tree/master) for users with custom sequences. It employs the trained model from [github.com/kslin/miRNA_models/tree/master/cnn](https://github.com/kslin/miRNA_models/tree/master/cnn) described in [this paper](https://www.biorxiv.org/content/10.1101/414763v1), predicting relative KD values and defining a KdModelList object. The result is a .csv file with the number of target sites of each type found in the 3'UTR region for each protein-coding transcripts (grouped into 8mer, 7mer, 6mer and non-canonical), including the number of canonical and non-canonical target sites in other regions (ORF), as well as the repression value (result of applying findSeedMatches() function with a KdModelList object).

Here is an example of how to use the script employing the system command:

	system('Rscript script_scanmir.R C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM azoos_isomir_subset.csv 200 result_scanmir_azoos C:/Users/marta/anaconda3/python.exe')

The first argument is the directory where the previously mentioned files are located. The file azoos_isomir_subset.csv, provided as the second argument, is a CSV file with two columns. The first column contains the names of miRNA/isomiR, while the second column contains their respective sequences. The third argument represents the desired number of target transcripts for each miRNA/isomiR. The script will select targets with the lowest repression values. The fourth argument is the name of the resulting .csv file with the analysis outcomes. As an optional parameter, you can include the absolute path to the Python version you wish to use, provided as the fifth argument. 

 ### automated_mirdb_mirna.py

This script implements a Selenium WebDriver to automate searches on [miRDB - MicroRNA Target Prediction and Functional Study Database](http://mirdb.org/). Its purpose is to automate the prediction of mRNA targets using a user-supplied list of miRNA/isomiR.

Here is an example of how to use the script employing the system command:

	system('python automated_mirdb_mirna.py data_clash_objetivo1.txt result_mirdb Human', wait=TRUE)

The arguments include a FASTA file with miRNA/isomiR sequences, a name for output .csv file and species (the options are Human, Rat, Mouse, Chicken and Dog).

 ### clash_scanmir_mirdb_reproducibility.R

Compare the results obtained for 5 miRNA from experimental, scanMiR, and miRDB analyses. Consider factors such as the scanned regions (3'UTR, CDS, and/or 5'UTR), the inclusion of protein-coding transcripts, the involvement of seed region and the type of target site identified (8mer, 7mer, 6mer, etc.).

 ### differential_expression_analysis_isotargets.R

Apply script_scanmir to the data from [Larriba et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37245055/). Assig a unique code to each set of target mRNA, forming isoTargets groups (canonical miRNA and isomiR regulating the same mRNA share the same code). Then, pooling the reads of all canonical miRNAs and isomiRs within each isoTargets group across all samples for the final differential expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (differential expression analysis at the isoTargets level).

