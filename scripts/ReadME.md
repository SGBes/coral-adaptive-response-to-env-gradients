1. [Instructions to install ANGSD](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/install_angsd.txt)

2. [Set up linux environment](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/bash_1_set%20up%20environment%20and%20directories.txt)

3. [Generate reference genome](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/bash_2_reference%20genome.txt)

4. Download sequence reads from SRA run selector, generate fixation index (Fst) and gene count data
   * linux script: "bash_3_(env gradient).txt"
   * [generate metadata for linux script and DESeq2](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/processing_fastq.R)
   * [check for clones](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/remove_clones.R)

5. [Convert Fst data to per gene](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/Fst_toPerGene.R)

6. [Use DESeq2 to calculate log-fold change from gene count data](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/DESeq2_logFoldChange.R)

7. [Multivariate analysis and plots](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/multivariate_analysis.R)

8. [Include general stress response - Red Module from Dixon et al., 2020 - in gene set overlap analysis for gene expression](https://github.com/SGBes/coral-adaptive-response-to-env-gradients/blob/main/scripts/multivariate_analysis_RedModule.R)
