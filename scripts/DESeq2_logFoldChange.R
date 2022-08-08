# set working directory
setwd("~/Desktop/analysis")

# load library
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(ape)
library(DESeq2)
library(cowplot)
rm(list=ls())

## heat stress ####

# upload read counts from htseq
counts = read.delim("feature_counts_out_heat_stress.txt", comment.char="#")
row.names(counts)=counts$Geneid
counts = counts[,7:ncol(counts)]
head(counts)
dim(counts)

# look at read count sums 
dim(counts)
tots = apply(counts, 2, sum)
hist(tots)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))

# remove genes with mean read count below cutoff
cut = 3 # this if the gene count cutoff
cc = counts
means = apply(cc,1,mean)
# number of genes that have a mean below the cutoff
table(means>cut)
counts = cc[means>cut,]

# BUILD TRAIT DATA TABLE
# table with rownames equivalent to the column names in the counts data frame
# can have as many columns as sample attributes

# upload and reformat sample info
# rows as sample names and traits as columns

# this is the metadata csv
tdat = read.csv("Metadata_heat_stress.csv", header = T)[c('num','treat')]
tdat = tdat %>% 
  mutate(num = sub("$", ".bam", num))
rownames(tdat) = tdat$num
head(tdat)

# subset for the samples with RNAseq data
coldata = tdat %>% 
  filter(num %in% colnames(counts))
dim(coldata)

# check the rownames and colnames match up
sum(rownames(coldata)==colnames(counts))==ncol(counts)

# make treatment factors
# note: make sure "control" is the 1st level in the treatment factor
str(coldata)
coldata$treat = as.factor(coldata$treat)
coldata$treat = relevel(coldata$treat, "MV")
str(coldata)
levels(coldata$treat)

#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
# set up input matrix for DESeq
ddsHTSeq=DESeqDataSetFromMatrix(counts,
                                colData = coldata,
                                # this formula is just based on 1 treatment which
                                # is the environmental characteristic
                                design = formula(~treat))

# run DESeq
dds = DESeq(ddsHTSeq)

# log 2 fold change between these the control and treatment environmental parameter
resultsNames(dds)
# extract desired comparison
res = results(dds,name="treat_HV_vs_MV")

# retreiving summary of results
# default p-val = 0.1
res01 = results(dds, alpha = 0.1, name="treat_HV_vs_MV")
summary(res01)

# How many adjusted p-values were less than 0.01?
sum(res01$padj < 0.1, na.rm=TRUE)

# volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-7,7),
                ylim = c(0,15),
                legendPosition = 'bottom')

# extracting log-fold change
# LFC >0 -- up expression
# LFC <0 -- down expression
res
res.df = as.data.frame(row.names(res))
res.df$LFC = res$log2FoldChange
res.df = res.df %>% 
  dplyr::rename("gene"="row.names(res)")

write.csv(res.df,file="DESeq2_LFC_heat_stress.csv")

## co2 dobu ####

# upload read counts from htseq
counts = read.delim("feature_counts_out_co2_dobu.txt", comment.char="#")
row.names(counts)=counts$Geneid
counts = counts[,7:ncol(counts)]
head(counts)
dim(counts)

# look at read count sums 
dim(counts)
tots = apply(counts, 2, sum)
hist(tots)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))

# remove genes with mean read count below cutoff
cut = 3 # this if the gene count cutoff
cc = counts
means = apply(cc,1,mean)
# number of genes that have a mean below the cutoff
table(means>cut)
counts = cc[means>cut,]

# BUILD TRAIT DATA TABLE
# table with rownames equivalent to the column names in the counts data frame
# can have as many columns as sample attributes

# upload and reformat sample info
# rows as sample names and traits as columns

# this is the metadata csv
tdat = read.csv("Metadata_co2_dobu.csv", header = T)[c('num','treat')]
tdat = tdat %>% 
  mutate(num = sub("$", ".bam", num)) %>% 
  arrange(num)
rownames(tdat) = tdat$num
head(tdat)

# subset for the samples with RNAseq data
coldata = tdat %>% 
  filter(num %in% colnames(counts))
dim(coldata)

# check the rownames and colnames match up
sum(rownames(coldata)==colnames(counts))==ncol(counts)

# make treatment factors
# note: make sure "control" is the 1st level in the treatment factor
str(coldata)
coldata$treat = as.factor(coldata$treat)
coldata$treat = relevel(coldata$treat, "Control")
str(coldata)
levels(coldata$treat)

#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
# set up input matrix for DESeq
ddsHTSeq=DESeqDataSetFromMatrix(counts,
                                colData = coldata,
                                # this formula is just based on 1 treatment which
                                # is the environmental characteristic
                                design = formula(~treat))

# run DESeq
dds = DESeq(ddsHTSeq)

# log 2 fold change between these the control and treatment environmental parameter
resultsNames(dds)
# extract desired comparison
res = results(dds,name="treat_Seep_vs_Control")

# retreiving summary of results
# default p-val = 0.1
res01 = results(dds, alpha = 0.1, name="treat_Seep_vs_Control")
summary(res01)

# How many adjusted p-values were less than 0.01?
sum(res01$padj < 0.1, na.rm=TRUE)

# volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-7,7),
                ylim = c(0,15),
                legendPosition = 'bottom')

# extracting log-fold change
# LFC >0 -- up expression
# LFC <0 -- down expression
res
res.df = as.data.frame(row.names(res))
res.df$LFC = res$log2FoldChange
res.df = res.df %>% 
  dplyr::rename("gene"="row.names(res)")

write.csv(res.df,file="DESeq2_LFC_co2_dobu.csv")

## co2 upa-upasina ####

# upload read counts from htseq
counts = read.delim("feature_counts_out_co2_upaUpasina.txt", comment.char="#")
row.names(counts)=counts$Geneid
counts = counts[,7:ncol(counts)]
head(counts)
dim(counts)

# look at read count sums 
dim(counts)
tots = apply(counts, 2, sum)
hist(tots)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))

# remove genes with mean read count below cutoff
cut = 3 # this if the gene count cutoff
cc = counts
means = apply(cc,1,mean)
# number of genes that have a mean below the cutoff
table(means>cut)
counts = cc[means>cut,]

# BUILD TRAIT DATA TABLE
# table with rownames equivalent to the column names in the counts data frame
# can have as many columns as sample attributes

# upload and reformat sample info
# rows as sample names and traits as columns

# this is the metadata csv
tdat = read.csv("Metadata_co2_upaUpasina.csv", header = T)[c('num','treat')]
tdat = tdat %>% 
  mutate(num = sub("$", ".bam", num)) %>% 
  arrange(num)
rownames(tdat) = tdat$num
head(tdat)

# subset for the samples with RNAseq data
coldata = tdat %>% 
  filter(num %in% colnames(counts))
dim(coldata)

# check the rownames and colnames match up
sum(rownames(coldata)==colnames(counts))==ncol(counts)

# make treatment factors
# note: make sure "control" is the 1st level in the treatment factor
str(coldata)
coldata$treat = as.factor(coldata$treat)
coldata$treat = relevel(coldata$treat, "Control")
str(coldata)
levels(coldata$treat)

#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
# set up input matrix for DESeq
ddsHTSeq=DESeqDataSetFromMatrix(counts,
                                colData = coldata,
                                # this formula is just based on 1 treatment which
                                # is the environmental characteristic
                                design = formula(~treat))

# run DESeq
dds = DESeq(ddsHTSeq)

# log 2 fold change between these the control and treatment environmental parameter
resultsNames(dds)
# extract desired comparison
res = results(dds,name="treat_Seep_vs_Control")

# retreiving summary of results
# default p-val = 0.1
res01 = results(dds, alpha = 0.1, name="treat_Seep_vs_Control")
summary(res01)

# How many adjusted p-values were less than 0.01?
sum(res01$padj < 0.1, na.rm=TRUE)

# volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-7,7),
                ylim = c(0,15),
                legendPosition = 'bottom')

# extracting log-fold change
# LFC >0 -- up expression
# LFC <0 -- down expression
res
res.df = as.data.frame(row.names(res))
res.df$LFC = res$log2FoldChange
res.df = res.df %>% 
  dplyr::rename("gene"="row.names(res)")

write.csv(res.df,file="DESeq2_LFC_co2_upaUpasina.csv")

## symbiont ####

# upload read counts from htseq
counts = read.delim("feature_counts_out_symbiont.txt", comment.char="#")
row.names(counts)=counts$Geneid
counts = counts[,7:ncol(counts)]
head(counts)
dim(counts)

# look at read count sums 
dim(counts)
tots = apply(counts, 2, sum)
hist(tots)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))

# remove genes with mean read count below cutoff
cut = 3 # this if the gene count cutoff
cc = counts
means = apply(cc,1,mean)
# number of genes that have a mean below the cutoff
table(means>cut)
counts = cc[means>cut,]

# BUILD TRAIT DATA TABLE
# table with rownames equivalent to the column names in the counts data frame
# can have as many columns as sample attributes

# upload and reformat sample info
# rows as sample names and traits as columns

# this is the metadata csv
tdat = read.csv("Metadata_symbiont.csv", header = T)[c('num','treat')]
tdat = tdat %>% 
  mutate(num = sub("$", ".bam", num))
rownames(tdat) = tdat$num
head(tdat)

# subset for the samples with RNAseq data
coldata = tdat %>% 
  filter(num %in% colnames(counts))
dim(coldata)

# select columns in count data
counts = counts %>% 
  select(tdat$num)

# check the rownames and colnames match up
sum(rownames(coldata)==colnames(counts))==ncol(counts)

# make treatment factors
# note: make sure "control" is the 1st level in the treatment factor
str(coldata)
coldata$treat = as.factor(coldata$treat)
coldata$treat = relevel(coldata$treat, "D")
str(coldata)
levels(coldata$treat)

#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
# set up input matrix for DESeq
ddsHTSeq=DESeqDataSetFromMatrix(counts,
                                colData = coldata,
                                # this formula is just based on 1 treatment which
                                # is the environmental characteristic
                                design = formula(~treat))

# run DESeq
dds = DESeq(ddsHTSeq)

# log 2 fold change between these the control and treatment environmental parameter
resultsNames(dds)
# extract desired comparison
res = results(dds,name="treat_C_vs_D")

# retreiving summary of results
# default p-val = 0.1
res01 = results(dds, alpha = 0.1, name="treat_C_vs_D")
summary(res01)

# How many adjusted p-values were less than 0.01?
sum(res01$padj < 0.1, na.rm=TRUE)

# volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-7,7),
                ylim = c(0,15),
                legendPosition = 'bottom')

# extracting log-fold change
# LFC >0 -- up expression
# LFC <0 -- down expression
res
res.df = as.data.frame(row.names(res))
res.df$LFC = res$log2FoldChange
res.df = res.df %>% 
  dplyr::rename("gene"="row.names(res)")

write.csv(res.df,file="DESeq2_LFC_symbiont.csv")

## latitude ####

# upload read counts from htseq
counts = read.delim("feature_counts_out_latitude.txt", comment.char="#")
row.names(counts)=counts$Geneid
counts = counts[,7:ncol(counts)]
head(counts)
dim(counts)

# look at read count sums 
dim(counts)
tots = apply(counts, 2, sum)
hist(tots)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))

# remove genes with mean read count below cutoff
cut = 3 # this if the gene count cutoff
cc = counts
means = apply(cc,1,mean)
# number of genes that have a mean below the cutoff
table(means>cut)
counts = cc[means>cut,]

# BUILD TRAIT DATA TABLE
# table with rownames equivalent to the column names in the counts data frame
# can have as many columns as sample attributes

# upload and reformat sample info
# rows as sample names and traits as columns

# this is the metadata csv
tdat = read.csv("Metadata_latitude.csv", header = T)[c('num','treat','Transplant')]
tdat = tdat %>% 
  mutate(num = sub("$", ".bam", num)) %>% 
  arrange(num)
rownames(tdat) = tdat$num
head(tdat)

# alphabetically order column names in counts dataset
counts = counts %>% 
  select(order(colnames(counts)))

# subset for the samples with RNAseq data
coldata = tdat %>% 
  filter(num %in% colnames(counts))
dim(coldata)

# check the rownames and colnames match up
sum(rownames(coldata)==colnames(counts))==ncol(counts)

# make treatment factors
# note: make sure "control" is the 1st level in the treatment factor
str(coldata)
coldata$treat = as.factor(coldata$treat)
coldata$treat = relevel(coldata$treat, "K")
str(coldata)
levels(coldata$treat)

#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
# set up input matrix for DESeq
ddsHTSeq=DESeqDataSetFromMatrix(counts,
                                colData = coldata,
                                # this formula is just based on 1 treatment which
                                # is the environmental characteristic
                                design = formula(~treat+Transplant))

# run DESeq
dds = DESeq(ddsHTSeq)

# log 2 fold change between these the control and treatment environmental parameter
resultsNames(dds)
# extract desired comparison
res = results(dds,name="treat_O_vs_K")

# retreiving summary of results
# default p-val = 0.1
res01 = results(dds, alpha = 0.1, name="treat_O_vs_K")
summary(res01)

# How many adjusted p-values were less than 0.01?
sum(res01$padj < 0.1, na.rm=TRUE)

# volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-7,7),
                ylim = c(0,15),
                legendPosition = 'bottom')

# extracting log-fold change
# LFC >0 -- up expression
# LFC <0 -- down expression
res
res.df = as.data.frame(row.names(res))
res.df$LFC = res$log2FoldChange
res.df = res.df %>% 
  dplyr::rename("gene"="row.names(res)")

write.csv(res.df,file="DESeq2_LFC_latitude.csv")
