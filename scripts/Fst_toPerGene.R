# set working directory
setwd("~/Desktop/analysis")

# load packages
library(dplyr)
library(tidyverse)

# read gene regions
genes=read.table("Amillepora_genes.tab")
names(genes)=c("contig","start","end","gene")

# add gene annotations, in the same order as genes.txt
gnames=read.table("Amillepora_gnames.tab",sep="\t")
names(gnames)=c("gene","protein")
genes=merge(genes,gnames,by="gene",all.x=T)
genes$protein=as.character(genes$protein)
genes$protein[is.na(genes$protein)]="unknown"
nrow(genes[genes$protein!="unknown",])

# input Fst file from TACC
fst01=read.table("OK_latitude.fst")
names(fst01)=c("contig","pos","a","b")

fst01$contig=as.character(fst01$contig)

# remove zero-only (invariant) bases
fst01[,3:4]=round(fst01[,3:4],3)
ch01=apply(fst01[,3:4],1,sum)
chh=(ch01>0)
table(chh)
# e.g.:
# FALSE   TRUE 
# 25068   335411 

fst01=fst01[chh,]
head(fst01)

# Fst density
jitter=0
f1=fst01$a/(fst01$b+1e-3)+jitter
plot(density(f1),ylim=c(0,250),xlab="Fst",main="",mgp=c(2.3,1,0),bty="n")

# computing weighted Fst per gene
i=1;gfst01=c();ns01=0
pb=txtProgressBar(0,nrow(genes))
for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  #sub=subset(fst01,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
  # to add 1kb to each end of gene coordinates in order to capture untranslated RNA regions:
  sub=subset(fst01,contig==genes$contig[i] & pos>=((genes$start[i])-1000) & pos<=((genes$end[i])+1000))
  if (is.null(sub[1,1]) | sum(sub$b)==0) { 
    gfst01=append(gfst01,NA)
  } else {
    gfst01=append(gfst01,sum(sub$a)/sum(sub$b))
    ns01=ns01+nrow(sub)
  }
}

# density plots of per-gene Fst
plot(density(na.omit(gfst01)))

# check amount of genes with Fst values
sum(!is.na(gfst01))
# check amount of genes without Fst values (i.e., NA)
sum(is.na(gfst01))

# adding results to genes table, saving
genes$fst01=gfst01
head(genes[genes$protein!="unknown",])
save(genes,file="PerGeneFst_latitude.RData")
