# set working directory
setwd("~/Desktop/analysis")

# load packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(vegan)

## heat stress ####

# list of bam files
bams=read.table("bams_heat_stress.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_heat_stress.ibsMat"))
dimnames(ma)=list(bams,bams)

# original PCA
# grouping for PCA
sites = read.csv("SraRunTable_heat_stress_updated.csv")
sites = left_join(data.frame(Run=bams), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$pool)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()


# removing clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

# pruning clonal replicates (if necessary)
cutclones=0.2
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
  pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
ma1[1:3,1:3]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8)

# use to subset bams to keep in tacc
write.table(paste(pruned,".bam",sep=""),file="bams_noclones_heat_stress",quote=F, col.names=F, row.names=F)

# grouping for PCA
sites = read.csv("SraRunTable_heat_stress_updated.csv")
sites = left_join(data.frame(Run=pruned), # Reorder data frame
                       sites,
                       by = "Run")
populations = as.vector(sites$pool)

# unconstrained ordination (PCoA)
pp1=capscale(ma1~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()

## heat stress [clone check] ####

# list of bam files
bams=read.table("bam_clone_heat_stress.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_clone_heat_stress.ibsMat"))
dimnames(ma)=list(bams,bams)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# checking clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

## co2: dobu ####

# list of bam files
bams=read.table("bams_co2_dobu.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_co2_dobu.ibsMat"))
dimnames(ma)=list(bams,bams)

# original PCA
# grouping for PCA
sites = read.csv("SraRunTable_co2.txt")
sites = left_join(data.frame(Run=bams), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$geo_loc_name)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()


# removing clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

# pruning clonal replicates (if necessary)
cutclones=0.13
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
  pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
ma1[1:3,1:3]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8)

# use to subset bams to keep in tacc
write.table(paste(pruned,".bam",sep=""),file="bams_noclones_co2_dobu",quote=F, col.names=F, row.names=F)

# grouping for PCA
sites = read.csv("SraRunTable_co2.txt")
sites = left_join(data.frame(Run=pruned), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$geo_loc_name)

# unconstrained ordination (PCoA)
pp1=capscale(ma1~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()

## co2: dobu [clone check] ####

# list of bam files
bams=read.table("bam_clone_co2_dobu.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_clone_co2_dobu.ibsMat"))
dimnames(ma)=list(bams,bams)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# checking clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

## co2: upa-upasina ####

# list of bam files
bams=read.table("bams_co2_upaUpasina.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_co2_upaUpasina.ibsMat"))
dimnames(ma)=list(bams,bams)

# original PCA
# grouping for PCA
sites = read.csv("SraRunTable_co2.txt")
sites = left_join(data.frame(Run=bams), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$geo_loc_name)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()


# removing clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

# didn't edit script below

# pruning clonal replicates (if necessary)
cutclones=0.13
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
  pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
ma1[1:3,1:3]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8)

# use to subset bams to keep in tacc
#write.table(paste(pruned,".bam",sep=""),file="bams_noclones_co2_dobu",quote=F, col.names=F, row.names=F)

# grouping for PCA
sites = read.csv("SraRunTable_co2.txt")
sites = left_join(data.frame(Run=pruned), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$geo_loc_name)

# unconstrained ordination (PCoA)
pp1=capscale(ma1~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()

## symbiont ####

# list of bam files
bams=read.table("bams_updated_symbiont.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_symbiont.ibsMat"))
dimnames(ma)=list(bams,bams)

# original PCA
# grouping for PCA
sites = read.csv("SraRunTable_symbiont_updated.csv")
sites = left_join(data.frame(Run=bams), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$clade)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()


# removing clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

# didn't edit script below

# pruning clonal replicates (if necessary)
cutclones=0.13
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
  pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
ma1[1:3,1:3]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8)

# use to subset bams to keep in tacc
#write.table(paste(pruned,".bam",sep=""),file="bams_noclones_co2_dobu",quote=F, col.names=F, row.names=F)

# grouping for PCA
sites = read.csv("SraRunTable_co2.txt")
sites = left_join(data.frame(Run=pruned), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$geo_loc_name)

# unconstrained ordination (PCoA)
pp1=capscale(ma1~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()

## latitude ####

# list of bam files
bams=read.table("bams_latitude.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_latitude.ibsMat"))
dimnames(ma)=list(bams,bams)

# original PCA
# grouping for PCA
sites = read.csv("SraRunTable_latitude_updated.csv")
sites = left_join(data.frame(Run=bams), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$Location)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()


# removing clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)

# didn't edit script below

# pruning clonal replicates (if necessary)
cutclones=0.13
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
  pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
ma1[1:3,1:3]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8)

# use to subset bams to keep in tacc
#write.table(paste(pruned,".bam",sep=""),file="bams_noclones_co2_dobu",quote=F, col.names=F, row.names=F)

# grouping for PCA
sites = read.csv("SraRunTable_co2.txt")
sites = left_join(data.frame(Run=pruned), # Reorder data frame
                  sites,
                  by = "Run")
populations = as.vector(sites$geo_loc_name)

# unconstrained ordination (PCoA)
pp1=capscale(ma1~1)
# simple plot:
plot(pp1)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp1$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp1$CA$eig/sum(pp1$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp1$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

# extracting "scores" table, to plot ; coordinates of the points in the new rotated space
axes2plot=c(1,2) # which PCAs to plot
scores1=data.frame(pp1$CA$u[,axes2plot])

# ggplot
ggplot(scores1,aes(scores1[,1],scores1[,2],color=populations)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores1)[1])+ylab(names(scores1)[2])+coord_equal()

## latitude [clone check] ####

# list of bam files
bams=read.table("bam_clone_latitude.qc")[,1]
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("zz8IBS_clone_latitude.ibsMat"))
dimnames(ma)=list(bams,bams)

# unconstrained ordination (PCoA)
pp1=capscale(ma~1)
# simple plot:
plot(pp1)

# checking clones...

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")
str(hc)

# clustering of samples by IBS (great to detect clones or closely related individuals)
plot(hc,cex=0.8)
