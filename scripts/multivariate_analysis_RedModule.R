# set working directory
setwd("~/Desktop/analysis")

# load packages
library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)
library(gplots)
library(cowplot)
library(UpSetR)
library(patchwork)
library(viridis)
library(ggrepel)

# set ggplot theme
theme_set(theme_cowplot())



## LFC: data ####

# heat stress
genes = read.csv("DESeq2_LFC_heat_stress.csv")
heat_stress = genes %>%
  select(gene,LFC)

# co2 dobu
genes = read.csv("DESeq2_LFC_co2_dobu.csv")
co2_dobu = genes %>%
  select(gene,LFC)

# co2 upa-upasina
genes = read.csv("DESeq2_LFC_co2_upaUpasina.csv")
co2_upaUpasina = genes %>%
  select(gene,LFC)

# symbiont
genes = read.csv("DESeq2_LFC_symbiont.csv")
symbiont = genes %>%
  select(gene,LFC)

# latitude
genes = read.csv("DESeq2_LFC_latitude.csv")
latitude <- genes %>%
  select(gene,LFC)

# bind all contrasts into a single dataframe
list = list("heat_stress"=heat_stress,"co2_dobu"=co2_dobu,"co2_upaUpasina"=co2_upaUpasina,"symbiont"=symbiont,"latitude"=latitude)
LFC = list %>% reduce(full_join, by = "gene")
colnames(LFC)[2:6] = as.vector(names(list))

# name rows as gene names
LFC = LFC %>%
  tibble::column_to_rownames('gene')

# remove LFC with Nan
LFC = LFC[complete.cases(LFC),]

LFC.orig = LFC

# standardize
LFC = scale(LFC)



## Red Module: data ####

# You can extract the Red Module hub genes (to analyze overlaps) by simply doing this in R:
load("moduleAssignment.Rdata")
reds = geneModuleMembership %>% 
  select(MMred) %>% 
  rownames_to_column(var = "gene")

nrow(reds)


## LFC (up regulated): venn diagram of ranked genes ####

# rename LFC dataframe
LFC2 = LFC

nrow(LFC)

LFC2 = LFC2 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  inner_join(reds, by="gene") %>%
  column_to_rownames(var = "gene")

nrow(LFC2)

dim(LFC) # original LFC gene count
dim(LFC2) # gene count with Red Module hub genes


# rank genes from original LFC data
for (i in 1:ncol(LFC2)){
  LFC2[,i]=rank(LFC2[,i])
}

LFC2full=LFC2

# 0.875 takes top 12.5% genes
LFC2[LFC2 < (0.875*nrow(LFC2))] = 0 

# keep gene rows with rank sum > 0
LFC22=LFC2[apply(LFC2,1,sum)>0,]

# format for venn diagram
LFC02 = LFC22
LFC02[LFC02 > 0] = 1
vv.real=data.frame(venn(LFC02, show.plot=F)[,1])

# randomized ranking of genes
VV=list()
for(i in 1:5000){ # 5,000 runs
  print(i)
  LFC2r=LFC2full # randomized ranks
  for (i in 1:ncol(LFC2r)){
    LFC2r[,i]=sample(LFC2r[,i])
  }
  
  # 0.875 takes top 12.5% genes
  LFC2r[LFC2r < (0.875*nrow(LFC2r))] = 0
  
  # keep gene rows with rank sum > 0
  LFC22r=LFC2r[apply(LFC2r,1,sum)>0,]
  
  # format for venn diagram
  LFC02r = LFC22r
  LFC02r[LFC02r > 0] = 1
  vv=data.frame(venn(LFC02r, show.plot=F)[,1])
  VV=append(VV,vv)
}
VD=do.call(cbind,VV)

# venn diagram
LFC.venn = LFC22
LFC.venn[LFC.venn > 0] = 1
LFC.venn.r=LFC22r
LFC.venn.r[LFC.venn.r > 0] = 1

# original data
vv=venn(LFC.venn, show.plot=TRUE) # error is for displaying venn diagram, continue

# randomized
vv.r=venn(LFC.venn.r, show.plot=TRUE)



## LFC (up regulated): upset plot ####

# probability of seeing certain number of genes for a contrast combo
ps=c();upper.tail=c()
for(i in 1:nrow(vv.real)) {
  samp=VD[i,]
  test=vv.real[i,1]
  plot(density(samp),xlim=c(mean(samp)-10*sd(samp),mean(samp)+10*sd(samp)))
  abline(v=test)
  if(test>mean(samp)) {
    ps=append(ps,sum(samp>=test)/length(samp))
    upper.tail=c(upper.tail,1)
  } else {
    ps=append(ps,sum(samp<=test)/length(samp))
    upper.tail=c(upper.tail,0)
  }
}

# this is the dataframe to plot
vv.real=cbind(vv.real,pval=ps,upper.tail)

# see if p-values for upset analyses would stand correction
# for multiple comparisons
pval.adj = p.adjust(vv.real$pval, method="BH")

vv.real=cbind(vv.real,pval.adj)

LFC.venn_v2 = LFC.venn
LFC.venn_v2 = LFC.venn_v2 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate_if(is.numeric, as.integer)

## Create the metadata object first
sets = names(LFC.venn_v2[2:7])
metadata = as.data.frame(cbind(sets, names(LFC.venn_v2[2:7])))
valz = c("z","y","x","w","v","u")
metadata = cbind(metadata, valz)
metadata$valz = as.character(metadata$valz)


plt.LFC.up = UpSetR::upset(LFC.venn_v2,
                           nsets=6,
                           sets=names(LFC.venn_v2[2:7]),
                           mainbar.y.max=850,
                           nintersects = NA,
                           order.by = c("freq"),
                           mb.ratio = c(0.6, 0.4),
                           point.size = 3,
                           text.scale = c(2, 2, 0, 0, 2, 0.8),
                           mainbar.y.label="Gene interaction size",
                           set.metadata=list(data=metadata,
                                             plots=list(list(type="matrix_rows",
                                                             column = "valz", colors = c(z = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[1],
                                                                                         y = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[2],
                                                                                         x = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[3],
                                                                                         w = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[4],
                                                                                         v = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[5],
                                                                                         u = "red"
                                                             ),
                                                             alpha = 0.5))))
plt.LFC.up


vv
vv.LFC.up = vv

vv.pvals = as.data.frame(view(vv))
vv.pvals = vv.pvals %>% 
  rownames_to_column(var="combo")

vv.real_v2 = vv.real

vv.real_v2 = vv.real_v2 %>% 
  dplyr::rename("geneCount" = "venn.LFC02..show.plot...F....1.") %>% 
  rownames_to_column(var="combo")

vv.pvals = vv.real_v2 %>% 
  left_join(vv.pvals, by="combo")
view(vv.pvals)
vv.pvals.LFC.up = vv.pvals

write.csv(vv.pvals.LFC.up, file="geneInteractions_pValues_LFC_up_mmRed.csv")



## LFC (down regulated): venn diagram of ranked genes ####

# rename LFC dataframe
LFC3 = LFC

# add Red Module gene hub
LFC3 = LFC3 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  inner_join(reds, by="gene") %>% 
  column_to_rownames(var = "gene")

# rank genes from original LFC data
for (i in 1:ncol(LFC3)){
  LFC3[,i]=rank(LFC3[,i])
}

LFC3full=LFC3

# 0.125 takes bottom 12.5% genes
LFC3[LFC3 > (0.125*nrow(LFC3))] = 0 

# keep gene rows with rank sum > 0
LFC33=LFC3[apply(LFC3,1,sum)>0,]

# format for venn diagram
LFC03 = LFC33
LFC03[LFC03 > 0] = 1
vv.real=data.frame(venn(LFC03, show.plot=F)[,1])

# randomized ranking of genes
VV=list()
for(i in 1:5000){ # 5,000 runs
  print(i)
  LFC3r=LFC3full # randomized ranks
  for (i in 1:ncol(LFC3r)){
    LFC3r[,i]=sample(LFC3r[,i])
  }
  
  # 0.875 takes top 12.5% genes
  LFC3r[LFC3r  > (0.125*nrow(LFC3r))] = 0
  
  # keep gene rows with rank sum > 0
  LFC33r=LFC3r[apply(LFC3r,1,sum)>0,]
  
  # format for venn diagram
  LFC03r = LFC33r
  LFC03r[LFC03r > 0] = 1
  vv=data.frame(venn(LFC03r, show.plot=F)[,1])
  VV=append(VV,vv)
}
VD=do.call(cbind,VV)

# venn diagram
LFC.venn = LFC33
LFC.venn[LFC.venn > 0] = 1
LFC.venn.r=LFC33r
LFC.venn.r[LFC.venn.r > 0] = 1

# original data
vv=venn(LFC.venn, show.plot=TRUE)

# randomized
vv.r=venn(LFC.venn.r, show.plot=TRUE)



## LFC (down regulated): upset plot ####

# probability of seeing certain number of genes for a contrast combo
ps=c();upper.tail=c()
for(i in 1:nrow(vv.real)) {
  samp=VD[i,]
  test=vv.real[i,1]
  plot(density(samp),xlim=c(mean(samp)-10*sd(samp),mean(samp)+10*sd(samp)))
  abline(v=test)
  if(test>mean(samp)) {
    ps=append(ps,sum(samp>=test)/length(samp))
    upper.tail=c(upper.tail,1)
  } else {
    ps=append(ps,sum(samp<=test)/length(samp))
    upper.tail=c(upper.tail,0)
  }
}

# this is the dataframe to plot
vv.real=cbind(vv.real,pval=ps,upper.tail)

# see if p-values for upset analyses would stand correction
# for multiple comparisons
pval.adj = p.adjust(vv.real$pval, method="BH")

vv.real=cbind(vv.real,pval.adj)

LFC.venn_v2 = LFC.venn
LFC.venn_v2 = LFC.venn_v2 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate_if(is.numeric, as.integer)

## Create the metadata object first
sets = names(LFC.venn_v2[2:7])
metadata = as.data.frame(cbind(sets, names(LFC.venn_v2[2:7])))
valz = c("z","y","x","w","v","u")
metadata = cbind(metadata, valz)
metadata$valz = as.character(metadata$valz)


plt.LFC.down = UpSetR::upset(LFC.venn_v2,
                             nsets=6,
                             sets=names(LFC.venn_v2[2:7]),
                             mainbar.y.max=850,
                             mb.ratio = c(0.6, 0.4),
                             nintersects = NA,
                             order.by = c("freq"),
                             point.size = 3,
                             text.scale = c(2, 2, 0, 0, 2, 0.8),
                             mainbar.y.label="Gene interaction size",
                             set.metadata=list(data=metadata,
                                               plots=list(list(type="matrix_rows",
                                                               column = "valz", colors = c(z = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[1],
                                                                                           y = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[2],
                                                                                           x = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[3],
                                                                                           w = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[4],
                                                                                           v = scales::viridis_pal()(length(names(LFC.venn_v2[2:7])))[5],
                                                                                           u = "red"),
                                                               alpha = 0.5))))
plt.LFC.down


vv
vv.LFC.down = vv

vv.pvals = as.data.frame(view(vv))
vv.pvals = vv.pvals %>% 
  rownames_to_column(var="combo")

vv.real_v2 = vv.real

vv.real_v2 = vv.real_v2 %>% 
  dplyr::rename("geneCount" = "venn.LFC03..show.plot...F....1.") %>% 
  rownames_to_column(var="combo")

vv.pvals = vv.real_v2 %>% 
  left_join(vv.pvals, by="combo")
view(vv.pvals)
vv.pvals.LFC.down = vv.pvals

write.csv(vv.pvals.LFC.down, file="geneInteractions_pValues_LFC_down_mmRed.csv")
