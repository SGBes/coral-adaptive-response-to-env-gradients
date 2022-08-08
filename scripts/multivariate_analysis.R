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



## Fst: data ####

# heat stress
load("PerGeneFst_heat_stress.RData")
heat_stress = genes %>%
  select(gene,fst01)

# co2 dobu
load("PerGeneFst_co2_dobu.RData")
co2_dobu = genes %>%
  select(gene,fst01)

# co2 upa-upasina
load("PerGeneFst_co2_upaUpasina.RData")
co2_upaUpasina = genes %>%
  select(gene,fst01)

# symbiont
load("PerGeneFst_symbiont.RData")
symbiont = genes %>%
  select(gene,fst01)

# latitude
load("PerGeneFst_latitude.RData")
latitude = genes %>%
  select(gene,fst01)

# bind all contrasts into a single dataframe
list = list("heat_stress"=heat_stress,"co2_dobu"=co2_dobu,"co2_upaUpasina"=co2_upaUpasina,"symbiont"=symbiont,"latitude"=latitude)
Fst = list %>% reduce(full_join, by = "gene")
colnames(Fst)[2:6] = as.vector(names(list))

# name rows as gene names
Fst = Fst %>%
  tibble::column_to_rownames('gene')

# remove Fst with Nan
Fst = Fst[complete.cases(Fst),]

## Fst: analysis ####

# rename Fst dataframe
fst0 = Fst

# set all Fst < 0 equal to 0
fst0[fst0 < 0] = 0

# keep gene rows with Fst sum > 0
fst00 = fst0[apply(fst0, 1, sum) > 0,]

# check how many genes were removed
dim(fst0) # original
dim(fst00) # with genes removed

# Bray-Curtis -
# only matches species (here, genes) with abundance > 0
fst0.bc = sqrt(vegdist(t(fst00)))

# PCoA -
# like RDA but can use non-Eudclidean dissimilarity indices (e.g., Bray-Curtis distance)
fst0.cap = capscale(fst0.bc ~ 1, comm=t(fst00))

# eigen values
screeplot(fst0.cap)

print(eigenvals(fst0.cap))
#eigvals0.name = as.vector(names(eigenvals(fst0.cap)))
eigvals0.name = c("PC1", "PC2", "PC3", "PC4")
eigvals0.value = as.vector(eigenvals(fst0.cap))

eigvals0.sum = sum(eigvals0.value)


eigvals0 = as.data.frame(cbind(eigvals0.name, eigvals0.value))
eigvals0 = eigvals0 %>%
  dplyr::rename(PC = eigvals0.name, value = eigvals0.value) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(value = (value/eigvals0.sum))

plt.eig0 = ggplot(data=eigvals0, aes(x=PC, y=value)) +
  geom_bar(stat="identity", fill="steelblue") + 
  scale_y_continuous(limit=c(0, 0.3), breaks=seq(0, 0.3, by=0.1)) +
  xlab("Principal component") +
  ylab("Proportion of variance") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))
plt.eig0


# select PCs
#ch = c(3,4)
# plot
#ordiplot(fst0.cap, display="sites", type="t", choices=ch)
#points(scores(fst0.cap, choices=ch)$species*5)


# extract scores with scaling=2, where angles between contrast arrows really reflect correlation
gene.scores00 = data.frame(scores(fst0.cap,scaling=2,display="species",choices=c(1,2)))
contrast.scores00 = data.frame(scores(fst0.cap,scaling=2,display="sites",choices=c(1,2)))

# rescale contrast scores for better plotting
contrast.scores00 = max(abs(gene.scores00))*contrast.scores00/max(abs(contrast.scores00))

print(row.names(contrast.scores00))
# this is for color coding the contrast vectors
Environment=c("Daily heat stress","CO2 seep: Dobu","CO2 seep: Upa-Upasina","Different algal symbiont","Latitudinal gradient")

print(row.names(contrast.scores00))
print(Environment)
contrast.scores00 = contrast.scores00 %>% 
  cbind(Environment) %>% 
  mutate(env = factor(Environment, levels=Environment))

# plot
plt.Fst = ggplot() +
  geom_point(data=gene.scores00, mapping=aes(x=MDS1, y=MDS2), alpha=0.3, color="grey60", size=2) +
  labs(color="Environmental Stressor") +
  geom_segment(data=contrast.scores00,aes(x=0,y=0,xend=MDS1,yend=MDS2,color=env),arrow=arrow(length=unit(0.5,"cm"),angle=20),lineend="round",linejoin="round",size=2) +
  scale_color_viridis(discrete=TRUE) +
  geom_hline(yintercept=0,linetype='dotted') +
  geom_vline(xintercept=0,linetype='dotted') +
  #labs(fill="Gene Count (log10-scale)") +
  xlab("PC1") +
  ylab("PC2") +
  coord_equal() +
  theme(legend.position = c(0.37, 0.15)) +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=14),
        #legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'), #change legend key height
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size
plt.Fst



## Fst: identify significant genes ####

# load gene annotation
Amillepora_gnames=read.delim("Amillepora_gnames.tab", header=FALSE, comment.char="#")
Amillepora_trinotate_annotation_report=read.delim("Amillepora_trinotate_annotation_report.xls", header=FALSE, comment.char="#")

# uniprot annotation
uniprot_annotation = read.table("C:/Users/Sofia/Desktop/analysis/amilsym_to_hsa_uniprot_sprot.2.txt", quote="\"", comment.char="")
uniprot_annotation = uniprot_annotation %>% 
  rename(gene=V1, uniprotAnnot=V2) %>% 
  mutate(gene = str_to_title(gene))

# Fst - PC1
gene.scores00 = gene.scores00 %>% 
  rownames_to_column(var="gene")

gene.scores00.PC1.names = gene.scores00 %>% 
  filter(MDS1 > 0.07) %>% 
  filter(MDS2 < 0) %>% 
  left_join(Amillepora_trinotate_annotation_report, by=c("gene"="V1")) %>% 
  left_join(uniprot_annotation, by=c("gene"))

# see description of gene
view(gene.scores00.PC1.names)

plt.Fst.PC1 = ggplot() +
  geom_point(data=gene.scores00, mapping=aes(x=MDS1, y=MDS2), alpha=0.3, color="grey60", size=2) +
  geom_label_repel(data=gene.scores00.PC1.names, aes(x=MDS1, y=MDS2, label=gene), min.segment.length=0.01) +
  labs(color="Environmental Stressor") +
  geom_segment(data=contrast.scores00, aes(x=0,y=0,xend=MDS1,yend=MDS2,color=env),arrow=arrow(length=unit(0.5,"cm"),angle=20),lineend="round",linejoin="round",size=2) +
  scale_color_viridis(discrete=TRUE) +
  #geom_text(data=gene.scores00, aes(x=MDS1, y=MDS2, label=gene)) +
  geom_hline(yintercept=0,linetype='dotted') +
  geom_vline(xintercept=0,linetype='dotted') +
  #labs(fill="Gene Count (log10-scale)") +
  xlab("PC1") +
  ylab("PC2") +
  coord_equal() +
  theme(legend.position = c(0.37, 0.15)) +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=14),
        #legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'), #change legend key height
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size
plt.Fst.PC1



## Fst: plot ####

plt.Fst.PC1 +
  inset_element(plt.eig0,
                left = 0.7,
                bottom = 0,
                right = 1,
                top = 0.3)



## Fst: ranked genes ####

# rename Fst dataframe
fst2 = Fst

# rank genes from original Fst data
for (i in 1:ncol(fst2)){
  fst2[,i]=rank(fst2[,i])
}
fst0full=fst2

# 0.75 takes top 25% genes
fst2[fst2 < (0.75*nrow(fst2))] = 0 

# keep gene rows with rank sum > 0
fst22=fst2[apply(fst2,1,sum)>0,]

# format for venn diagram
fst02 = fst22
fst02[fst02 > 0] = 1
vv.real=data.frame(venn(fst02, show.plot=F)[,1])

# randomized ranking of genes
VV=list()
for(i in 1:5000){ # 5,000 runs
  print(i)
  fst2r=fst0full # randomized ranks
  for (i in 1:ncol(fst2r)){
    fst2r[,i]=sample(fst2r[,i])
  }
  
  # 0.75 takes top 25% genes
  fst2r[fst2r < (0.75*nrow(fst2r))] = 0
  
  # keep gene rows with Fst sum > 0
  fst22r=fst2r[apply(fst2r,1,sum)>0,]
  
  # format for venn diagram
  fst02r=fst22r
  fst02r[fst02r > 0] = 1
  vv=data.frame(venn(fst02r, show.plot=F)[,1])
  VV=append(VV,vv)
}
VD=do.call(cbind,VV)

# venn diagram
fst.venn = fst22
fst.venn[fst.venn > 0] = 1
fst.venn.r=fst22r
fst.venn.r[fst.venn.r > 0] = 1

# original data
vv=venn(fst.venn, show.plot=TRUE)

# randomized
vv.r=venn(fst.venn.r, show.plot=TRUE)

## Fst: upset plot ####

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

fst.venn_v2 = fst.venn
colnames(fst.venn_v2) = Environment
fst.venn_v2 = fst.venn_v2 %>% 
  rownames_to_column(var="gene") %>% 
  mutate_if(is.numeric, as.integer)

## Create the metadata object first
sets = names(fst.venn_v2[2:6])
metadata = as.data.frame(cbind(sets, Environment))
names(metadata) = c("sets", "Environment")
valz = c("z","y","x","w","v")
metadata = cbind(metadata, valz)
metadata$valz = as.character(metadata$valz)


UpSetR::upset(fst.venn_v2,
              nsets=5,
              sets=Environment,
              mb.ratio = c(0.6, 0.4),
              point.size = 3,
              text.scale = c(2, 2, 0, 0, 2, 0),
              mainbar.y.label="Gene interaction size",
              set.metadata=list(data=metadata,
                                plots=list(list(type="matrix_rows",
                                                column = "valz", colors = c(z = scales::viridis_pal()(length(Environment))[1],
                                                                            y = scales::viridis_pal()(length(Environment))[2],
                                                                            x = scales::viridis_pal()(length(Environment))[3],
                                                                            w = scales::viridis_pal()(length(Environment))[4],
                                                                            v = scales::viridis_pal()(length(Environment))[5]),
                                                alpha = 0.5))))


vv
vv.Fst = vv

vv.pvals = as.data.frame(view(vv))
vv.pvals = vv.pvals %>% 
  rownames_to_column(var="combo")

vv.real_v2 = vv.real

vv.real_v2 = vv.real_v2 %>% 
  dplyr::rename("geneCount" = "venn.fst02..show.plot...F....1.") %>% 
  rownames_to_column(var="combo")

vv.pvals = vv.real_v2 %>% 
  left_join(vv.pvals, by="combo")
view(vv.pvals)
vv.pvals.Fst = vv.pvals

write.csv(vv.pvals.Fst, file="geneInteractions_pValues_Fst.csv")



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



## LFC: plot ####

# use vegan's RDA (=PCA) to do matrix algebra
pc=rda(LFC~1, scale=TRUE)

# eigen values
screeplot(pc)

print(eigenvals(pc))
eigvals1.name = as.vector(names(eigenvals(pc)))
eigvals1.value = as.vector(eigenvals(pc))

eigvals1.sum = sum(eigvals1.value)

eigvals1 = as.data.frame(cbind(eigvals1.name, eigvals1.value))
eigvals1 = eigvals1 %>%
  dplyr::rename(PC = eigvals1.name, value = eigvals1.value) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(value = (value/eigvals1.sum))

print(eigvals1)
eigvals1$PC = c("PC1", "PC2", "PC3", "PC4", "PC5")

plt.eig1 = ggplot(data=eigvals1, aes(x=PC, y=value)) +
  geom_bar(stat="identity", fill="steelblue") + 
  scale_y_continuous(limit=c(0, 0.25), breaks=seq(0, 0.25, by=0.05)) +
  xlab("Principal component") +
  ylab("Proportion of variance") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))
plt.eig1


# select PCs
#ch = c(3,4)
# plot
#ordiplot(pc, display="species", type="t", choices=ch)
#points(scores(pc, choices=ch)$sites*5)


# extract scores with scaling=2, where angles between contrast arrows really reflect correlation
gene.scores11 = data.frame(scores(pc,scaling=2,display="sites",choices=c(1,2)))
contrast.scores11 = data.frame(scores(pc,scaling=2,display="species",choices=c(1,2)))

# rescale contrast scores for better plotting
contrast.scores11 = max(abs(gene.scores11))*contrast.scores11/max(abs(contrast.scores11))

print(row.names(contrast.scores11))
# this is for color coding the contrast vectors
Environment=c("Daily heat stress","CO2 seep: Dobu","CO2 seep: Upa-Upasina","Different algal symbiont","Latitudinal gradient")

print(row.names(contrast.scores11))
print(Environment)
contrast.scores11 = contrast.scores11 %>% 
  cbind(Environment) %>% 
  mutate(env = factor(Environment, levels=Environment))

# plot
plt.LFC = ggplot() +
  geom_point(data=gene.scores11, mapping=aes(x=PC1, y=PC2), alpha=0.3, color="grey60", size=2) +
  labs(color="Environmental Stressor") +
  geom_segment(data=contrast.scores11,aes(x=0,y=0,xend=PC1,yend=PC2,color=env),arrow=arrow(length=unit(0.5,"cm"),angle=20),lineend="round",linejoin="round",size=2) +
  scale_color_viridis(discrete=TRUE) +
  geom_hline(yintercept=0,linetype='dotted') +
  geom_vline(xintercept=0,linetype='dotted') +
  #labs(fill="Gene Count (log10-scale)") +
  xlab("PC1") +
  ylab("PC2") +
  coord_equal() +
  theme(legend.position = c(0.03, 0.85)) +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=14),
        #legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'), #change legend key height
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size
plt.LFC


# combined LFC plot and eigenvalues plot
plt.LFC +
  inset_element(plt.eig1,
                left = 0.0,
                bottom = 0,
                right = 0.3,
                top = 0.3)



## LFC (up regulated): ranked genes ####

# rename LFC dataframe
LFC2 = LFC

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
vv=venn(LFC.venn, show.plot=TRUE)

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
colnames(LFC.venn_v2) = Environment
LFC.venn_v2 = LFC.venn_v2 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate_if(is.numeric, as.integer)

## Create the metadata object first
sets = names(LFC.venn_v2[2:6])
metadata = as.data.frame(cbind(sets, Environment))
names(metadata) = c("sets", "Environment")
valz = c("z","y","x","w","v")
metadata = cbind(metadata, valz)
metadata$valz = as.character(metadata$valz)


plt.LFC.up = UpSetR::upset(LFC.venn_v2,
                           nsets=5,
                           sets=Environment,
                           mainbar.y.max=850,
                           mb.ratio = c(0.6, 0.4),
                           point.size = 3,
                           text.scale = c(2, 2, 0, 0, 2, 0),
                           mainbar.y.label="Gene interaction size",
                           set.metadata=list(data=metadata,
                                             plots=list(list(type="matrix_rows",
                                                             column = "valz", colors = c(z = scales::viridis_pal()(length(Environment))[1],
                                                                                         y = scales::viridis_pal()(length(Environment))[2],
                                                                                         x = scales::viridis_pal()(length(Environment))[3],
                                                                                         w = scales::viridis_pal()(length(Environment))[4],
                                                                                         v = scales::viridis_pal()(length(Environment))[5]),
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

write.csv(vv.pvals.LFC.up, file="geneInteractions_pValues_LFC_up.csv")



## LFC (down regulated): ranked genes ####

# rename LFC dataframe
LFC3 = LFC

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
colnames(LFC.venn_v2) = Environment
LFC.venn_v2 = LFC.venn_v2 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate_if(is.numeric, as.integer)

## Create the metadata object first
sets = names(LFC.venn_v2[2:6])
metadata = as.data.frame(cbind(sets, Environment))
names(metadata) = c("sets", "Environment")
valz = c("z","y","x","w","v")
metadata = cbind(metadata, valz)
metadata$valz = as.character(metadata$valz)


plt.LFC.down = UpSetR::upset(LFC.venn_v2,
                             nsets=5,
                             sets=Environment,
                             mainbar.y.max=850,
                             mb.ratio = c(0.6, 0.4),
                             point.size = 3,
                             text.scale = c(2, 2, 0, 0, 2, 0),
                             mainbar.y.label="Gene interaction size",
                             set.metadata=list(data=metadata,
                                               plots=list(list(type="matrix_rows",
                                                               column = "valz", colors = c(z = scales::viridis_pal()(length(Environment))[1],
                                                                                           y = scales::viridis_pal()(length(Environment))[2],
                                                                                           x = scales::viridis_pal()(length(Environment))[3],
                                                                                           w = scales::viridis_pal()(length(Environment))[4],
                                                                                           v = scales::viridis_pal()(length(Environment))[5]),
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

write.csv(vv.pvals.LFC.down, file="geneInteractions_pValues_LFC_down.csv")



## Fst and LFC: gene overlap ####

vv.Lst.Fst = attr(vv.Fst, "intersections")
vv.Lst.LFC.up = attr(vv.LFC.up, "intersections")
vv.Lst.LFC.down = attr(vv.LFC.down, "intersections")



## Fst vs LFC up-regulated

# match list orders
vv.Lst.LFC.up = vv.Lst.LFC.up[order(match(vv.Lst.LFC.up, vv.Lst.Fst))]

comp.Fst.LFCup = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("comparison", "fst", "shared", "lfc.up", "fst.length", "lfc.up.length"))
comp.shared.Fst.LFCup = setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("comparison", "gene"))

for (i in 1:length(vv.Lst.Fst)){
  print(i)
  fst.temp = unlist(vv.Lst.Fst[[i]], use.names=FALSE)
  lfc.temp = unlist(vv.Lst.LFC.up[[i]], use.names=FALSE)
  
  # comparison
  comp.name = names(vv.Lst.Fst[i])
  
  # genes in both
  gene.union = length(intersect(fst.temp, lfc.temp))
  gene.union.names = intersect(fst.temp, lfc.temp)
  
  if (gene.union > 0){
    
    for (x in 1:length(gene.union.names)){
      
      geneAdd = gene.union.names[x]
      
      geneVectorAdd = data.frame(comp.name, geneAdd)
      colnames(geneVectorAdd) = c("comparison", "gene")
      
      comp.shared.Fst.LFCup = rbind(comp.shared.Fst.LFCup, geneVectorAdd)
      
    }
    
  }
  
  # genes in fst only
  gene.fst = length(setdiff(fst.temp, lfc.temp))
  gene.fst.length = length(fst.temp)
  
  # genes in lfc only
  gene.lfc = length(setdiff(lfc.temp, fst.temp))
  gene.lfc.length = length(lfc.temp)
  
  vectorAdd = data.frame(comp.name, gene.fst, gene.union, gene.lfc, gene.fst.length, gene.lfc.length)
  colnames(vectorAdd) = c("comparison", "fst", "shared", "lfc.up", "fst.length", "lfc.up.length")
  
  comp.Fst.LFCup = rbind(comp.Fst.LFCup, vectorAdd)
  
}



## Fst vs LFC down-regulated

# match list orders
vv.Lst.LFC.down = vv.Lst.LFC.down[order(match(vv.Lst.LFC.down, vv.Lst.Fst))]

comp.Fst.LFCdown = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("comparison", "fst", "shared", "lfc.down", "fst.length", "lfc.down.length"))
comp.shared.Fst.LFCdown = setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("comparison", "gene"))

for (i in 1:length(vv.Lst.Fst)){
  print(i)
  fst.temp = unlist(vv.Lst.Fst[[i]], use.names=FALSE)
  lfc.temp = unlist(vv.Lst.LFC.down[[i]], use.names=FALSE)
  
  # comparison
  comp.name = names(vv.Lst.Fst[i])
  
  # genes in both
  gene.union = length(intersect(fst.temp, lfc.temp))
  gene.union.names = intersect(fst.temp, lfc.temp)
  
  if (gene.union > 0){
    
    for (x in 1:length(gene.union.names)){
      
      geneAdd = gene.union.names[x]
      
      geneVectorAdd = data.frame(comp.name, geneAdd)
      colnames(geneVectorAdd) = c("comparison", "gene")
      
      comp.shared.Fst.LFCdown = rbind(comp.shared.Fst.LFCdown, geneVectorAdd)
      
    }
    
  }
  
  # genes in fst only
  gene.fst = length(setdiff(fst.temp, lfc.temp))
  gene.fst.length = length(fst.temp)
  
  # genes in lfc only
  gene.lfc = length(setdiff(lfc.temp, fst.temp))
  gene.lfc.length = length(lfc.temp)
  
  vectorAdd = data.frame(comp.name, gene.fst, gene.union, gene.lfc, gene.fst.length, gene.lfc.length)
  colnames(vectorAdd) = c("comparison", "fst", "shared", "lfc.down", "fst.length", "lfc.down.length")
  
  comp.Fst.LFCdown = rbind(comp.Fst.LFCdown, vectorAdd)
  
}



## Fst and LFC: plot of shared genes ####

comparisonDF = comp.Fst.LFCdown %>% 
  left_join(comp.Fst.LFCup, by="comparison", suffix=c(".down", ".up")) %>% 
  select(comparison, fst.length.down, lfc.down.length, lfc.up.length, shared.down, shared.up) %>% 
  dplyr::rename("Fst"="fst.length.down",
                "Positive log-fold change" = "lfc.up.length",
                "Negative log-fold change" = "lfc.down.length",
                "Fst : Positive log-fold change" = "shared.up",
                "Fst : Negative log-fold change" = "shared.down") %>% 
  pivot_longer(-comparison, names_to = "value", values_to = "count")

comparisonDF$value = factor(comparisonDF$value, levels=c("Negative log-fold change",
                                                         "Fst : Negative log-fold change",
                                                         "Fst",
                                                         "Fst : Positive log-fold change",
                                                         "Positive log-fold change"))

comparisonDF$comparison = factor(comparisonDF$comparison, levels=(sort(unique(comparisonDF$comparison), decreasing=TRUE)))


comp_Fst_LFC = ggplot(data=comparisonDF, aes(x=value, y=comparison)) +
  geom_tile(color="black", fill="white") +
  geom_text(aes(label=count), color = "black", size = 3) +
  ylab("Comparison") +
  xlab("Dataset") +
  labs(fill="Shared genes") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=14),
        #legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(0.6, 'cm'), #change legend key height
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size
comp_Fst_LFC



## LFC (up regulated) and Fst: significant shared genes ####

# randomize Fst values in dataset
VV = list()
VV_v2 = list()
for(i in 1:1000){ # 1,000 runs
  
  print(i)
  
  fst90r=fst0 # randomized Fst values
  for (i in 1:ncol(fst90r)){
    fst90r[,i]=sample(fst90r[,i])
  }
  
  # rank genes from Fst data
  for (i in 1:ncol(fst90r)){
    fst90r[,i]=rank(fst90r[,i])
  }
  
  # 0.75 takes top 25% genes
  fst90r[fst90r < (0.75*nrow(fst90r))] = 0 
  
  # keep gene rows with rank sum > 0
  fst90r=fst90r[apply(fst90r,1,sum)>0,]
  
  # format for venn diagram
  fst99 = fst90r
  fst99[fst99 > 0] = 1
  vv = data.frame(venn(fst99, show.plot=F)[,1])
  vv_temp = venn(fst99, show.plot=F)
  vv_v2 = list(attr(vv_temp, "intersections"))
  
  VV = append(VV,vv)
  
  
  VV_v2 = append(VV_v2, vv_v2)
  
}
VD = do.call(cbind,VV)
vv.r.fst = VV

vv_v2.r.fst = VV_v2


# randomize LFC values in dataset - subset for up-regulated
VV = list()
VV_v2 = list()
for(i in 1:1000){ # 1,000 runs
  
  print(i)
  
  lfc90r=LFC # randomized LFC values
  for (i in 1:ncol(lfc90r)){
    lfc90r[,i]=sample(lfc90r[,i])
  }
  
  # rank genes from LFC data
  for (i in 1:ncol(lfc90r)){
    lfc90r[,i]=rank(lfc90r[,i])
  }
  
  # 0.875 takes top 12.5% genes
  lfc90r[lfc90r < (0.875*nrow(lfc90r))] = 0 
  
  # keep gene rows with rank sum > 0
  lfc90r=lfc90r[apply(lfc90r,1,sum)>0,]
  
  # format for venn diagram
  lfc99 = lfc90r
  lfc99[lfc99 > 0] = 1
  vv = data.frame(venn(lfc99, show.plot=F)[,1])
  vv_temp = venn(fst99, show.plot=F)
  vv_v2 = list(attr(vv_temp, "intersections"))
  
  VV = append(VV,vv)
  VV_v2 = append(VV_v2, vv_v2)
  
}
VC=do.call(cbind,VV)
vv.r.lfc.up = VV

vv_v2.r.lfc.up = VV_v2


comp.shared.Fst.LFCup.r = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("comparison", "gene", "cycle"))

# shared genes from randomization
for (i in 1:length(vv_v2.r.fst)){
  
  print(i)
  
  vv.fst = vv_v2.r.fst[[i]]
  vv.lfc = vv_v2.r.lfc.up[[i]]
  
  toDelete = setdiff(names(vv.fst), names(vv.lfc))
  
  if (length(toDelete) > 0){
    
    for (z in 1:length(toDelete)){
      
      if (toDelete[z] %in% names(vv.fst)){
        
        vv.fst = vv.fst[names(vv.fst) != toDelete[z]]
        
      }
      
      if (toDelete[z] %in% names(vv.lfc)){
        
        vv.lfc = vv.lfc[names(vv.lfc) != toDelete[z]]
        
      }
      
    }
    
  }
  
  for (x in 1:length(vv.fst)){
    
    comp.name = names(vv.fst[x])
    
    fst.temp = unlist(vv.fst[[x]], use.names=FALSE)
    lfc.temp = unlist(vv.lfc[[x]], use.names=FALSE)
    
    gene.union.names = intersect(fst.temp, lfc.temp)
    
    if (length(gene.union.names) > 0){
      
      for (y in 1:length(gene.union.names)){
        
        geneVectorAdd = data.frame(comp.name, gene.union.names[y], i)
        colnames(geneVectorAdd) = c("comparison", "gene", "cycle")
        
        comp.shared.Fst.LFCup.r = rbind(comp.shared.Fst.LFCup.r, geneVectorAdd)
        
      }
    }
  }
}

# randomized shared genes dataframe
comp.shared.Fst.LFCup.r_v2 = comp.shared.Fst.LFCup.r %>% 
  group_by(comparison, cycle) %>% 
  tally()


print(names(attr(vv.Fst, "intersections")))
venn.names = c("none",
               "latitude",
               "symbiont",
               "symbiont:latitude",
               "co2_upaUpasina",
               "co2_upaUpasina:latitude",
               "co2_upaUpasina:symbiont",
               "co2_upaUpasina:symbiont:latitude",
               "co2_dobu",
               "co2_dobu:latitude",
               "co2_dobu:symbiont",
               "co2_dobu:symbiont:latitude",
               "co2_dobu:co2_upaUpasina",
               "co2_dobu:co2_upaUpasina:latitude",
               "co2_dobu:co2_upaUpasina:symbiont",
               "co2_dobu:co2_upaUpasina:symbiont:latitude",
               "heat_stress",
               "heat_stress:latitude",
               "heat_stress:symbiont",
               "heat_stress:symbiont:latitude",
               "heat_stress:co2_upaUpasina",
               "heat_stress:co2_upaUpasina:latitude",                 
               "heat_stress:co2_upaUpasina:symbiont",
               "heat_stress:co2_upaUpasina:symbiont:latitude",
               "heat_stress:co2_dobu",
               "heat_stress:co2_dobu:latitude",                        
               "heat_stress:co2_dobu:symbiont",
               "heat_stress:co2_dobu:symbiont:latitude",
               "heat_stress:co2_dobu:co2_upaUpasina",
               "heat_stress:co2_dobu:co2_upaUpasina:latitude",
               "heat_stress:co2_dobu:co2_upaUpasina:symbiont",
               "heat_stress:co2_dobu:co2_upaUpasina:symbiont:latitude")

venn.values = as.vector(vv.pvals.Fst$combo)

venn.info = as.data.frame(cbind(venn.names, venn.values))
venn.info = venn.info %>% filter(venn.values != "00000")


comp.Fst.LFCup_v2 = comp.Fst.LFCup %>% 
  select(comparison, shared) %>% 
  left_join(venn.info, by=c("comparison"="venn.names"))

comp.shared.Fst.LFCup.r_v3 = comp.shared.Fst.LFCup.r_v2 %>% 
  pivot_wider(names_from = "cycle", values_from="n") %>% 
  replace(is.na(.), 0)

comp.shared.Fst.LFCup.r_v3 = subset(comp.shared.Fst.LFCup.r_v3, (comp.shared.Fst.LFCup.r_v3$comparison %in% comp.Fst.LFCup_v2$comparison))
comp.shared.Fst.LFCup.r_v3 = comp.shared.Fst.LFCup.r_v3[order(match(comp.shared.Fst.LFCup.r_v3$comparison, comp.Fst.LFCup_v2$comparison)), ]

comp.Fst.LFCup_v2 = subset(comp.Fst.LFCup_v2, (comp.Fst.LFCup_v2$comparison %in% comp.shared.Fst.LFCup.r_v3$comparison))

comp.shared.Fst.LFCup.r_v3 = comp.shared.Fst.LFCup.r_v3 %>% 
  column_to_rownames(var="comparison")


# probability of seeing certain shared number of genes for a contrast combo
ps=c()
for(i in 1:nrow(comp.Fst.LFCup_v2)) {
  samp=as.numeric(comp.shared.Fst.LFCup.r_v3[i,])
  test=as.numeric(comp.Fst.LFCup_v2[i,2])
  #  if(test>mean(samp)) {
  ps=append(ps,sum(samp>=test)/length(samp))
  #  } else {
  #    ps=append(ps,sum(samp<=test)/length(samp))
  #  }
}

# this is the dataframe to plot
comp.Fst.LFCup_v2=cbind(comp.Fst.LFCup_v2,pval=ps)
view(comp.Fst.LFCup_v2)

write.csv(comp.Fst.LFCup_v2, file="sharedGenes_pValues_Fst_LFC_upregulated.csv")



## LFC (down regulated) and Fst: significant shared genes ####

# randomize Fst values in dataset
VV = list()
VV_v2 = list()
for(i in 1:1000){ # 1,000 runs
  
  print(i)
  
  fst90r=fst0 # randomized Fst values
  for (i in 1:ncol(fst90r)){
    fst90r[,i]=sample(fst90r[,i])
  }
  
  # rank genes from Fst data
  for (i in 1:ncol(fst90r)){
    fst90r[,i]=rank(fst90r[,i])
  }
  
  # 0.75 takes top 25% genes
  fst90r[fst90r < (0.75*nrow(fst90r))] = 0 
  
  # keep gene rows with rank sum > 0
  fst90r=fst90r[apply(fst90r,1,sum)>0,]
  
  # format for venn diagram
  fst99 = fst90r
  fst99[fst99 > 0] = 1
  vv = data.frame(venn(fst99, show.plot=F)[,1])
  vv_temp = venn(fst99, show.plot=F)
  vv_v2 = list(attr(vv_temp, "intersections"))
  
  VV = append(VV,vv)
  
  
  VV_v2 = append(VV_v2, vv_v2)
  
}
VD_x3 = do.call(cbind,VV)
vv.r.fst_x3 = VV

vv_v2.r.fst_x3 = VV_v2


# randomize LFC values in dataset - subset for down-regulated
VV = list()
VV_v2 = list()
for(i in 1:1000){ # 1,000 runs
  
  print(i)
  
  lfc90r=LFC # randomized LFC values
  for (i in 1:ncol(lfc90r)){
    lfc90r[,i]=sample(lfc90r[,i])
  }
  
  # rank genes from LFC data
  for (i in 1:ncol(lfc90r)){
    lfc90r[,i]=rank(lfc90r[,i])
  }
  
  # 0.125 takes bottom 12.5% genes
  lfc90r[lfc90r > (0.125*nrow(lfc90r))] = 0
  
  # keep gene rows with rank sum > 0
  lfc90r=lfc90r[apply(lfc90r,1,sum)>0,]
  
  # format for venn diagram
  lfc99 = lfc90r
  lfc99[lfc99 > 0] = 1
  vv = data.frame(venn(lfc99, show.plot=F)[,1])
  vv_temp = venn(fst99, show.plot=F)
  vv_v2 = list(attr(vv_temp, "intersections"))
  
  VV = append(VV,vv)
  VV_v2 = append(VV_v2, vv_v2)
  
}
VC_x3=do.call(cbind,VV)
vv.r.lfc.down_x3 = VV

vv_v2.r.lfc.down_x3 = VV_v2


comp.shared.Fst.LFCdown.r = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("comparison", "gene", "cycle"))

# shared genes from randomization
for (i in 1:length(vv_v2.r.fst_x3)){
  
  print(i)
  
  vv.fst = vv_v2.r.fst_x3[[i]]
  vv.lfc = vv_v2.r.lfc.down_x3[[i]]
  
  toDelete = setdiff(names(vv.fst), names(vv.lfc))
  
  if (length(toDelete) > 0){
    
    for (z in 1:length(toDelete)){
      
      if (toDelete[z] %in% names(vv.fst)){
        
        vv.fst = vv.fst[names(vv.fst) != toDelete[z]]
        
      }
      
      if (toDelete[z] %in% names(vv.lfc)){
        
        vv.lfc = vv.lfc[names(vv.lfc) != toDelete[z]]
        
      }
      
    }
    
  }
  
  for (x in 1:length(vv.fst)){
    
    comp.name = names(vv.fst[x])
    
    fst.temp = unlist(vv.fst[[x]], use.names=FALSE)
    lfc.temp = unlist(vv.lfc[[x]], use.names=FALSE)
    
    gene.union.names = intersect(fst.temp, lfc.temp)
    
    if (length(gene.union.names) > 0){
      
      for (y in 1:length(gene.union.names)){
        
        geneVectorAdd = data.frame(comp.name, gene.union.names[y], i)
        colnames(geneVectorAdd) = c("comparison", "gene", "cycle")
        
        comp.shared.Fst.LFCdown.r = rbind(comp.shared.Fst.LFCdown.r, geneVectorAdd)
        
      }
    }
  }
}

# randomized shared genes dataframe
comp.shared.Fst.LFCdown.r_v2 = comp.shared.Fst.LFCdown.r %>% 
  group_by(comparison, cycle) %>% 
  tally()


print(names(attr(vv.Fst, "intersections")))
venn.names = c("none",
               "latitude",
               "symbiont",
               "symbiont:latitude",
               "co2_upaUpasina",
               "co2_upaUpasina:latitude",
               "co2_upaUpasina:symbiont",
               "co2_upaUpasina:symbiont:latitude",
               "co2_dobu",
               "co2_dobu:latitude",
               "co2_dobu:symbiont",
               "co2_dobu:symbiont:latitude",
               "co2_dobu:co2_upaUpasina",
               "co2_dobu:co2_upaUpasina:latitude",
               "co2_dobu:co2_upaUpasina:symbiont",
               "co2_dobu:co2_upaUpasina:symbiont:latitude",
               "heat_stress",
               "heat_stress:latitude",
               "heat_stress:symbiont",
               "heat_stress:symbiont:latitude",
               "heat_stress:co2_upaUpasina",
               "heat_stress:co2_upaUpasina:latitude",                 
               "heat_stress:co2_upaUpasina:symbiont",
               "heat_stress:co2_upaUpasina:symbiont:latitude",
               "heat_stress:co2_dobu",
               "heat_stress:co2_dobu:latitude",                        
               "heat_stress:co2_dobu:symbiont",
               "heat_stress:co2_dobu:symbiont:latitude",
               "heat_stress:co2_dobu:co2_upaUpasina",
               "heat_stress:co2_dobu:co2_upaUpasina:latitude",
               "heat_stress:co2_dobu:co2_upaUpasina:symbiont",
               "heat_stress:co2_dobu:co2_upaUpasina:symbiont:latitude")

venn.values = as.vector(vv.pvals.Fst$combo)

venn.info = as.data.frame(cbind(venn.names, venn.values))
venn.info = venn.info %>% filter(venn.values != "00000")


comp.Fst.LFCdown_v2 = comp.Fst.LFCdown %>% 
  select(comparison, shared) %>% 
  left_join(venn.info, by=c("comparison"="venn.names"))

comp.shared.Fst.LFCdown.r_v3 = comp.shared.Fst.LFCdown.r_v2 %>% 
  pivot_wider(names_from = "cycle", values_from="n") %>% 
  replace(is.na(.), 0)

comp.shared.Fst.LFCdown.r_v3 = subset(comp.shared.Fst.LFCdown.r_v3, (comp.shared.Fst.LFCdown.r_v3$comparison %in% comp.Fst.LFCdown_v2$comparison))
comp.shared.Fst.LFCdown.r_v3 = comp.shared.Fst.LFCdown.r_v3[order(match(comp.shared.Fst.LFCdown.r_v3$comparison, comp.Fst.LFCdown_v2$comparison)), ]

comp.Fst.LFCdown_v2 = subset(comp.Fst.LFCdown_v2, (comp.Fst.LFCdown_v2$comparison %in% comp.shared.Fst.LFCdown.r_v3$comparison))

comp.shared.Fst.LFCdown.r_v3 = comp.shared.Fst.LFCdown.r_v3 %>% 
  column_to_rownames(var="comparison")


# probability of seeing certain shared number of genes for a contrast combo
ps=c()
for(i in 1:nrow(comp.Fst.LFCdown_v2)) {
  samp=as.numeric(comp.shared.Fst.LFCdown.r_v3[i,])
  test=as.numeric(comp.Fst.LFCdown_v2[i,2])
  #  if(test>mean(samp)) {
  ps=append(ps,sum(samp>=test)/length(samp))
  #  } else {
  #    ps=append(ps,sum(samp<=test)/length(samp))
  #  }
}

# this is the dataframe to plot
comp.Fst.LFCdown_v2=cbind(comp.Fst.LFCdown_v2,pval=ps)
view(comp.Fst.LFCdown_v2)

write.csv(comp.Fst.LFCdown_v2, file="sharedGenes_pValues_Fst_LFC_downregulated")
