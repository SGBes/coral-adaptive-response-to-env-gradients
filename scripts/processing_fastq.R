# set working directory
setwd("~/Desktop/analysis")

# load packages
library(dplyr)
library(tidyverse)

## heat stress ####

# SRA metadata
SraRunTable_heat_stress = read.csv("SraRunTable_heat_stress.txt")

heatStress_SRA = SraRunTable_heat_stress %>% 
  mutate(geo_loc_name = gsub("\\\\", "", geo_loc_name)) %>% # remove escape character
  filter(LibrarySource == "TRANSCRIPTOMIC") %>% # keep transcriptomic biosamples
  filter(Host == "Acropora hyacinthus (HE)") %>% # keep biosamples with same coral species within different reef pools
  filter(geo_loc_name %in% c("American Samoa: Pool 400, Ofu", "American Samoa: Pool 300, Ofu")) # keep biosamples collected in reef pools

heatStress_SRA_names = heatStress_SRA %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".fastq", Run)) # add ".fastq" to end of string

# use this to select for fastq files to keep in bash directory
write.csv(heatStress_SRA_names, file="SRAnames_heat_stress.csv", row.names=FALSE)
# remove column name afterwards

# create population lists
write.csv(SraRunTable_heat_stress, file="SraRunTable_heat_stress_updated.csv")
# use supplemental information to add pool metadata

SraRunTable_heat_stress_updated = read.csv("SraRunTable_heat_stress_updated.csv")

# metadata
meta_heat_stress = SraRunTable_heat_stress_updated %>% 
  select(Run, pool) %>% 
  filter(pool %in% c("MV", "HV")) %>% 
  dplyr::rename(num=Run,treat=pool)

write.csv(meta_heat_stress, file="Metadata_heat_stress.csv")

# moderately variable pool
PoolMV = SraRunTable_heat_stress_updated %>% 
  filter(pool == "MV")

PoolMV_names = PoolMV %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

# highly variable pool
PoolHV = SraRunTable_heat_stress_updated %>% 
  filter(pool == "HV")

PoolHV_names = PoolHV %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

write.csv(PoolMV_names, file="SRAnames_heat_stress_pop_MV.csv")
write.csv(PoolHV_names, file="SRAnames_heat_stress_pop_HV.csv")
# remove column and row names afterwards

## co2 ####

# SRA metadata
SraRunTable_co2 = read.csv("SraRunTable_co2.txt")

# Site: Dobu
co2_dobu_SRA = SraRunTable_co2 %>% 
  mutate(geo_loc_name = gsub("\\\\", "", geo_loc_name)) %>% # remove escape character
  filter(geo_loc_name %in% c("Papua New Guinea: Milne Bay Province, Dobu-Control", "Papua New Guinea: Milne Bay Province, Dobu-Seep")) # keep biosamples at Dobu site

co2_dobu_SRA_names = co2_dobu_SRA %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".fastq", Run)) # add ".fastq" to end of string

write.csv(co2_dobu_SRA_names, file="SRAnames_co2_dobu.csv", row.names=FALSE)
# remove column name afterwards

# metadata
meta_co2_dobu = co2_dobu_SRA %>% 
  select(Run, geo_loc_name) %>% 
  mutate(geo_loc_name = if_else(geo_loc_name=="Papua New Guinea: Milne Bay Province, Dobu-Control", "Control", "Seep")) %>% 
  dplyr::rename(num=Run, treat=geo_loc_name)

write.csv(meta_co2_dobu, file="Metadata_co2_dobu.csv")

# create population lists

# co2 dobu seep
co2_dobu_seep = co2_dobu_SRA %>% 
  filter(geo_loc_name == "Papua New Guinea: Milne Bay Province, Dobu-Seep")

co2_dobu_seep_names = co2_dobu_seep %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

# co2 dobu control
co2_dobu_control = co2_dobu_SRA %>% 
  filter(geo_loc_name == "Papua New Guinea: Milne Bay Province, Dobu-Control")

co2_dobu_control_names = co2_dobu_control %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

write.csv(co2_dobu_seep_names, file="SRAnames_co2_dobu_pop_seep.csv")
write.csv(co2_dobu_control_names, file="SRAnames_co2_dobu_pop_control.csv")
# remove column and row names afterwards

# Site: Upa-Upasina
co2_upaUpasina_SRA = SraRunTable_co2 %>% 
  mutate(geo_loc_name = gsub("\\\\", "", geo_loc_name)) %>% # remove escape character
  filter(geo_loc_name %in% c("Papua New Guinea: Milne Bay Province, Upa-Upasina-seep", "Papua New Guinea: Milne Bay Province, Upa-Upasina-control")) # keep biosamples at Dobu site

co2_upaUpasina_SRA_names = co2_upaUpasina_SRA %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".fastq", Run)) # add ".fastq" to end of string

write.csv(co2_upaUpasina_SRA_names, file="SRAnames_co2_upaUpasina.csv", row.names=FALSE)
# remove column name afterwards

# metadata
meta_co2_upaUpasina = co2_upaUpasina_SRA %>% 
  select(Run, geo_loc_name) %>% 
  mutate(geo_loc_name = if_else(geo_loc_name=="Papua New Guinea: Milne Bay Province, Upa-Upasina-control", "Control", "Seep")) %>% 
  dplyr::rename(num=Run, treat=geo_loc_name)

write.csv(meta_co2_upaUpasina, file="Metadata_co2_upaUpasina.csv")

# create population lists

# co2 upa-upasina seep
co2_upaUpasina_seep = co2_upaUpasina_SRA %>% 
  filter(geo_loc_name == "Papua New Guinea: Milne Bay Province, Upa-Upasina-seep")

co2_upaUpasina_seep_names = co2_upaUpasina_seep %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

# co2 upa-upasina control
co2_upaUpasina_control = co2_upaUpasina_SRA %>% 
  filter(geo_loc_name == "Papua New Guinea: Milne Bay Province, Upa-Upasina-control")

co2_upaUpasina_control_names = co2_upaUpasina_control %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

write.csv(co2_upaUpasina_seep_names, file="SRAnames_co2_upaUpasina_pop_seep.csv")
write.csv(co2_upaUpasina_control_names, file="SRAnames_co2_upaUpasina_pop_control.csv")
# remove column and row names afterwards

## symbiont ####

# SRA metadata
SraRunTable_symbiont = read.csv("SraRunTable_symbiont.txt")

# for concatenating genotyping replicates with samtools
symbiont_SRA = SraRunTable_symbiont %>% 
  select(Run, Sample.Name) %>% 
  separate(Sample.Name, into = c('Sample.Name', 'Sample.Clone'), sep = -1, convert = TRUE) %>% 
  pivot_wider(names_from=Sample.Clone, values_from=Run)

symbiont_SRA_concat2 = symbiont_SRA %>% 
  filter_all(any_vars(is.na(.))) %>% 
  mutate(c = coalesce(c,b)) %>% # check which columns are missing sample names
  select(-c(b)) %>% 
  select(a, c, Sample.Name)

symbiont_SRA_concat3 = symbiont_SRA %>% 
  drop_na() %>% 
  select(a, b, c, Sample.Name)

write.csv(symbiont_SRA_concat2, file="symbiont_forConcatenate_n2.csv", col.names=FALSE, row.names=FALSE)
write.csv(symbiont_SRA_concat3, file="symbiont_forConcatenate_n3.csv", col.names=FALSE, row.names=FALSE)

# need updated list of bams for tacc to exclude samples with too few reads
# these samples were excluded from original study

symbiont_SRA_updated = symbiont_SRA %>% 
  mutate(clade = "C") %>% 
  arrange(Sample.Name)
# last 4 Wilkie samples are clade D
symbiont_SRA_updated[(nrow(symbiont_SRA_updated)-3):nrow(symbiont_SRA),5] = "D"
symbiont_SRA_updated = symbiont_SRA_updated %>% 
  filter(!Sample.Name %in% c("o8","o4","oM1")) %>% 
  filter(grepl("w",Sample.Name)) %>% 
  select(Sample.Name,clade) %>% 
  dplyr::rename(Run=Sample.Name)

symbiont_SRA_updated_names = symbiont_SRA_updated %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

write.csv(symbiont_SRA_updated_names, file="SRAnames_symbiont_updated.csv")
write.csv(symbiont_SRA_updated, file="SraRunTable_symbiont_updated.csv")

# metadata
symbiont_SRA_updated = symbiont_SRA_updated %>% 
  dplyr::rename(num=Run, treat=clade)

write.csv(symbiont_SRA_updated, file="Metadata_symbiont.csv")

# create population lists

# clade C
symbiont_c = symbiont_SRA_updated %>% 
  filter(clade == "C")

symbiont_c_names = symbiont_c %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

# clade D
symbiont_d = symbiont_SRA_updated %>% 
  filter(clade == "D")

symbiont_d_names = symbiont_d %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

write.csv(symbiont_c_names, file="SRAnames_symbiont_pop_c.csv")
write.csv(symbiont_d_names, file="SRAnames_symbiont_pop_d.csv")


## latitude ####

# SRA metadata
SraRunTable_latitude = read.csv("SraRunTable_latitude.txt")

latitude_SRA = SraRunTable_latitude %>% 
  filter(LibrarySource == "TRANSCRIPTOMIC")

latitude_SRA_names = latitude_SRA %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".fastq", Run)) # add ".fastq" to end of string

write.csv(latitude_SRA_names, file="SRAnames_latitude.csv", row.names=FALSE)
# remove column name afterwards

# for concatenating genotyping replicates with samtools
latitude_SRA_concat = latitude_SRA %>% 
  select(Run, Genotype) %>% 
  group_by(Genotype) %>% 
  summarize(concat = str_c(Run, collapse="; ")) %>% 
  separate(concat, c("Run1", "Run2")) %>% 
  select(Run1, Run2, Genotype)

latitude_SRA_concat1 = latitude_SRA_concat %>% 
  filter_all(any_vars(is.na(.))) %>% 
  mutate(Run1 = coalesce(Run1,Run2)) %>% # check which columns are missing sample names
  select(Run1, Genotype)

latitude_SRA_concat2 = latitude_SRA_concat %>% 
  drop_na() %>% 
  select(Run1, Run2, Genotype)

write.csv(latitude_SRA_concat1, file="latitude_forConcatenate_n1.csv", col.names=FALSE, row.names=FALSE)
write.csv(latitude_SRA_concat2, file="latitude_forConcatenate_n2.csv", col.names=FALSE, row.names=FALSE)

# create updated SRA metadata
latitude_SRA_concat_updated = latitude_SRA_concat %>% 
  mutate(Run = Genotype) %>% 
  separate(Genotype, into = c('Location', 'Genotype'), sep = 1, convert = TRUE) %>% 
  select(Run, Location) %>% 
  mutate(Location = if_else(Location=="K", "Keppel", "Orpheus"))

write.csv(latitude_SRA_concat_updated, file="SraRunTable_latitude_updated.csv")

latitude_SRA_concat_updated = latitude_SRA_concat_updated %>% 
  dplyr::rename(num=Run, treat=Location)

# metadata
latitude_SRA_updated = latitude_SRA %>% 
  select(Run, Genotype, treatment) %>% 
  separate(Genotype, into = c('Location', 'Genotype'), sep = 1, convert = TRUE) %>% 
  separate(treatment, into = c('Origin', 'Transplant'), sep = "placed at ") %>% 
  select(-c(Genotype, Origin)) %>% 
  dplyr::rename(num=Run, treat=Location)

write.csv(latitude_SRA_updated, file="Metadata_latitude.csv")

# create population lists

# orpheus
latitude_orpheus = latitude_SRA_concat %>% 
  mutate(Run = Genotype) %>% 
  separate(Genotype, into = c('Location', 'Genotype'), sep = 1, convert = TRUE) %>% 
  select(Run, Location, Genotype) %>% 
  filter(Location=="O")
  
latitude_orpheus_names = latitude_orpheus %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

# keppel
latitude_keppel = latitude_SRA_concat %>% 
  mutate(Run = Genotype) %>% 
  separate(Genotype, into = c('Location', 'Genotype'), sep = 1, convert = TRUE) %>% 
  select(Run, Location, Genotype) %>% 
  filter(Location=="K")

latitude_keppel_names = latitude_keppel %>% 
  select(Run) %>% 
  mutate(Run = sub("$", ".bam", Run)) # add ".bam" to end of string

write.csv(latitude_orpheus_names, file="SRAnames_latitude_pop_orpheus.csv")
write.csv(latitude_keppel_names, file="SRAnames_latitude_pop_keppel.csv")
