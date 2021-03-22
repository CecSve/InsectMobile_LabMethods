setwd("H:/Documents/Insektmobilen/PhD Courses/Stats for BioScience II") # error, but set wd manually

#### set environment ####
library(tidyverse)
library(tidyselect)
library(stringr)
library(data.table)
library(vegan)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(dplyr)

#### load & subset data ####
# metadata
metadata <- read.table("metadata_novaseqrun_mod_mergedsizes.txt", sep = "\t", header = TRUE, row.names = 1)
metadata <- metadata %>% filter(Primer == 'fwh') %>% filter(Sample_type %in% c('ethanol','carnet'))  # select samples that are amplified with the fwh primer

# the extraction comparison needs to be based on the same samples, so the corresponding ethanol and non-destructive samples needs to be extracted
etoh <- metadata %>% filter(Sample_type == 'ethanol')
etoh <- etoh[-c(18:22), ] # remove the samples that have been extracted on 5 탅
etoh <-
  etoh %>% mutate(
    SampleID = recode(
      SampleID,
      "P74.1A_50킠" = "P74.1A",
      "P67.2A_50킠" = "P67.2A",
      "P66.2A_50킠" = "P66.2A",
      "P87.2B_50킠" = "P87.2B",
      "P87.1B_50킠" = "P87.1B"
    )
  ) # change the names to not contain volume input

etoh <- droplevels(etoh)
samples <- etoh$SampleID # extract the samples names for comparison with non-destructive extracted samples

nd <- metadata[(grep("P176.2A|P176.2B|P184.2B|P21.2A|P25.1A|P3.2A|P34.2A|P34.2B|P35.2A|P36.2B|P65.2B|P66.1B|P74.1A|P67.2A|P66.2A|P87.2B|P87.1B", metadata$SampleID)), ] # select samples that partially match samples

nd <- nd %>% filter(Sample_type == "carnet") # more samples because the have been split into size fractions
metadata <- merge(nd, etoh, all = TRUE) # combine for one dataset
metadata <-
  metadata %>% mutate(
    SampleID_nosize = recode(
      SampleID_nosize,
      "P74.1A_50킠" = "P74.1A",
      "P67.2A_50킠" = "P67.2A",
      "P66.2A_50킠" = "P66.2A",
      "P87.2B_50킠" = "P87.2B",
      "P87.1B_50킠" = "P87.1B"
    )
  ) # change the names to not contain volume input

metadata <- metadata %>% select(PCRID, SampleID_nosize, Sample_type, Extraction_method, Purification_method, Primer)
metadata <- unique(metadata)

# one outlier occurs in the mds, so that sample will be excluded
metadata <-  metadata[!metadata$SampleID_nosize %in% "P184.2B",]
metadata <- droplevels(metadata)
samples <- metadata$PCRID # create a factor for subsetting the otu table

# sequence data
# first we will include the otu table
lulified_fwh_nochim <-
  readRDS("raw-data/lulified_fwh_nochim.RDS") # read in the lulufied RDS file

otus <-
  lulified_fwh_nochim[["curated_table"]] # extract the otutable

# some samples are split into two for the non-destructive method if they have a large size fraction, so those samples needs to be combined as a 'total' sample for the extraction method comparison
# It's not completely transparent which samples should be combined, but I have organised them in a spreadsheet so I can assemble them based on column names
otus$P176.2B <-
  otus$IM18_279 + otus$IM18_278

otus$P21.2A <-
  otus$IM18_282 + otus$IM18_281

otus$P25.1A <-
  otus$IM18_284 + otus$IM18_283

otus$P34.2A <-
  otus$IM18_287 + otus$IM18_286

otus$P34.2B <-
  otus$IM18_289 + otus$IM18_288

otus$P36.2B <-
  otus$IM18_292 + otus$IM18_291

otus$P65.2B <-
  otus$IM18_294 + otus$IM18_293

otus$P66.2A <-
  otus$IM18_341 + otus$IM18_340

otus$P74.1A <-
  otus$IM18_338 + otus$IM18_337

# remove the individual samples before proceeding
otus <-
  select(
    otus,-c(
      IM18_278,
      IM18_279,
      IM18_281,
      IM18_282,
      IM18_283,
      IM18_284,
      IM18_286,
      IM18_287,
      IM18_288,
      IM18_289,
      IM18_291,
      IM18_292,
      IM18_293,
      IM18_294,
      IM18_340,
      IM18_341,
      IM18_337,
      IM18_338
    )
  )


otus <- otus[, (names(otus) %in% samples)]

otus <-
  otus[apply(otus[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

# metadata rownames need to match the column headers in the otu data

occurences <-
  rownames(otus) # save occurence IDs so they can be used to subset taxonomy

# taxonomy data
tax <-
  read.table(
    "raw-data/blastresult_GBIF_lulufied_fwh_novaseq.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )
taxonomy <- tax[(rownames(tax) %in% occurences),]
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package
taxonomy <-
  taxonomy %>% rownames_to_column('otu') %>% filter(identity >= 99) %>% filter(class == "Insecta") %>% column_to_rownames('otu') # choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames

# the filtrated taxonomy occurences needs to be removed from otus
occurences <- rownames(taxonomy)
otus <- otus[(rownames(otus) %in% occurences),]
otus <- droplevels(otus)

tax <- tax[(rownames(tax) %in% occurences),] # subset original loaded taxonomy
tax <- tax %>% select(classification) # format needed for phyloseq

### applying stats from sampling method paper instead - based on presence absence data ###############

# get stats for taxonomic assignment
taxonomy %>% group_by(order) %>% tally()
apply(taxonomy, 2, function(x){ length(which(!is.na(x))) })

otus <- decostand(otus, "pa") # make otutable into presence absence

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(otus, var = "otuid")
taxonomy <- rownames_to_column(taxonomy, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

data <- left_join(otus, taxonomy, by = "otuid")

# working with reshape2 to wrangle data from wide to long format
str(data)
data <- data %>% select(-c(marker, identity, bitScore, expectValue,matchType, scientificName, sequence))
data <- column_to_rownames(data, var = "otuid")
longdata <- melt(data) 
longdata <- longdata[longdata$value!=0,]
test <- longdata %>%
  group_by(variable) %>%
  summarise(richness=sum(value))
#longdata <- aggregate(. ~ family + variable, data = longdata, FUN = sum) # if values in the insect group column and the sample (variable) column are identical, then sum up how many unique otus there were in the sample from that insect group (richness)

# join metadata and longdata
plot_data <- merge(test, metadata, by.x = "variable", by.y = "PCRID")
plot_data <- as_tibble(plot_data)
plot_data <- plot_data %>% dplyr::rename(sampleid = variable, oturichness = richness)

table(plot_data$Sample_type)
str(plot_data)

hist(plot_data$oturichness)
hist(log(plot_data$oturichness))

### plot richness ################

data <- plot_data %>% select(-c(sampleid, Sample_type, Purification_method)) # remove columns before widening data

# generate columns for each sample types oturichness
data <- data %>%
  pivot_wider(names_from = Extraction_method, values_from = oturichness, values_fill = 0)

# Convert order as a grouping variable
#data$order <- as.factor(data$order)

# basic plot
ggscatter(
  data,
  x = "nondestructive",
  y = "ethanol",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "darkgrey", fill = "lightgray")
) +  stat_cor(method = "pearson") 

ggscatter(data, x = "nondestructive", y = "ethanol",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "darkgrey", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coeff.args = list(method = "pearson")
) + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) #+ xscale("log", .format = T) + yscale("log", .format = T)

# scaling does not improve
qqnorm(log(data$nondestructive))
qqnorm(log(data$ethanol))
