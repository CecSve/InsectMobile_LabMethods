# organizing data for insektmobilen lab method manuscript. 
# sequencing platform: illumnia novaseq
# primers pairs: fwh (CO1), Art (CO1), INS01 (16S)
library(tidyverse)
# begin with the fwh dataset and then for art and ins as well
setwd("H:/Documents/Insektmobilen/PhD Courses/Stats for BioScience II")

# first we will include the otu table
lulified_fwh_nochim <- readRDS("raw-data/lulified_fwh_nochim.RDS") # read in the lulufied RDS file

fwh_otutable <- lulified_fwh_nochim[["curated_table"]]

names(fwh_otutable) = gsub(pattern = "X*", replacement = "", x = names(fwh_otutable)) # remove the Xs that have emerged in the column headers that are numerical

# subset data to create datasets for the different purification methods
fwh_qq_otus <- fwh_otutable[, grepl("_QQ$", names(fwh_otutable))] # make a dataframe for qiaquick (manual) purified samples
fwh_qs_otus <- fwh_otutable[, grepl("_QS$", names(fwh_otutable))] # make a dataframe for qiasymphony (automatic) purified samples

#names(samples_keep) = gsub(pattern = "_Q.*",replacement = "", x = names(samples_keep)) # we need to rename the sample IDs so they don't contain _QS and can match with our other tables

# Ins dataset
lulified_ins_nochim <- readRDS("raw-data/lulified_ins_nochim.RDS") 
ins_otutable <- lulified_ins_nochim[["curated_table"]]
names(ins_otutable) = gsub(pattern = "X*", replacement = "", x = names(ins_otutable)) # remove the Xs that have emerged in the column headers that are numerical

# somehow the fastq.gz file are included and should be removed
ins_otutable <- ins_otutable %>% select(-contains("filt.fastq.gz"))

# subset data to create datasets for the different purification methods
ins_qq_otus <- ins_otutable[, grepl("_QQ$", names(ins_otutable))] # make a dataframe for qiaquick (manual) purified samples
ins_qs_otus <- ins_otutable[, grepl("_QS$", names(ins_otutable))] # make a dataframe for qiasymphony (automatic) purified samples

# removing rows with zeros
ins_qq_otus_omit <-
  ins_qq_otus[apply(ins_qq_otus[, -1], 1, function(x)
    ! all(x == 0)), ] # remove rows that contain only zeros (OTUs that are not present in our subset samples)

ins_qs_otus_omit <-
  ins_qs_otus[apply(ins_qs_otus[, -1], 1, function(x)
    ! all(x == 0)), ]

fwh_qq_otus_omit <-
  fwh_qq_otus[apply(fwh_qq_otus[, -1], 1, function(x)
    ! all(x == 0)), ] # remove rows that contain only zeros (OTUs that are not present in our subset samples)

fwh_qs_otus_omit <-
  fwh_qs_otus[apply(fwh_qs_otus[, -1], 1, function(x)
    ! all(x == 0)), ]

#histograms
hist(t(log(ins_qq_otus_omit)))
hist(t(log(ins_qs_otus_omit)))
hist(t(log(fwh_qq_otus_omit)))
hist(t(log(fwh_qs_otus_omit)))

#qqnorm
qqnorm(t(log(ins_qq_otus_omit)), ylim=c(0,11), xlim=c(1.9,4.5))
qqnorm(t(log(ins_qs_otus_omit)), ylim=c(0,11), xlim=c(1.9,4.5))
qqnorm(t(log(fwh_qq_otus_omit)), ylim=c(0,12), xlim=c(1.9,4))
qqnorm(t(log(fwh_qs_otus_omit)), ylim=c(0,12), xlim=c(1.9,4))

# art dataset
lulified_tab <- readRDS("raw-data/lulified_art_nochim.RDS")

otutable <- lulified_tab[["curated_table"]]
names(otutable) = gsub(pattern = "X*", replacement = "", x = names(otutable)) # remove the Xs that have emerged in the column headers that are numerical

# subset data to create datasets for the different purification methods
art_qq_otus <- otutable[, grepl("_QQ$", names(otutable))] # make a dataframe for qiaquick (manual) purified samples
art_qs_otus <- otutable[, grepl("_QS$", names(otutable))] # make a dataframe for qiasymphony (automatic) purified samples

# removing rows/observations with zeros
art_qq_otus_omit <-
  art_qq_otus[apply(art_qq_otus[, -1], 1, function(x)
    ! all(x == 0)), ] # remove rows that contain only zeros (OTUs that are not present in our subset samples)

art_qs_otus_omit <-
  art_qs_otus[apply(art_qs_otus[, -1], 1, function(x)
    ! all(x == 0)), ]

# vizualise art data
hist(t(log(art_qq_otus_omit)))
hist(t(log(art_qs_otus_omit)))

qqnorm(t(log(art_qq_otus_omit)), ylim=c(0,13), xlim=c(1.9,4))
qqnorm(t(log(art_qs_otus_omit)), ylim=c(0,12), xlim=c(1.9,4))

# write output
write.table(art_qq_otus_omit, file = "cleaned-data/art_qq_otutable.txt", col.names = NA, sep = "\t")
write.table(art_qs_otus_omit, file = "cleaned-data/art_qs_otutable.txt", col.names = NA, sep = "\t")
write.table(fwh_qq_otus_omit, file = "cleaned-data/fwh_qq_otutable.txt", col.names = NA, sep = "\t")
write.table(fwh_qs_otus_omit, file = "cleaned-data/fwh_qs_otutable.txt", col.names = NA, sep = "\t")
write.table(ins_qq_otus_omit, file = "cleaned-data/ins_qq_otutable.txt", col.names = NA, sep = "\t")
write.table(ins_qs_otus_omit, file = "cleaned-data/ins_qs_otutable.txt", col.names = NA, sep = "\t")
