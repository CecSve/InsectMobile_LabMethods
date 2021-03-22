setwd("H:/Documents/Insektmobilen/PhD Courses/Stats for BioScience II") # error, but set wd manually

#### set environment ####
library(tidyverse)
library(tidyselect)
library(stringr)
library(data.table)
library(ggrepel)
library(reshape2)
library(dplyr)

### metadata #################
metadata <- read.table("cleaned-data/metadata_novaseq_purification_primer_comparison.txt", sep = "\t", header = TRUE, row.names = 1)
metadata <- metadata %>% filter(purification == 'qiasymphony') # only select samples purified with qiasymphony
table(metadata$primer)

samples <- read.table("cleaned-data/PID_SampleID_PCRID3_qubit.txt", sep = "\t", header = TRUE, row.names = 1)
samples <- rownames_to_column(samples, var = "PCRID")

data <- merge(samples, metadata, by = "PCRID")
table(data$primer) # 61 samples for each primer to compare

# remove blanks, tests of nets and etoh blanks
remove.list <- paste(c("EtOH", "ExBlank", "Net_"), collapse = '|')
data <- data %>% filter(!grepl(remove.list, PID))
data <- data %>% 
  filter(!grepl('50µL', SampleID))
data <- data %>% 
  filter(!grepl('5µL', SampleID))

table(data$primer)
table(data$SampleID)
data <- droplevels(data)

# make sure samples match between primer pairs
art <- data %>% filter(primer == "art") #%>% select(PCRID, SampleID) #%>% rename("PCRID_art"=PCRID)
fwh <- data %>% filter(primer == "fwh") #%>% select(PCRID, SampleID) #%>% rename("PCRID_fwh"=PCRID)
ins <- data %>% filter(primer == "ins") #%>% select(PCRID, SampleID) #%>% rename("PCRID_ins"=PCRID)

setdiff(fwh$PCRID, ins$PCRID)
setdiff(ins$PCRID, fwh$PCRID)
setdiff(art$SampleID, fwh$SampleID)

### ins ASVs ################################################
lulified_nochim <-
  readRDS("raw-data/lulified_ins_nochim.RDS") # read in the lulufied RDS file

asvs_ins <-
  lulified_nochim[["curated_table"]] # extract the otutable

asvs_ins <- asvs_ins[colnames(asvs_ins) %in% ins$PCRID]
rowSums(asvs_ins) # some zeros can be found
colSums(asvs_ins)
asvs_ins <- asvs_ins[, colSums(asvs_ins != 0) > 0]
asvs_ins <-
  asvs_ins[apply(asvs_ins[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (ASVs that are not present in the subsetted samples)

### ins taxonomy data #######################################
tax <-
  read.table(
    "raw-data/ins_taxonomy&samples_final_nofastqc.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )
str(tax)
taxonomy <- tax[(rownames(tax) %in% rownames(asvs_ins)),]
taxonomy_ins <-
  taxonomy %>% rownames_to_column('otu') %>% filter(pident >= 99) %>% filter(class == "Insecta") %>% select(otu, pident, phylum, class, order, family, genus, species, alternatives, redundancy) %>% column_to_rownames('otu') # choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames

occurences <- rownames(taxonomy_ins)
tax <- tax[(rownames(tax) %in% occurences),] # subset original loaded taxonomy
tax_ins <- tax %>% 
  unite(classification, kingdom, phylum, class, order, family, genus, species, sep = "_", remove = FALSE)

tax_ins <- tax_ins %>% select(classification) # format needed for phyloseq

# the filtrated taxonomy occurences needs to be removed from asvs
asvs_ins <- asvs_ins[(rownames(asvs_ins) %in% occurences),]
rowSums(asvs_ins)
colSums(asvs_ins)
asvs_ins <- droplevels(asvs_ins)

### fwh ASVs ################################################
lulified_nochim <-
  readRDS("raw-data/lulified_fwh_nochim.RDS") # read in the lulufied RDS file

asvs_fwh <-
  lulified_nochim[["curated_table"]] # extract the otutable

asvs_fwh <- asvs_fwh[colnames(asvs_fwh) %in% fwh$PCRID]

### art ASVs ################################################
lulified_nochim <-
  readRDS("raw-data/lulified_art_nochim.RDS") # read in the lulufied RDS file

asvs_art <-
  lulified_nochim[["curated_table"]] # extract the otutable

asvs_art <- asvs_art[colnames(asvs_art) %in% art$PCRID]

# ins only has 37 variables after subset, so we need to subset fwh and art based on those samples
asvs_art <- asvs_art[colnames(asvs_art) %in% colnames(asvs_ins)]
asvs_fwh <- asvs_fwh[colnames(asvs_fwh) %in% colnames(asvs_ins)]

rowSums(asvs_fwh) # some zeros can be found
colSums(asvs_fwh)

asvs_fwh <-
  asvs_fwh[apply(asvs_fwh[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (ASVs that are not present in the subsetted samples)

rowSums(asvs_art) # some zeros can be found
colSums(asvs_art)

asvs_art <-
  asvs_art[apply(asvs_art[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (ASVs that are not present in the subsetted samples)

### fwh taxonomy data #######################################
tax <-
  read.table(
    "raw-data/blastresult_GBIF_lulufied_fwh_novaseq.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )
taxonomy <- tax[(rownames(tax) %in% rownames(asvs_fwh)),]
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package
taxonomy_fwh <-
  taxonomy %>% rownames_to_column('otu') %>% filter(identity >= 99) %>% filter(class == "Insecta") %>% select(otu, marker, identity, bitScore, expectValue, matchType, scientificName, phylum, class, order, family, genus, species, sequence) %>%  column_to_rownames('otu') # choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames

# the filtrated taxonomy occurences needs to be removed from asvs
occurences <- rownames(taxonomy_fwh)
asvs_fwh <- asvs_fwh[(rownames(asvs_fwh) %in% occurences),]
asvs_fwh <- droplevels(asvs_fwh)

tax <- tax[(rownames(tax) %in% occurences),] # subset original loaded taxonomy
tax_fwh <- tax %>% select(classification) # format needed for phyloseq

### art taxonomy data #######################################
tax <-
  read.table(
    "raw-data/art_lib1-9_lulufied_blastresult.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )
taxonomy <- tax[(rownames(tax) %in% rownames(asvs_art)),]
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package
taxonomy_art <-
  taxonomy %>% rownames_to_column('otu') %>% filter(identity >= 99) %>% filter(class == "Insecta" )%>% select(otu, marker, identity, bitScore, expectValue, matchType, scientificName, phylum, class, order, family, genus, species, sequence) %>% column_to_rownames('otu') # choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames

# the filtrated taxonomy occurences needs to be removed from asvs
occurences <- rownames(taxonomy_art)
asvs_art <- asvs_art[(rownames(asvs_art) %in% occurences),]
asvs_art <- droplevels(asvs_art)

tax <- tax[(rownames(tax) %in% occurences),] # subset original loaded taxonomy
tax_art <- tax %>% select(classification) # format needed for phyloseq

### make a new ID column in metadata and asv tables #############
# metadata
art <- art %>% group_by(PCRID, SampleID, primer) %>% mutate(PCRID_primer = paste(PCRID, primer, sep = "_"))
ins <- ins %>% group_by(PCRID, SampleID, primer) %>% mutate(PCRID_primer = paste(PCRID, primer, sep = "_"))
fwh <- fwh %>% group_by(PCRID, SampleID, primer) %>% mutate(PCRID_primer = paste(PCRID, primer, sep = "_"))

# reorder columns
colnames(art)
art <- art[,c(9,1,2,3,4,5,6,7,8)]
ins <- ins[,c(9,1,2,3,4,5,6,7,8)]
fwh <- fwh[,c(9,1,2,3,4,5,6,7,8)]

metadata <- rbind(art, ins, fwh)

# rename columns in asv tables
asvs_art <- asvs_art[colnames(asvs_art) %in% art$PCRID]
asvs_ins <- asvs_ins[colnames(asvs_ins) %in% ins$PCRID]
asvs_fwh <- asvs_fwh[colnames(asvs_fwh) %in% fwh$PCRID]

colnames(asvs_fwh) <- fwh$PCRID_primer[match(names(asvs_fwh),fwh$PCRID)]
colnames(asvs_ins) <- ins$PCRID_primer[match(names(asvs_ins),ins$PCRID)]
colnames(asvs_art) <- art$PCRID_primer[match(names(asvs_art),art$PCRID)]

# make rownames to column prior to merge
listdfs <- list(asvs_art, asvs_fwh, asvs_ins)

temp <- lapply(listdfs, function(x){
  x %>% rownames_to_column(var = 'asvid')
})

names(temp) <- c("asvs_art", "asvs_fwh", "asvs_ins")
asvs_art <- temp$asvs_art
asvs_fwh <- temp$asvs_fwh
asvs_ins <- temp$asvs_ins

test <- merge(asvs_art, asvs_fwh, by = "asvid", all = T)
asvtable <- merge(test, asvs_ins, by = "asvid", all = T)
#test2 <- semi_join(test, asvs_ins, by = "asvid") are there any asv matches? --> no

# there's not the same amount of samples in the asv table as in the metadata table
remove.list <- setdiff(metadata$PCRID_primer, colnames(asvtable))
metadata <- metadata[!metadata$PCRID_primer %in% remove.list, ]

# taxonomy
listdfs <- list(taxonomy_art, taxonomy_fwh, taxonomy_ins)

temp <- lapply(listdfs, function(x){
  x %>% rownames_to_column(var = 'asvid')
})

names(temp) <- c("taxonomy_art", "taxonomy_fwh", "taxonomy_ins")
taxonomy_art <- temp$taxonomy_art
taxonomy_fwh <- temp$taxonomy_fwh
taxonomy_ins <- temp$taxonomy_ins

str(taxonomy_art)
test <- merge(taxonomy_art, taxonomy_fwh, all = T)
taxonomytable <- merge(test, taxonomy_ins, all = T)
#test2 <- semi_join(taxonomy_art, taxonomy_fwh, by = "asvid") #are there any asv matches? --> no
taxonomytable$marker[is.na(taxonomytable$marker)] <- "16S" # set marker to 16S for NAs

# tax
listdfs <- list(tax_art, tax_fwh, tax_ins)

temp <- lapply(listdfs, function(x){
  x %>% rownames_to_column(var = 'asvid')
})

names(temp) <- c("tax_art", "tax_fwh", "tax_ins")
tax_art <- temp$tax_art
tax_fwh <- temp$tax_fwh
tax_ins <- temp$tax_ins

str(tax_art)
test <- merge(tax_art, tax_fwh, all = T)
taxtable <- merge(test, tax_ins, all = T)
#test2 <- semi_join(tax_art, tax_fwh, by = "asvid") #are there any asv matches? --> no

# Now the complete metadata, taxonomy and asv table are ready for analysis
write.table(metadata, file = "cleaned-data/metadata_primer_comparison.txt", sep = "\t", col.names = T, row.names = F)

asvtable[is.na(asvtable)] <- 0 # make NAs into zeros
write.table(asvtable, file = "cleaned-data/asvtable_primer_comparison.txt", sep = "\t", col.names = T, row.names = F)

write.table(taxonomytable, file = "cleaned-data/taxonomy_primer_comparison.txt", sep = "\t", col.names = T, row.names = F)

write.table(taxtable, file = "cleaned-data/tax_primer_comparison.txt", sep = "\t", col.names = T, row.names = F)
