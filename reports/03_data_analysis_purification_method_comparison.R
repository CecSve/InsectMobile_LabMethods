---
  title: "Evaluating methodological steps in bulk insect DNA metabarcoding for large-scale biodiversity studies"
subtitle: 'Data analysis for DNA purification method comparison'
author: "Cecilie S. Svenningsen"
date: "28/02/2020"
output:
  pdf_document: default
html_document: default
editor_options: 
  chunk_output_type: console
---


setwd("H:/Documents/Insektmobilen/PhD Courses/Stats for BioScience II")

#### set environment ####
library(tidyverse)
library(tidyselect)
library(stringr)
library(data.table)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ampvis2)
library(hrbrthemes)
library(ggpubr)
library(reshape2)
library(emmeans)
library(indicspecies)
library(lme4)
library(lmerTest)
library(MuMIn)

#### load & subset data ####
# metadata
metadata <- read.table("metadata_novaseqrun_mod_mergedsizes.txt", sep = "\t", header = TRUE, row.names = 1)

metadata <- metadata %>% filter(Primer == 'fwh') %>% filter(Purification_method %in% c('qiasymphony','qiaquick')) %>% filter(Sample_type %in% c('carnet', 'malaise'))  # select samples that are purified with qiasymphony + qiaquick and amplified with the fwh primer and from real samples (not contamination tests, pcr negatives, blanks etc)

### merging total samples from size sorted samples ###########
keep <- metadata %>% select(PCRID, SampleID_nosize, Purification_method)
metadata_unique <- metadata %>% distinct(SampleID, Purification_method, .keep_all = TRUE)

# sequence data
lulified_fwh_nochim <-
  readRDS("raw-data/lulified_fwh_nochim.RDS") # read in the lulufied RDS file

otus <-
  lulified_fwh_nochim[["curated_table"]] # extract the otutable

t.asvs <- t(otus)
t.asvs <- as.data.frame(t.asvs) %>% rownames_to_column(var = "PCRID") 
test <- left_join(keep, t.asvs, by = "PCRID")
test <- test %>% drop_na(SampleID_nosize)

test2 <- test %>% dplyr::select(-PCRID) %>% group_by(SampleID_nosize, Purification_method) %>% summarise_all(list(sum))

test3 <- merge(metadata_unique, test2, by=c("SampleID_nosize", "Purification_method")) 
asvs <- test3 %>% select(PCRID, 8:4987)
metadata <- test3 %>% select(PCRID, SampleID_nosize, Purification_method)

# now equal sample sizes should be obtained for further analysis
table(metadata$Purification_method) 
qq <- metadata %>% filter(Purification_method == "qiaquick")
qs <- metadata %>% filter(Purification_method == "qiasymphony")

setdiff(qs$SampleID_nosize, qq$SampleID_nosize) 
setdiff(qq$SampleID_nosize, qs$SampleID_nosize)

# match samples
samples <- qq$SampleID_nosize
qs <- qs[(qs$SampleID_nosize) %in% samples, ]
samples <- qs$SampleID_nosize
qq <- qq[(qq$SampleID_nosize) %in% samples, ]

data <- rbind(qs,qq)
table(data$Purification_method)

# two outliers occur in the mds, so the samples will be excluded ----
data <-  data[(!data$PCRID %in% c("IM17_44_QQ", "IM17_44_QS", "IM17_32_QQ", "IM17_32_QS")),]

data <- droplevels(data)
samples <- data$PCRID # create a factor for subsetting the otu table

# sequence data
otus <- asvs # based on merged size sorted samples for total sample read counts
otus[is.na(otus)] <- 0 # some NAs were introduced - recode them to zeros
otus <- otus[(otus$PCRID %in% samples), ]
str(otus)
rownames(otus) <- NULL
otus <- otus %>% column_to_rownames(var = "PCRID") 
otus <- t(otus)
otus <- as.data.frame(otus)
rowSums(otus)

otus <-
  otus[apply(otus[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

occurences <-
  rownames(otus) # save occurence IDs so they can be used to subset taxonomy

# taxonomy data
tax <-
  read.table(
    "raw-data/blastresult_length_corrected_fwh.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )
taxonomy <- tax[(rownames(tax) %in% occurences),]
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package
taxonomy <-
  taxonomy %>% rownames_to_column('otu') %>% filter(identity >= 99) %>% filter(phylum == "Arthropoda") %>% column_to_rownames('otu') # choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames

# the filtrated taxonomy occurences needs to be removed from otus
occurences <- rownames(taxonomy)
otus <- otus[(rownames(otus) %in% occurences),]
colSums(otus)
otus <- otus[,colSums(otus) > 0] # some samples are empty and will be removed
otus <- droplevels(otus)

tax <- tax[(rownames(tax) %in% occurences),] # subset original loaded taxonomy
#tax <- tax %>% select(classification) # format needed for phyloseq

samples <- colnames(otus)
metadata <- data[(data$PCRID) %in% samples, ] # some samples were dropped along the filtration steps and needs to be removed from metadata
metadata <- droplevels(metadata)
table(metadata$Purification_method)

# make sure the same samples are used for comparison between both methods
samplesqs <- metadata %>% filter(Purification_method == "qiasymphony") %>% select(PCRID)
samplesqs$PCRID = gsub(pattern = "_Q.*",
                      replacement = "",
                      x = samplesqs$PCRID) # we need to rename the sample IDs so they don't contain _QS and can match with our other

samplesqq <- metadata %>% filter(Purification_method == "qiaquick") %>% select(PCRID)
samplesqq$PCRID = gsub(pattern = "_Q.*",
                     replacement = "",
                     x = samplesqq$PCRID) # we need to rename the sample IDs so they don't contain _QS and can match with our other

setdiff(samplesqq$PCRID, samplesqs$PCRID)
setdiff(samplesqs$PCRID, samplesqq$PCRID)
metadata <-  metadata[(!metadata$PCRID %in% c("IM17_12_QQ")),]
metadata <-  metadata[(!metadata$PCRID %in% c("IM17_5_QS", "IM17_13_QS", "IM17_19_QS", "IM17_17_QS")),]
table(metadata$Purification_method) 

# even sample sizes
qq <- metadata %>% filter(Purification_method == "qiaquick")
qs <- metadata %>% filter(Purification_method == "qiasymphony")

setdiff(qs$SampleID_nosize, qq$SampleID_nosize) 
setdiff(qq$SampleID_nosize, qs$SampleID_nosize)

# match samples
samples <- qq$SampleID_nosize
qs <- qs[(qs$SampleID_nosize) %in% samples, ]
samples <- qs$SampleID_nosize
qq <- qq[(qq$SampleID_nosize) %in% samples, ]

metadata <- rbind(qs,qq)
samples <- metadata$PCRID

# subset otus
otus <- otus[, (names(otus) %in% samples)]
rowSums(otus)

otus <-
  otus[apply(otus[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

occurences <-
  rownames(otus) # save occurence IDs so they can be used to subset taxonomy

taxonomy <- taxonomy[(rownames(taxonomy) %in% occurences), ]
tax <- tax[(rownames(tax) %in% occurences),]# subset original loaded taxonomy

tax <- tax %>% dplyr::select(classification) # format needed for phyloseq

metadata <- metadata %>% separate(col = PCRID, into = c("year", "sample", "purification_abr"), sep = "\\_")

table(metadata$purification_abr)

#### data analysis ####
# community composition ----
tasvs <- t(otus) # samples as rows, taxa as columns
#bray.data<-vegdist(tasvs, Type = "bray") # 
# the Bray-Curtis index is based on abundance data, while the Sorensen index is based on presence/absence data. Both indices have similarity and dissimilarity (or distance) versions. 
tasvs <- decostand(tasvs, method = "pa") # transform to presence absence
sorensen.data <- designdist(tasvs, "(A+B-2*J)/(A+B)") # for presence absence
set.seed(42)
meta_mds <- metaMDS(sorensen.data, k = 2, trymax = 100)
plot(meta_mds)
meta_mds$points
stressplot(meta_mds)

class(meta_mds)
methods(class="metaMDS")

Y <- scores(meta_mds, display="sites")
all(rownames(tasvs) == rownames(Y)) # make sure that rows of Y are in the same order as rows of norm.data

Y <- data.frame(Y)
Y <- Y %>% rownames_to_column(var = "PCRID")

metadata <- metadata %>% unite("PCRID", year:purification_abr, sep = "_", remove = FALSE)
data <- merge(metadata, Y, by = "PCRID")
str(data)
data$sample <- as.numeric(data$sample)

tpa <- as.data.frame(tasvs)
tpa$richness <- rowSums(tpa) # get richness for each sample

richnessdata <- tpa %>% rownames_to_column(var = "PCRID") %>% select(PCRID, richness)
data <- merge(data, richnessdata, by = "PCRID") # now theres a response value for richness and community composition which can be used for further analysis

shannon <- diversity(tpa[-1], "shannon")
pilous_evenness <- shannon/log(specnumber(tpa))

shannon <- as.data.frame(shannon) %>% rownames_to_column(var = "PCRID")
pilous <- as.data.frame(pilous_evenness) %>% rownames_to_column(var = "PCRID")
data <- merge(data, shannon, by = "PCRID")
data <- merge(data, pilous, by = "PCRID")


custom_col = c("#2e420d","#e031e7", "#47a900", "#4d47f1", "#008c34", "#c177ff","#256200", "#172d9c",  "#ff7f14","#0171cb", "#cf1a00", "#006fb6", "#ff5148","#84a6ff", "#c6a93c",  "#860059", "#7aafe5",  "#8d4200", "#ef89ce", "#512f5c")

# plot all samples with normal confidence ellipse by extraction method (0.95)
nmds <- ggplot(data, aes(x=NMDS1, y=NMDS2, size=richness, colour = SampleID_nosize, group = SampleID_nosize)) + geom_point(show.legend = F) + geom_text_repel(aes(label=purification_abr), show.legend = FALSE) + 
  theme_minimal() + scale_size_continuous(range = c(4,8)) + labs(colour = "Sample ID", size = "Sample \nrichnes") + geom_line(aes(colour = SampleID_nosize), size = 1, show.legend = F) + scale_colour_manual(values = custom_col)

save_plot("plots/nmds_purification.png", nmds, base_height = 8, base_width = 12)

### indicator analysis #########
abund <- tasvs # based on presence absence
method <- data$purification_abr

# analysis with indicator value - explained in De Cáceres et al. (2010) - the accounts for unequal sample sizes
indval <-  multipatt(abund, method, func = "IndVal.g", duleg = TRUE, control = how(nperm=999)) # number of permutations affect the precision of the p-value
summary(indval, indvalcomp=TRUE)

# Adonis test
set.seed(42)
adonis(sorensen.data ~ Purification_method, data = data, strata = metadata$SampleID, permutations = 999)

# ampvis data load ----
# join otus and taxonomy
otuspa <- decostand(otus, "pa") 
otus <- rownames_to_column(otus, var = "otuid")
otuspa <- rownames_to_column(otuspa, var = "otuid")
taxonomy <- taxonomy %>% select(kingdom, phylum, order, family, genus, species)
taxonomy <- rownames_to_column(taxonomy, var = "otuid")
taxonomy$class <- NA
taxonomy <- taxonomy %>% select(otuid, kingdom, phylum, class, order, family, genus, species)
otutable <- left_join(otus, taxonomy, by = "otuid")
otutable <- left_join(otuspa, taxonomy, by = "otuid") # choose between pa or rel abundance
# revert otuids back to rownames if needed
otus <- column_to_rownames(otus, var = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames

# metadata must have sample data in the first column and column classes matter 
str(data)
#metadata <- rownames_to_column(metadata, var = "otuid")

# load amp data
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

amp <- amp_load(otutable = otutable, metadata = metadata)

# heatmap vizulaization ----
# convert read abundances into percentages
#psn <- transform_sample_counts(ps, function(x) x/sum(x) * 100)

# heatmap 
pur_heatmap<- amp_heatmap(
  data = amp,
  group_by = "sample",
  facet_by = "purification_abr",
  tax_aggregate = "Family",
  tax_add = "Order",
  plot_values = FALSE,
  color_vector = c("mintcream","darkblue"),
  plot_colorscale = "sqrt",
  tax_empty = "remove",
  tax_show = 30,
)+ guides(fill = guide_legend(title = "Most frequently \nobserved taxa \n(% per sample)", reverse = T)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_plot("plots/heatmap_purification_relfreq.png", pur_heatmap, base_height = 6, base_width = 10)

# venn diagram: https://madsalbertsen.github.io/ampvis2/reference/amp_venn.html 
amp_venn(amp,
         group_by = "Purification_method",
         cut_a = 0.1,
         cut_f = 1) # core set to 1, standard is 80 but in this case it is more interesting what the overlap is

# boxplot/point plot
amp_boxplot(
  amp,
  group_by = "purification_abr",
  sort_by = "mean",
  plot_type = "boxplot",
  point_size = 3,
  tax_aggregate = "Family",
  tax_add = "Order",
  tax_show = 20,
  tax_empty = "remove",
  plot_log = TRUE,
  plot_flip = TRUE
) + ylab('Abundance  \n(log transformed)') + guides(colour = guide_legend(title = "Purification method", reverse = T)) + scale_color_manual(values = c("#2EDE28", "#DE1282"))

# ordination: https://madsalbertsen.github.io/ampvis2/reference/amp_ordinate.html 
ord <- amp_ordinate(
  amp,
  type = "NMDS",
  distmeasure = "bray",
  transform = "none",
  sample_color_by = "purification_abr",
  sample_colorframe = TRUE
)

ord + scale_color_manual(values = c("#E69F00", "#56B4E9")) + scale_fill_manual(values = c("#E69F00", "#56B4E9"))

# richness analysis ----
# richness indices - calculated for pa data
amp_alphadiv <-
  amp_alphadiv(
    amp,
    measure = c("observed", "shannon", "simpson", "invsimpson"),
    richness = TRUE
  )

otus <- rownames_to_column(otus, var = "otuid")
otutable <- left_join(otus, taxonomy, by = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)
amp_abun <- amp_load(otutable = otutable, metadata = data)

# rarefaction curve
rarecurve <-
  amp_rarecurve(
    amp_abun,
    facet_by = "Purification_method",
    stepsize = 50,
    color_by = "SampleID_nosize",
    facet_scales = "free"
  )

custom_final28 = c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E", "#DDAA77")

rarfac <- rarecurve + scale_color_manual(values = custom_col) + geom_line(size = 2, na.rm = TRUE) + scale_x_continuous(
  labels = scales::number,
  limits = c(0, 50000),
  breaks = seq(from = 0, to = 50000, by = 1000)
) + guides(colour=guide_legend(ncol=8)) + ylab("Number of observed ASVs") + ylim(0, 260) + theme_minimal()+ theme(legend.position = "none", axis.text.x = element_text(angle = 90), legend.title = element_blank()) 

save_plot("plots/rarecurve_purificationmethod.png", rarfac, base_height = 8, base_width = 12)

# plot richness indices on log10 transformed y-axis
amp_alphadiv <- amp_alphadiv %>% dplyr::rename(ObservedASVs = ObservedOTUs)

richplot <- amp_alphadiv %>% tidyr::gather("id", "Richness", c(8:9, 12)) %>% ggplot(., aes(Purification_method, Richness, fill = Purification_method)) + geom_boxplot(show.legend = F) + facet_wrap(~ id, ncol = 5) + scale_y_log10() + scale_fill_manual(values = c("#BFBFBB", "#80807D")) + theme_pubclean() + labs(fill = "Purification method") + xlab("Purification method") + ylab("Richness (log-transformed)")

save_plot("plots/richind_purificationmethod.png", richplot, base_height = 6, base_width = 10)

### eveness #######
hist(data$shannon)
hist(sqrt(data$pilous_evenness))
qqnorm(data$shannon)
qqnorm(data$pilous_evenness)
qqnorm(log(data$pilous_evenness))

# H0: richness in ND = richness in ethanol
# two-sided test
# assume non-equal variances
t.test(pilous_evenness~Purification_method, mu = 0, alt = "two.sided", conf = 0.95, var.eq = F, paired = T, data = data)

data %>% ggplot + geom_boxplot(aes(Purification_method, pilous_evenness)) # so more variation in abundances between different taxa within ND samples = more complex patterns detected

group_by(data, Purification_method) %>%
  summarise(
    count = n(),
    mean = mean(pilous_evenness, na.rm = TRUE),
    median = median(pilous_evenness, na.rm = TRUE),
    sd = sd(pilous_evenness, na.rm = TRUE)
  )

# summary statistics
group_by(amp_alphadiv, Purification_method) %>%
  summarise(
    count = n(),
    mean = mean(ObservedASVs, na.rm = TRUE),
    median = median(ObservedASVs, na.rm = TRUE),
    sd = sd(ObservedASVs, na.rm = TRUE)
  )
# Logistic regression where random effect is accounted for
# get stats for taxonomic assignment
taxonomy %>% group_by(family) %>% tally()
apply(taxonomy, 2, function(x){ length(which(!is.na(x))) })

### richness model #######
qqnorm(data$richness)
qqnorm(log(data$richness))
hist(data$richness)
hist(log(data$richness))

purificationnmodel <- lmer(log(richness)~ Purification_method + (1|SampleID_nosize), data = data)
summary(purificationnmodel)
r.squaredGLMM(purificationnmodel)

ems <- emmeans(purificationnmodel, specs = pairwise ~ Purification_method, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)
