---
  title: "Evaluating methodological steps in bulk insect DNA metabarcoding for large-scale biodiversity studies"
subtitle: 'Data analysis for DNA extraction method comparison'
author: "Cecilie S. Svenningsen"
date: "19/02/2020"
output:
  pdf_document: default
html_document: default
editor_options: 
  chunk_output_type: console
---

#### set environment ####
library(tidyverse)
library(tidyselect)
library(stringr)
library(data.table)
library(vegan)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(ampvis2)
library(hrbrthemes)
library(ggpubr)
library(reshape2)
library(lme4)
library(lmerTest)
library(emmeans)
library(cowplot)
library(MuMIn)
library(indicspecies)

#### load & subset data ####
# metadata
metadata <- read.table("metadata_novaseqrun_mod_mergedsizes.txt", sep = "\t", header = TRUE, row.names = 1)
metadata <- metadata %>% filter(Primer == 'fwh') %>% filter(Sample_type %in% c('ethanol','carnet'))  # select samples that are purified with qiasymphony and amplified with the fwh primer

# the extraction comparison needs to be based on the same samples, so the corresponding ethanol and non-destructive samples needs to be extracted
etoh <- metadata %>% filter(Sample_type == 'ethanol')
etoh <- etoh[-c(18:22), ] # remove the samples that have been extracted on 5 µl
etoh <-
  etoh %>% mutate(
    SampleID = recode(
      SampleID,
      "P74.1A_50µL" = "P74.1A",
      "P67.2A_50µL" = "P67.2A",
      "P66.2A_50µL" = "P66.2A",
      "P87.2B_50µL" = "P87.2B",
      "P87.1B_50µL" = "P87.1B"
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
      "P74.1A_50µL" = "P74.1A",
      "P67.2A_50µL" = "P67.2A",
      "P66.2A_50µL" = "P66.2A",
      "P87.2B_50µL" = "P87.2B",
      "P87.1B_50µL" = "P87.1B"
    )
  ) # change the names to not contain volume input

metadata <- metadata %>% select(PCRID, SampleID_nosize, Sample_type, Extraction_method, Purification_method, Primer)
metadata <- unique(metadata)

# one outlier occurs in the mds, so that sample will be excluded
metadata <-  metadata[!metadata$SampleID_nosize %in% "P184.2B",]
#metadata <- metadata %>% mutate(SampleID = recode(SampleID, "P176.2AS" =	"P176.2A", "P176.2BS" =	"P176.2B", "P176.2BL" =	"P176.2B", "P184.2BS" =	"P184.2B", "P21.2AS" =	"P21.2A", "P21.2AL" =	"P21.2A", "P25.1AS" =	"P25.1A"))
metadata <- droplevels(metadata)
samples <- metadata$PCRID # create a factor for subsetting the otu table

# sequence data
# first we will include the otu table
lulified_fwh_nochim <-
  readRDS("raw-data/lulified_fwh_nochim.RDS") # read in the lulufied RDS file

otus <-
  lulified_fwh_nochim[["curated_table"]] # extract the otutable

# samples are split into two for the non-destructive method if they have a large size fraction, so those samples needs to be combined as a 'total' sample for the extraction method comparison
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
otus <- droplevels(otus)

tax <- tax[(rownames(tax) %in% occurences),] # subset original loaded taxonomy
tax <- tax %>% select(classification) # format needed for phyloseq

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

meta_mds$stress
Y <- scores(meta_mds, display="sites")
all(rownames(tasvs) == rownames(Y)) # make sure that rows of Y are in the same order as rows of norm.data

Y <- data.frame(Y)
Y <- Y %>% rownames_to_column(var = "PCRID") 

data <- merge(metadata, Y, by = "PCRID")

tpa <- as.data.frame(tasvs)
tpa$richness <- rowSums(tpa) # get richness for eack sample

richnessdata <- tpa %>% rownames_to_column(var = "PCRID") %>% select(PCRID, richness)
data <- merge(data, richnessdata, by = "PCRID") # now theres a response value for richness and community composition which can be used for further analysis

# Adonis test
set.seed(42)
perm <- adonis(
  sorensen.data ~ Extraction_method,
  data = metadata,
  strata = metadata$SampleID_nosize
) # This output tells us that our adonis test is significant so we can reject the null hypothesis that our two extraction methods have the same centroid. Marginal effect of extraction - most effect of sample variation (not surprising)

perm
hist(perm$f.perms)

extraction_cols <- c("#e68dd1","#00723d","#98004c", "#e39d39", "#8c88bd", "#ff4843","#005a88", "#881600", "#47aaff", "#ba002f","#0051cf", "#b8ae38", "#93008d", "#6f9900", "#ff2fdd","#f98d7d", "#c746f8")

nmds <- ggplot(data, aes(x=NMDS1, y=NMDS2, size=richness, colour = SampleID_nosize, group = SampleID_nosize)) + geom_point(show.legend = F) + geom_text_repel(aes(label=replace(Extraction_method, Extraction_method == "nondestructive", "ND")), show.legend = FALSE) + 
  theme_minimal() + scale_size_continuous(range = c(4,8)) + labs(colour = "Sample ID", size = "Sample \nrichnes") + geom_line(aes(colour = SampleID_nosize), size = 1, show.legend = F) + scale_colour_manual(values = extraction_cols)

save_plot("plots/nmds_extraction.png", nmds, base_height = 8, base_width = 12)

### indicator analysis #########
abund <- tasvs # based on presence absence
method <- data$Extraction_method

# analysis with indicator value - explained in De C?ceres et al. (2010) - the accounts for unequal sample sizes
indval <-  multipatt(abund, method, func = "IndVal.g", duleg = TRUE, control = how(nperm=999)) # number of permutations affect the precision of the p-value
summary(indval, indvalcomp=TRUE)

# ampvis data load ----
# join otus and taxonomy
otuspa <- decostand(otus, "pa") # make otutable into presence absence, not needed for rarecurve
otuspa <- rownames_to_column(otuspa, var = "otuid") 
taxonomy <- taxonomy %>% select(kingdom, phylum, order, family, genus, species)
taxonomy <- rownames_to_column(taxonomy, var = "otuid")
taxonomy$class <- NA
taxonomy <- taxonomy %>% select(otuid, kingdom, phylum, class, order, family, genus, species)
otus <- rownames_to_column(otus, var = "otuid") 
otuspa <- rownames_to_column(otuspa, var = "otuid")
otutable <- left_join(otuspa, taxonomy, by = "otuid") # choose otuspa if you want pa instead
# revert otuids back to rownames if needed
otuspa <- column_to_rownames(otuspa, var = "otuid")
otus <- column_to_rownames(otus, var = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames

# metadata must have sample data in the first column and column classes matter 
str(metadata)
#metadata <- rownames_to_column(metadata, var = "otuid")

metadata$exmet <- NULL
metadata$exmet[metadata$Extraction_method=="nondestructive"]<-'ND'
metadata$exmet[metadata$Extraction_method=="ethanol"]<-'ethanol'

# load amp data
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

amp <- amp_load(otutable = otutable, metadata = metadata)

# heatmap vizulaization ----
# convert read abundances into percentages
#psn <- transform_sample_counts(ps, function(x) x/sum(x) * 100)

library(RColorBrewer)
(brewer.pal(7,"RdPu"))

# heatmap (normalise gives the proportion of the sample)
ex_heatmap<- amp_heatmap(
  data = amp,
  group_by = "SampleID_nosize",
  facet_by = "exmet",
  normalise = T,
  tax_aggregate = "Family",
  tax_add = "Order",
  plot_values = F,
  plot_colorscale = "sqrt",
  color_vector = c("grey90", "#FA9FB5", "#DD3497", "#7A0177"),
  tax_empty = "remove",
  tax_show = 30,
) + guides(fill = guide_legend(title = "Most frequently \nobserved taxa \n(% per sample)", reverse = T)) +theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

save_plot("plots/heatmap_extractionmethod_relfreq.png", ex_heatmap, base_height = 6, base_width = 10)

mean(rowSums(tpa))
median(rowSums(tpa))

# venn diagram: https://madsalbertsen.github.io/ampvis2/reference/amp_venn.html 
amp_venn(amp,
         group_by = "exmet",
         cut_a = 0.1, # Abundance cutoff in percent. OTU's below this abundance are excluded from the analysis. (default: 0.1)
         cut_f = 1) # Frequency cutoff in percent. OTU's within the top cut_f of the reads are considered a "core" OTU. (default: 80) 

# boxplot/point plot
amp_boxplot(
  amp,
  group_by = "exmet",
  sort_by = "mean",
  plot_type = "boxplot",
  point_size = 3,
  tax_aggregate = "Family",
  tax_add = "Order",
  tax_show = 20,
  tax_empty = "remove",
  plot_log = TRUE,
  plot_flip = TRUE
) + ylab('Abundance  \n(log transformed)') + guides(colour = guide_legend(title = "Exraction method", reverse = T)) + scale_color_manual(values = c("#2EDE28", "#DE1282"))

# ordination: https://madsalbertsen.github.io/ampvis2/reference/amp_ordinate.html 
ord <- amp_ordinate(
  amp,
  type = "NMDS",
  distmeasure = "bray",
  transform = "none",
  sample_color_by = "Sample_type",
  sample_colorframe = TRUE
)

ord + scale_color_manual(values = c("#E69F00", "#56B4E9")) + scale_fill_manual(values = c("#E69F00", "#56B4E9"))

# richness analysis ----
# richness indices
amp_alphadiv <-
  amp_alphadiv(
    amp,
    measure = c("observed", "shannon", "simpson", "invsimpson"),
    richness = TRUE
  )

### rarefaction curve on read abundance ############
# join otus and taxonomy
otus <- rownames_to_column(otus, var = "otuid")
taxonomy <- rownames_to_column(taxonomy, var = "otuid")
taxonomy <- taxonomy %>% select(otuid, kingdom, phylum, class, order, family, genus, species)
otutable <- left_join(otus, taxonomy, by = "otuid")
# revert otuids back to rownames if needed
otus <- column_to_rownames(otus, var = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames

# metadata must have sample data in the first column and column classes matter 
str(metadata)
#metadata <- rownames_to_column(metadata, var = "otuid")

# load amp data
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

amp <- amp_load(otutable = otutable, metadata = metadata)

# rarefaction curve
rarecurve <-
  amp_rarecurve(
    amp,
    facet_by = "exmet",
    stepsize = 100,
    color_by = "SampleID_nosize",
    facet_scales = "free"
  )

rarfac <- rarecurve + scale_color_manual(values = extraction_cols) + geom_line(size = 2, na.rm = TRUE) + scale_x_continuous(
  labels = scales::number,
  limits = c(0, 50000),
  breaks = seq(from = 0, to = 50000, by = 2000)
) + guides(colour=guide_legend(ncol=8)) + ylab("Number of observed ASVs") + ylim(0, 180) + theme_minimal()+ theme(legend.position = "none", axis.text.x = element_text(angle = 90), legend.title = element_blank())

save_plot("plots/rarecurve_extractionmethod.png", rarfac, base_height = 6, base_width = 10)

# plot richness indices on log10 transformed y-axis
amp_alphadiv <- amp_alphadiv %>% rename(ObservedASVs = ObservedOTUs)

amp_alphadiv$exmet <- NULL
amp_alphadiv$exmet[amp_alphadiv$Extraction_method=="nondestructive"]<-'ND'
amp_alphadiv$exmet[amp_alphadiv$Extraction_method=="ethanol"]<-'ethanol'

richplot <- amp_alphadiv %>% tidyr::gather("id", "Richness", c(8:9, 12)) %>% ggplot(., aes(exmet, Richness, fill = exmet)) + geom_boxplot(show.legend = F) + facet_wrap(~ id, ncol = 5) + scale_y_log10() + scale_fill_manual(values = c("#BFBFBB", "#80807D")) + theme_pubclean() + labs(fill = "Extraction method") + xlab("Extraction method") + ylab("Richness (log-transformed)")

save_plot("plots/richind_extractionmethod.png", richplot, base_height = 6, base_width = 10)

### eveness #######
#Evenness is a measure of how homogeneous or even a community or ecosystem is in terms of the abundances of its species. A community in which all species are equally common is considered even and has a high degree of evenness.

# Pilou evenness (J)	compares the actual diversity value (such as the Shannon-Wiener Index, H') to the maximum possible diversity value (when all species are equally common, Hmax=ln s where S is the total number of species).

# Pilou evenness (J) is constrained between 0 and 1.0 and the more variation in abundances between different taxa within the community, the lower J. Unfortunately, Pilou's J is highly dependent on sample size (since S - the estimated number of species is dependent on sampling effort) and is also highly sensitive to rare taxa. 

asvs <- tpa %>% rownames_to_column(var = "PCRID")
shannon <- diversity(tpa[-1], "shannon")
pilous_evenness <- shannon/log(specnumber(tpa))

shannon <- as.data.frame(shannon) %>% rownames_to_column(var = "PCRID")
pilous <- as.data.frame(pilous_evenness) %>% rownames_to_column(var = "PCRID")
divdata <- merge(data, shannon, by = "PCRID")
divdata <- merge(divdata, pilous, by = "PCRID")

hist(divdata$shannon)
hist(divdata$pilous_evenness)
qqnorm(divdata$shannon)
qqnorm(divdata$pilous_evenness)

# H0: richness in ND = richness in ethanol
# two-sided test
# assume non-equal variances
t.test(pilous_evenness~Extraction_method, mu = 0, alt = "two.sided", conf = 0.95, var.eq = F, paired = T, data = divdata)

divdata %>% ggplot + geom_boxplot(aes(Extraction_method, pilous_evenness)) # so more variation in abundances between different taxa within ND samples = more complex patterns detected

group_by(divdata, Extraction_method) %>%
  summarise(
    count = n(),
    mean = mean(pilous_evenness, na.rm = TRUE),
    median = median(pilous_evenness, na.rm = TRUE),
    sd = sd(pilous_evenness, na.rm = TRUE)
  )

### summary statistics ####
group_by(amp_alphadiv, Extraction_method) %>%
  summarise(
    count = n(),
    mean = mean(ObservedOTUs, na.rm = TRUE),
    median = median(ObservedOTUs, na.rm = TRUE),
    sd = sd(ObservedOTUs, na.rm = TRUE)
  )

# Logistic regression where random effect is accounted for
# get stats for taxonomic assignment
taxonomy %>% group_by(family) %>% tally()
apply(taxonomy, 2, function(x){ length(which(!is.na(x))) })

### richness ##############
divind <- divdata %>% select(PCRID, shannon, pilous_evenness)
ampdiv <- amp_alphadiv %>% select(PCRID, Chao1, exmet)
data <- merge(data, divind, by = "PCRID")
data <- merge(data, ampdiv, by = "PCRID")

hist(sqrt(data$richness))
qqnorm(sqrt(data$richness))
#hist(data$shannon)
#qqnorm((data$pilous_evenness))
#hist(data$Chao1)

extractionmodel <- lmer(sqrt(richness)~ exmet + (1|SampleID_nosize), data = data)
summary(extractionmodel)
r.squaredGLMM(extractionmodel)
plot(extractionmodel)
qqnorm(resid(extractionmodel)) # normal enough

ems <- emmeans(extractionmodel, specs = pairwise ~ exmet, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

t.test(sqrt(richness)~exmet, mu = 0, alt = "two.sided", conf = 0.95, var.eq = F, paired = T, data = data)

