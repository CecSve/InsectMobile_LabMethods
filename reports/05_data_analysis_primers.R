# data analysis primer performance Insektmobilen

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
library(phyloseq)
library(ampvis2)
library(nlme)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(cowplot)

### load data #################t
metadata <- read.table("cleaned-data/metadata_primer_comparison.txt", sep = "\t", header = TRUE)
asvs <- read.table("cleaned-data/asvtable_primer_comparison.txt", sep = "\t", header = TRUE)
setdiff(colnames(asvs), metadata$PCRID_primer)
setdiff(metadata$PCRID_primer, colnames(asvs))
table(metadata$primer)

taxonomy <- read.table("cleaned-data/taxonomy_primer_comparison.txt", sep = "\t", header = TRUE)
tax <- read.table("cleaned-data/tax_primer_comparison.txt", sep = "\t", header = TRUE)

# get stats for taxonomic assignment
taxonomy %>% group_by(order) %>% tally()
apply(taxonomy, 2, function(x){ length(which(!is.na(x))) })
1214/1339 
1088/1339 
871/1339 

### analysis ############################
asvs <- column_to_rownames(asvs, var = "asvid")
#pa_asvs <- decostand(asvs, "pa") # make otutable into presence absence

rowSums(asvs)
colSums(asvs)
min(rowSums(asvs))
min(colSums(asvs))

# community composition ----
bray.data<-vegdist(t(asvs), Type = "bray", binary = TRUE) 
meta_mds <- metaMDS(bray.data)
plot(meta_mds)

class(meta_mds)
methods(class="metaMDS")

Y <- scores(meta_mds, display="sites")
plot(Y, type="n")
text(Y[,1], Y[,2], rownames(Y), col="red")

Y <- data.frame(Y)
mdsdata <- Y %>% rownames_to_column(var = "PCRID_primer")
data <- merge(metadata, mdsdata, by = "PCRID_primer")

# add differences in distance column
data$distdif <- data$NMDS1 - data$NMDS2

Y$labels <- metadata$PCRID
Y$colour <- metadata$primer

custom_final19 = c( "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#771155","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C")

# plot all samples with normal confidence ellipse by extraction method (0.95)
ggplot(Y,
       aes(
         x = NMDS1,
         y = NMDS2,
         colour = colour,
         group = labels
       )) +
  geom_point(size = 5) +
  theme_minimal() + labs(colour = "Primer") + geom_line(aes(group = labels), colour= "darkgrey", alpha = 0.5, linetype = "dashed", size = 1, show.legend = F)+ geom_text_repel(aes(label = labels), show.legend = FALSE) + scale_colour_manual(values = c("#77AADD","#114477", "#1ED278"), labels=c("ART","INS", "FWH"))  + theme_classic()

### lmer on distdiff ################################
distmodel <- lmer(distdif~ primer + (1|PCRID), data = data)
summary(distmodel)
r.squaredGLMM(distmodel)

ems <- emmeans(distmodel, specs = pairwise ~ primer, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

ems_contrasts = pairs(ems) %>%
  as.data.frame()

ems_means <- summary(ems) %>%
  as.data.frame()

a<- ems_means %>% mutate(
  emmeans.primer = fct_relevel(
    emmeans.primer,
    "art",
    "fwh",
    "ins"
  )
) %>% ggplot(aes(x=emmeans.primer,y=emmeans.emmean,ymax=emmeans.upper.CL,ymin=emmeans.lower.CL,size=2))
#this defines the plot type
b<-a+geom_pointrange()
#this flips the co-ordinates so your x axis becomes your y and vice versa
c<-b+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)
#all of this gets rid of the grey grid and legends
e<-d+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
#this sets x and y axis titles
f<-e+ labs(x = "Primer pair\n", y = "", subtitle = "A: Community composition")
#this sets axis label size
g<-f+theme(axis.text.x = element_text(size = 12, colour = 'black')) +theme(axis.text.y = element_text(size = 12, colour = 'black', face = "italic"))
#this sets axis title size and there is your finished summary plot!
primerplot <- g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete() + theme(plot.subtitle = element_text(size = 20, face = "bold"))

save_plot("plots/emmeans_primermethod_communitycomp.png", primerplot, base_height = 4, base_width = 8)

# ampvis data load ----
# join otus and taxonomy
otuspa <- decostand(asvs, "pa") # make otutable into presence absence, not needed for rarecurve
otuspa <- rownames_to_column(otuspa, var = "otuid") 
taxonomy <- column_to_rownames(taxonomy, var = "asvid")
taxonomy$kingdom <- "Animalia"
taxonomy <- taxonomy %>% select(kingdom, phylum, order, family, genus, species)
taxonomy <- rownames_to_column(taxonomy, var = "otuid")
taxonomy$class <- NA
taxonomy <- taxonomy %>% select(otuid, kingdom, phylum, class, order, family, genus, species)
#otus <- rownames_to_column(otus, var = "otuid") 
otutable <- left_join(otuspa, taxonomy, by = "otuid") # choose otuspa if you want pa instead
# revert otuids back to rownames if needed
otuspa <- column_to_rownames(otuspa, var = "otuid")
#otus <- column_to_rownames(otus, var = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames

# metadata must have sample data in the first column and column classes matter 
str(metadata)
metadata <- column_to_rownames(metadata, var = "PCRID_primer")
metadata <- rownames_to_column(metadata, var = "otuid")

# load amp data
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

amp <- amp_load(otutable = otutable, metadata = metadata)

# heatmap vizulaization ----
# convert read abundances into percentages
#psn <- transform_sample_counts(ps, function(x) x/sum(x) * 100)

library(RColorBrewer)
(brewer.pal(7,"BuPu"))

# heatmap (normalise gives the proportion of the sample)
primer_heatmap<- amp_heatmap(
  data = amp,
  group_by = "SampleID",
  facet_by = "primer",
  normalise = T,
  tax_aggregate = "Family",
  tax_add = "Order",
  plot_values = F,
  plot_colorscale = "sqrt",
  color_vector = c("grey90", "#9EBCDA", "#8C6BB1", "#6E016B"),
  tax_empty = "remove",
  tax_show = 50,
) + guides(fill = guide_legend(title = "Relative abundance \n(% per sample)", reverse = T)) 

save_plot("plots/heatmap_primermethod_relabun.png", primer_heatmap, base_height = 12, base_width = 16)

# venn diagram: https://madsalbertsen.github.io/ampvis2/reference/amp_venn.html 
amp_venn(amp,
         group_by = "primer",
         cut_a = 0.1, # Abundance cutoff in percent. OTU's below this abundance are excluded from the analysis. (default: 0.1)
         cut_f = 1) # Frequency cutoff in percent. OTU's within the top cut_f of the reads are considered a "core" OTU. (default: 80) 

# richness analysis ----
# richness indices - calculated for pa data
amp_alphadiv <-
  amp_alphadiv(
    amp,
    measure = c("observed", "shannon", "simpson", "invsimpson"),
    richness = TRUE
  )

otus <- rownames_to_column(asvs, var = "otuid")
otutable <- left_join(otus, taxonomy, by = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)
amp_abun <- amp_load(otutable = otutable, metadata = metadata)

# rarefaction curve
rarecurve <-
  amp_rarecurve(
    amp_abun,
    facet_by = "primer",
    stepsize = 100,
    color_by = "SampleID",
    facet_scales = "free"
  )

custom_final28 = c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E", "#DDAA77")

rarecurve + scale_color_manual(values = custom_final28) + geom_line(size = 1.5, na.rm = TRUE) + scale_x_continuous(
  labels = scales::number,
  limits = c(0, 60000),
  breaks = seq(from = 0, to = 60000, by = 2000)
) + guides(colour=guide_legend(ncol=8)) + ylab("Number of observed ASVs") + ylim(0, 180) + theme_minimal()+ theme(legend.position = "bottom", axis.text.x = element_text(angle = 90), legend.title = element_blank()) + guides(colour = guide_legend(nrow = 2))

# plot richness indices on log10 transformed y-axis
amp_alphadiv %>% tidyr::gather("id", "Richness", 11:15) %>% ggplot(., aes(primer, Richness, fill = primer)) + geom_boxplot(show.legend = F) + facet_wrap(~ id, ncol = 5) + scale_y_log10() + scale_fill_manual(values = c("#EDF8FB", "#BFD3E6", "#9EBCDA")) + theme_pubclean() + labs(fill = "Primer choice") + xlab("Primer choice") + ylab("Richness (log-transformed)")

### summary statistics ####
group_by(amp_alphadiv, primer) %>%
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

taxonomy <- taxonomy %>% column_to_rownames(var = "otuid")
taxa <- taxonomy %>% select(family) # or choose other taxonomic levels (but remember morphology only has order level for all individuals)

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otuspa <- rownames_to_column(otuspa, var = "otuid")
taxa <- rownames_to_column(taxa, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

data <- left_join(otuspa, taxa, by = "otuid")
data <- column_to_rownames(data, var = "otuid")

# working with reshape2 to wrangle data from wide to long format
longdata <- melt(data) 
longdata <- longdata[longdata$value!=0,] # exclude perhaps, so zeros can be used as well for lmer
longdata <- aggregate(. ~ family + variable, data = longdata, FUN = sum) # if values in the insect family column and the sample (variable) column are identical, then sum up how many unique otus there were in the sample from that insect family (richness)
meta <- dplyr::rename(metadata,c('variable'='otuid'))
#meta <- metadata

str(longdata)
family_richness <- longdata %>% select(variable, value) %>% group_by(variable) %>% summarize(value = sum(value)) # now summarize the richness per sample 
plot_data <- left_join(meta, family_richness)
plot_data <- as_tibble(plot_data)
plot_data <- plot_data %>% dplyr::rename(sampleid = variable, asvrichness = value)

richness_extractionmethod <- plot_data %>% dplyr::group_by(primer) %>% dplyr::summarise(richness = sum(asvrichness)) # asvrichness by primer and insect family
plot_data %>% group_by(primer) %>% tally()

### lmer on richness data ####
library(nlme)
library(lme4)
library(lmerTest)

primermodel <- lmer(asvrichness~ primer + (1|SampleID), data = plot_data) # should the interaction be * r : ? And Bo wanted random effect of sample:family, but I have removed it since model does not converge, maybe due to something like this: https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer 
summary(primermodel)

ems <- emmeans(primermodel, specs = pairwise ~ primer, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

library(MuMIn)
r.squaredGLMM(primermodel)

ems_contrasts = pairs(ems) %>%
  as.data.frame()

ems_means <- summary(ems) %>%
  as.data.frame()

a<- ems_means %>% mutate(
  emmeans.Extraction_method = fct_relevel(
    emmeans.Extraction_method,
    "ethanol",
    "nondestructive"
  )
) %>% ggplot(aes(x=emmeans.Extraction_method,y=emmeans.emmean,ymax=emmeans.upper.CL,ymin=emmeans.lower.CL,size=2))
#this defines the plot type
b<-a+geom_pointrange()
#this flips the co-ordinates so your x axis becomes your y and vice versa
c<-b+coord_flip()+scale_size_area(max_size = 1.5)
#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)
#all of this gets rid of the grey grid and legends
e<-d+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
#this sets x and y axis titles
f<-e+ labs(x = "", y = "", subtitle = "B: ASV richness")
#this sets axis label size
g<-f+theme(axis.text.x = element_text(size = 12, colour = 'black')) +theme(axis.text.y = element_text(size = 12, colour = 'black', face = "italic"))
#this sets axis title size and there is your finished summary plot!
explot <- g+theme(axis.title.x = element_text(size = 15, colour = 'black'))+theme(axis.title.y = element_text(size = 15, colour = 'black')) + scale_x_discrete(labels = c(
  "ethanol" = "Ethanol",
  "nondestructive" = "ND lysis buffer"
)) + theme(plot.subtitle = element_text(size = 20, face = "bold"))

save_plot("plots/emmeans_extractionmethod_richness.png", explot, base_height = 4, base_width = 8)
