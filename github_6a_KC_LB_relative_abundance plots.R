
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(data.table)
library(RColorBrewer)

source("./Helper_Functions.R")

Laura.ps = readRDS("physeq1007")
Laura.ps <- Laura.ps %>% subset_taxa(Kingdom!= "Unassigned")

# isofemale line 10 (lab derived)
ps_10 <- prune_samples(sample_data(Laura.ps)$line == "10", Laura.ps)
ps_10 <- merge_samples(ps_10, "sampletype")
family_ps_10 <- MakeAbundanceDF(ps_10, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_10 <- subset(family_ps_10, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_10, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_10 <- rbind(family_ps_10,other.abund)
family_ps_10 %>% mutate(Line = "10", Location = "lab") -> family_ps_10

# make a function and do this for all
process_data <- function(line_num) {
  ps <- prune_samples(sample_data(Laura.ps)$line == line_num, Laura.ps)
  ps <- merge_samples(ps, "sampletype")
  ps <- MakeAbundanceDF(ps, taxRank = "Family", abundanceFilter = 0.05)
  ps <- subset(ps, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
  other.abund <- aggregate(Abundance ~ Sample, ps, sum)
  other.abund$Abundance <- 1- other.abund$Abundance
  other.abund$Phylum <- "Other"
  other.abund$Class <- "Other"
  other.abund$Order <- "Other"
  other.abund$Family <- "Other"
  ps <- rbind(ps,other.abund)
}

# run function for each isolfemale line 

# isofemale line 7 (lab derived)
family_ps_7 <- process_data("7")
family_ps_7 %>% mutate(Line = "7", Location = "lab") -> family_ps_7

# isofemale line 2 (field)
family_ps_2 <- process_data("2")
family_ps_2 %>% mutate(Line = "2", Location = "field") -> family_ps_2

# isofemale line 3 (field)
family_ps_3 <- process_data("3")
family_ps_3 %>% mutate(Line = "3", Location = "field") -> family_ps_3

# isofemale line 4 (field)
family_ps_4 <- process_data("4")
family_ps_4 %>% mutate(Line = "4", Location = "field") -> family_ps_4

# isofemale line 12 (lab)
family_ps_12 <- process_data("12")
family_ps_12 %>% mutate(Line = "12", Location = "lab") -> family_ps_12

# isofemale line 15 (lab)
family_ps_15 <- process_data("15")
family_ps_15 %>% mutate(Line = "15", Location = "lab") -> family_ps_15

# isofemale line 17 (lab)
family_ps_17 <- process_data("17")
family_ps_17 %>% mutate(Line = "17", Location = "lab") -> family_ps_17

# isofemale line 27 (lab)
family_ps_27 <- process_data("27")
family_ps_27 %>% mutate(Line = "27", Location = "lab") -> family_ps_27

# isofemale line 23 (lab)
family_ps_23 <- process_data("23")
family_ps_23 %>% mutate(Line = "23", Location = "lab") -> family_ps_23

# collate data for all isofemale line.
# next make on df that has all the info together, this should make plotting empty sections easier
df_family_abund_all <- rbind(family_ps_10, family_ps_7, family_ps_12, family_ps_15, family_ps_17, family_ps_2, family_ps_23, family_ps_27, family_ps_3, family_ps_4)
df_family_abund_all[is.na(df_family_abund_all)] <- "Other"

## plotting - all samples

# Assign colors
#myColors <- c("#A3DAC9", "#85CEB7", "#66C2A5", "#529B84", "#3D7463",
#	"#FDAF91", "#FC8D62", "#BBC6E0", "#A4B3D5", "#8DA0CB",
#	"#7180A2", "#55607A", "#384051", "#F1B9DB", "#ECA1CF",
#	"#E78AC3", "#B96E9C", "#8B5375", "#5C374E", "#DBEFBB",
#	"#CAE898", "#C1E487", "#B8E076", "#A6D854", "#85AD43",
#	"#74973B", "#648232", "#536C2A", "#425622", "#B3B3B3")

#names(myColors) <- c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
#	"Sphingobacteriaceae", "Weeksellaceae", "Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
#	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae", "Acetobacteraceae", "Beijerinckiaceae", 
#	"Caulobacteraceae", "Rhodobacteraceae", "Sphingomonadaceae", "Xanthobacteraceae", "Burkholderiaceae", 
#	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
#	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other")

myColors <- c("#A3DAC9", "#85CEB7", "#66C2A5", "#529B84", "#3D7463",
	"#BBC6E0", "#A4B3D5", "#8DA0CB",
	"#7180A2", "#55607A", "#384051",
	"#B96E9C", "#5C374E", "#DBEFBB",
	"#CAE898", "#C1E487", "#B8E076", "#A6D854", "#85AD43",
	"#74973B", "#648232", "#536C2A", "#425622", "#B3B3B3")

names(myColors) <- c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other")

# circle plot - offspring - larval water

Larval_water_abund <- df_family_abund_all[df_family_abund_all$Sample == "Larval Water",]
Larval_water_abund$Line <- factor(Larval_water_abund$Line, levels = c("7","10","12","15","17","23","27","2","3","4"))
Larval_water_abund$Family <- factor(Larval_water_abund$Family, levels=c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other"))

Larval_water_circle <- ggplot(Larval_water_abund,
                        aes(x = Line, y = Abundance, fill = Family)) +
  geom_col() +
  scale_fill_manual(values = myColors) +
  scale_y_continuous("", limits = c(-1, 1.1)) +
  coord_polar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())
Larval_water_circle

# circle plot - offspring - larvae

Larvae_abund <- df_family_abund_all[df_family_abund_all$Sample == "Offspring - Larvae",]
Larvae_abund$Line <- factor(Larvae_abund$Line, levels = c("7","10","12","15","17","23","27","2","3","4"))
Larvae_abund$Family <- factor(Larvae_abund$Family, levels=c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other"))

Larvae_circle <- ggplot(Larvae_abund,
                        aes(x = Line, y = Abundance, fill = Family)) +
  geom_col() +
  scale_fill_manual(values = myColors) +
  scale_y_continuous("", limits = c(-1, 1.1)) +
  coord_polar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())
Larvae_circle

# circle plot - offspring - adults

Adult_abund <- df_family_abund_all[df_family_abund_all$Sample == "Offspring - Adult",]
Adult_abund$Line <- factor(Adult_abund$Line, levels = c("7","10","12","15","17","23","27","2","3","4"))
Adult_abund$Family <- factor(Adult_abund$Family, levels=c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other"))

Adult_circle <- ggplot(Adult_abund,
                        aes(x = Line, y = Abundance, fill = Family)) +
  geom_col() +
  scale_fill_manual(values = myColors) +
  scale_y_continuous("", limits = c(-1, 1.1)) +
  coord_polar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())
Adult_circle

# mom ovary plots (non-circular)

Ovary_abund <- df_family_abund_all[df_family_abund_all$Sample == "Ovary",]
Ovary_abund$Line <- factor(Ovary_abund$Line, levels = c("7","10","12","15","17","23","27","2","3","4"))
Ovary_abund$Family <- factor(Ovary_abund$Family, levels=c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other"))

Ovary_linear <- ggplot(Ovary_abund,
                        aes(x = Line, y = Abundance, fill = Family)) +
  geom_col() +
  scale_fill_manual(values = myColors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())
Ovary_linear

# mom midgut plots (non-circular)

Midgut_abund <- df_family_abund_all[df_family_abund_all$Sample == "Midgut",]
Midgut_abund$Line <- factor(Midgut_abund$Line, levels = c("7","10","12","15","17","23","27","2","3","4"))
Midgut_abund$Family <- factor(Midgut_abund$Family, levels=c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other"))

Midgut_linear <- ggplot(Midgut_abund,
                        aes(x = Line, y = Abundance, fill = Family)) +
  geom_col() +
  scale_fill_manual(values = myColors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())
Midgut_linear


## circle plot - mothers

Laura.ps.mom <- subset_samples(Laura.ps, sampletype=="Ovary" | sampletype=="Midgut")

# isofemale line 10 (lab-derived)
ps_10_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "10", Laura.ps.mom)
family_ps_10_mom <- MakeAbundanceDF(ps_10_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_10_mom <- subset(family_ps_10_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_10_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_10_mom <- rbind(family_ps_10_mom,other.abund)
family_ps_10_mom %>% mutate(Line = "10", Location = "lab") -> family_ps_10_mom

# isofemale line 7 (lab-derived)
family_ps_7_mom <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(family_ps_7_mom) <- c("Sample", "Abundance", "Phylum", "Class", "Order", "Family")
family_ps_7_mom[1,1] <- 7
family_ps_7_mom[1,2] <- 0
family_ps_7_mom[1,3:6] <- "Other"
family_ps_7_mom %>% mutate(Line = "7", Location = "lab") -> family_ps_7_mom

# isofemale line 2 (field)
ps_2_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "2", Laura.ps.mom)
family_ps_2_mom <- MakeAbundanceDF(ps_2_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_2_mom <- subset(family_ps_2_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_2_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_2_mom <- rbind(family_ps_2_mom,other.abund)
family_ps_2_mom %>% mutate(Line = "2", Location = "lab") -> family_ps_2_mom

# isofemale line 3 (field)
ps_3_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "3", Laura.ps.mom)
ps_3_mom <- merge_samples(ps_3_mom, "line")
family_ps_3_mom <- MakeAbundanceDF(ps_3_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_3_mom <- subset(family_ps_3_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_3_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_3_mom <- rbind(family_ps_3_mom,other.abund)
family_ps_3_mom %>% mutate(Line = "3", Location = "lab") -> family_ps_3_mom

# isofemale line 4 (field)
ps_4_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "4", Laura.ps.mom)
ps_4_mom <- merge_samples(ps_4_mom, "line")
family_ps_4_mom <- MakeAbundanceDF(ps_4_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_4_mom <- subset(family_ps_4_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_4_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_4_mom <- rbind(family_ps_4_mom,other.abund)
family_ps_4_mom %>% mutate(Line = "4", Location = "lab") -> family_ps_4_mom

# isofemale line 12 (lab)
ps_12_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "12", Laura.ps.mom)
family_ps_12_mom <- MakeAbundanceDF(ps_12_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_12_mom <- subset(family_ps_12_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_12_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_12_mom <- rbind(family_ps_12_mom,other.abund)
family_ps_12_mom %>% mutate(Line = "12", Location = "lab") -> family_ps_12_mom

# isofemale line 15 (lab)
ps_15_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "15", Laura.ps.mom)
ps_15_mom <- merge_samples(ps_15_mom, "line")
family_ps_15_mom <- MakeAbundanceDF(ps_15_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_15_mom <- subset(family_ps_15_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_15_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_15_mom <- rbind(family_ps_15_mom,other.abund)
family_ps_15_mom %>% mutate(Line = "15", Location = "lab") -> family_ps_15_mom

# isofemale line 17 (lab)
ps_17_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "17", Laura.ps.mom)
family_ps_17_mom <- MakeAbundanceDF(ps_17_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_17_mom <- subset(family_ps_17_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_17_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_17_mom <- rbind(family_ps_17_mom,other.abund)
family_ps_17_mom %>% mutate(Line = "17", Location = "lab") -> family_ps_17_mom

# isofemale line 27 (lab)
ps_27_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "27", Laura.ps.mom)
ps_27_mom <- merge_samples(ps_27_mom, "line")
family_ps_27_mom <- MakeAbundanceDF(ps_27_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_27_mom <- subset(family_ps_27_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_27_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_27_mom <- rbind(family_ps_27_mom,other.abund)
family_ps_27_mom %>% mutate(Line = "27", Location = "lab") -> family_ps_27_mom

# isofemale line 23 (lab)
ps_23_mom <- prune_samples(sample_data(Laura.ps.mom)$line == "23", Laura.ps.mom)
ps_23_mom <- merge_samples(ps_23_mom, "line")
family_ps_23_mom <- MakeAbundanceDF(ps_23_mom, 
                                   taxRank = "Family", 
                                   abundanceFilter = 0.05)
family_ps_23_mom <- subset(family_ps_23_mom, select = c("Sample", "Abundance","Phylum","Class","Order","Family"))
other.abund <- aggregate(Abundance ~ Sample, family_ps_23_mom, sum)
other.abund$Abundance <- 1- other.abund$Abundance
other.abund$Phylum <- "Other"
other.abund$Class <- "Other"
other.abund$Order <- "Other"
other.abund$Family <- "Other"
family_ps_23_mom <- rbind(family_ps_23_mom,other.abund)
family_ps_23_mom %>% mutate(Line = "23", Location = "lab") -> family_ps_23_mom

# collate data for all isofemale lines
# next make on df that has all the info together, this should make plotting empty sections easier
df_family_abund_all_mom <- rbind(family_ps_10_mom, family_ps_12_mom, family_ps_15_mom, family_ps_17_mom, family_ps_2_mom, family_ps_7_mom, family_ps_23_mom, family_ps_27_mom, family_ps_3_mom, family_ps_4_mom)
df_family_abund_all_mom[is.na(df_family_abund_all_mom)] <- "Other"

## plotting - all samples
df_family_abund_all_mom$Line <- factor(df_family_abund_all_mom$Line, levels = c("7","10","12","15","17","23","27","2","3","4"))
df_family_abund_all_mom$Family <- factor(df_family_abund_all_mom$Family, levels=c("Corynebacteriaceae", "Dermacoccaceae", "Microbacteriaceae", "Micrococcaceae", "Promicromonosporaceae", 
	"Bacillaceae", "Brevibacillaceae", "Paenibacillaceae", 
	"Peptostreptococcaceae", "Staphylococcaceae", "Streptococcaceae",
	"Rhodobacteraceae", "Xanthobacteraceae", "Burkholderiaceae", 
	"Comamonadaceae", "Enterobacteriaceae", "Halomonadaceae", "Moraxellaceae", "Pseudomonadaceae", 
	"Rhodanobacteraceae", "Rhodocyclaceae", "Xanthomonadaceae", "Yersiniaceae", "Other"))

Mom_circle <- ggplot(df_family_abund_all_mom,
                        aes(x = Line, y = Abundance, fill = Family)) +
  geom_col() +
  scale_fill_manual(values = myColors) +
  scale_y_continuous("", limits = c(-1, 1.1)) +
  coord_polar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())
Mom_circle
