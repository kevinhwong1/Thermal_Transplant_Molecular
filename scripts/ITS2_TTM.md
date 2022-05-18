Thermal Transplant ITS2 Analysis
================
KH Wong

# **Setup**

Set up workspace, set options, and load required packages.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

# Load Packages

``` r
library("dplyr")
library("tidyverse") 
library("ggplot2")
library("janitor")
library("funrar")
```

# Load data

``` r
its2 <- read.delim("output/ITS2/its2_type_profiles/3_ThermalTrans_analysis_20211216T114913.profiles.absolute.abund_and_meta.txt", sep = "\t")
sample <- read.csv("data/ITS2/ITS2_Samplesheet.csv")
meta <- read.csv("data/ITS2/ITS2_meta.csv")

#write.csv(its2, "output/ITS2/its2_type_profiles/ITS2_Abs_Profiles_Edit.csv")

its2_edit <- read.csv("output/ITS2/its2_type_profiles/ITS2_Abs_Profiles_Edit.csv")
```

# Edit metadata

``` r
meta_2 <- merge(sample, meta, by = "Sample.ID")
meta_3 <- meta_2 %>% dplyr::select(Sample.ID, Coral.ID, History, LifeStage, Origin, Treatment, Transplant)
```

# Make dataframes for strain level

``` r
# Cleaning data frame
its2_edit_strain <- its2_edit[-c(1:6), ] #removing unncessary rows
its2_edit_strain <- its2_edit_strain[-c(1:2)] # Removing first column

its2_edit_strain_2 <- its2_edit_strain %>%
  row_to_names(row_number = 1) # Making a new column header

colnames(its2_edit_strain_2)[1] <- "Sample.ID" # renaming column 2 

its2_edit_strain_3 <- its2_edit_strain_2[-c(1)] # Removing first column
rownames(its2_edit_strain_3) <- its2_edit_strain_2[,1] #making the sample ids the row names 

its2_edit_strain_4 <- data.frame(sapply(its2_edit_strain_3, function(x) as.numeric(as.character(x)))) # making dataframe numeric
```

# Make relative abundance df at strain level and merging metadata

``` r
# Making the relative abundance dataframe 
its2_edit_strain_5 <- data.matrix(sapply(its2_edit_strain_4, as.numeric))
rel_its2_strain <- make_relative(its2_edit_strain_5)

#Transform Relative Abundance Matrix using sqrt
rel_its2_strain_sqrt <- sqrt(rel_its2_strain) #sqrt transform the matrix

rel_its2_strain_df <- as.data.frame(rel_its2_strain_sqrt)

# adding Fragment ID column
rel_its2_strain_df_2 <- cbind(rel_its2_strain_df, its2_edit_strain_2$Sample.ID) 
colnames(rel_its2_strain_df_2)[15] <- "Sample.ID" # renaming column

# Gathering data 
rel_its2_strain_df_3 <- rel_its2_strain_df_2 %>%
  gather(key=Strain, value=relabund, 1:14)

# Merging metadata
rel_its2_strain_meta <- merge(rel_its2_strain_df_3, meta_3, by = "Sample.ID")
rel_its2_strain_meta2 <- rel_its2_strain_meta %>% dplyr::select(-Sample.ID)

#Filtering for only samples transplanted to the Patch site and samples that match WGBS 
rel_its2_strain_meta_select <- rel_its2_strain_meta2 %>% 
  filter(Transplant == "Patch") %>%
  filter(Coral.ID != "P-3-A") %>%
  filter(Coral.ID != "P-16-A") %>%
  filter(Coral.ID != "R-14-A")

#Filtering for only strains present in this subset of data
rel_its2_strain_meta_select2 <- rel_its2_strain_meta_select %>% 
  filter(Strain != "B1374") %>%
  filter(Strain != "B2") %>%
  filter(Strain != "B23") %>%
  filter(Strain != "C17") %>%
  filter(Strain != "D6")

#Making relative abundance into percent
rel_its2_strain_meta_select2$relabund_per <- rel_its2_strain_meta_select2$relabund *100

#making strain a factor 
rel_its2_strain_meta_select2$Strain <- as.factor(rel_its2_strain_meta_select2$Strain)
```

# Plotting Relative Abundance

``` r
# library(RColorBrewer)
# library(colorRamps)
# library(ggthemes)
# library(scales)
# 
# rel_abun_bar <- ggplot(rel_its2_strain_meta_select2, aes(x = Coral.ID, y = relabund, fill = Strain)) + 
#   geom_bar(stat = "identity", aes(fill=Strain)) + 
#   labs(x="", y="Relative Abundance") +
#   scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(9))  +
# #  scale_fill_brewer(palette="Spectral") +
#   facet_wrap(~LifeStage+Origin+Treatment, scales= "free_x", nrow=1) +
#   theme_bw() + 
#   theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#         panel.grid.minor = element_blank(), axis.line = element_blank(), axis.text.x.bottom = element_text(angle = -45))
# rel_abun_bar
# 
# ggsave(filename="output/ITS2/Relative_Abundance_Barplot.png", plot=rel_abun_bar, dpi=300, width=10, height=5, units="in")
```

# Plotting Heatmap

``` r
relabun_heat<- rel_its2_strain_meta_select2 %>%
  mutate(Strain = fct_reorder(Strain,relabund_per, .desc=TRUE)) %>%
  ggplot( aes(x = Coral.ID, y = Strain)) + 
  geom_tile(aes(fill = relabund_per), color = "grey") +
  scale_fill_distiller(palette = "Reds", direction=100, labels = scales::label_percent(scale=1)) + 
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = -90),
        axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),) + 
  facet_grid(~LifeStage+Origin+Treatment, scales = "free") +
  labs(x="Sample ID", fill="Relative Abundance", y="Strain") 

relabun_heat
```

![](ITS2_TTM_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave(filename="output/ITS2/Relative_Abundance_Heatmap.png", plot=relabun_heat, dpi=300, width=10, height=5, units="in")
```

# PERMANOVA

``` r
library(tidyr)
library(reshape2)
library(ape)
library(vegan)

relabun_wide <- dcast(rel_its2_strain_meta_select2, Coral.ID+History+LifeStage+Origin+Treatment~Strain, value.var= "relabund") %>%
  drop_na()

# Bray Curtis caluclations 
dist <- vegdist(relabun_wide[c(6:ncol(relabun_wide))], method="bray")

# PERMANOVA with Bray Curtis
permanova_dist <- adonis(dist ~ Origin*Treatment*LifeStage, data = relabun_wide, method='eu', permutations=9999)
z_dist <- permanova_dist$aov.tab
z_dist
```

    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
    ## Origin                      1    0.4018 0.40179  1.6997 0.04877 0.1989   
    ## Treatment                   1    0.0148 0.01483  0.0628 0.00180 0.9675   
    ## LifeStage                   1    1.6995 1.69950  7.1895 0.20630 0.0058 **
    ## Origin:Treatment            1    0.0553 0.05535  0.2341 0.00672 0.7190   
    ## Origin:LifeStage            1    0.5079 0.50789  2.1485 0.06165 0.1247   
    ## Treatment:LifeStage         1    0.5538 0.55381  2.3428 0.06722 0.1020   
    ## Origin:Treatment:LifeStage  1    0.0409 0.04088  0.1729 0.00496 0.8263   
    ## Residuals                  21    4.9641 0.23639         0.60258          
    ## Total                      28    8.2382                 1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
capture.output(z_dist, file = "output/ITS2/PERMANOVA_ITS2_bray.csv")
```
