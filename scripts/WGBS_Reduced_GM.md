Thermal Transplant WGBS Reduced Data DMG Analysis
================
KH Wong

## Load Packages

<https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/GM.Rmd>

``` r
library(plotrix) 
library(ggplot2)
library(gridExtra)
library(seacarb) 
library(pheatmap)
library(tidyverse)
library(tidyr)
library(goseq)
library(genefilter)
library(gplots)
library(cowplot)
library(lsmeans)
library(data.table)
library(RColorBrewer)
library(randomcoloR)
#library(GSEABase)
library(ggpubr)
library(hrbrthemes)
library(viridis)
library(factoextra)
library(ropls)
library(mixOmics)
library(vegan)
library(plyr)
library(dplyr)
```

# Loading annotation file

``` r
Annot <- read.csv("../data/Past_genome/Past_annotation_20220310.csv", header=TRUE, na.string="NA", stringsAsFactors = FALSE, skip=0)

Annot[Annot == ""] <- NA

Annot.1 <- Annot %>%
  dplyr::select(gene, Length, InterPro_GO_Names, BLAST_GO_Names, SwissProt_GO_Names, TrEMBL_GO_Names) %>%
  mutate(InterPro_GO_Names = na_if(InterPro_GO_Names, "no GO terms")) %>%
  mutate(InterPro_GO_Names = na_if(InterPro_GO_Names, "no IPS match")) %>%
  unite("GO_Name", InterPro_GO_Names:TrEMBL_GO_Names, remove = FALSE, na.rm = TRUE) %>%
  dplyr::select(gene, Length, GO_Name)

Annot.2 <- Annot %>%
  dplyr::select(gene, InterPro_GO_IDs, BLAST_GO_IDs, SwissProt_GO_IDs, TrEMBL_GO_IDs) %>%
  mutate(InterPro_GO_IDs = na_if(InterPro_GO_IDs, "no GO terms")) %>%
  mutate(InterPro_GO_IDs = na_if(InterPro_GO_IDs, "no IPS match")) %>%
  unite("GO_IDs", InterPro_GO_IDs:TrEMBL_GO_IDs, remove = FALSE, na.rm = TRUE) %>%
 dplyr::select(gene, GO_IDs)

Annot.select <- merge(Annot.1, Annot.2, by = "gene")
Annot.select$GO_Name <- gsub("F:*","",Annot.select$GO_Name)
Annot.select$GO_Name <- gsub("P:*","",Annot.select$GO_Name)
Annot.select$GO_Name <- gsub("C:*","",Annot.select$GO_Name)
Annot.select$GO_Name <- gsub("_","; ",Annot.select$GO_Name)

Annot.select$GO_IDs <- gsub("F:*","",Annot.select$GO_IDs)
Annot.select$GO_IDs <- gsub("P:*","",Annot.select$GO_IDs)
Annot.select$GO_IDs <- gsub("C:*","",Annot.select$GO_IDs)
Annot.select$GO_IDs <- gsub("_","; ",Annot.select$GO_IDs)



#### ANNOTATION STATSITICS 

Annot.NA <- Annot %>%
  dplyr::select(gene, Length, InterPro_GO_Names, BLAST_GO_Names, SwissProt_GO_Names, TrEMBL_GO_Names) %>%
  mutate(InterPro_GO_Names = na_if(InterPro_GO_Names, "no GO terms")) %>%
  mutate(InterPro_GO_Names = na_if(InterPro_GO_Names, "no IPS match")) %>%
  mutate_all(na_if,"")

Annot.IPS <- Annot.NA %>% drop_na(InterPro_GO_Names) # 24089
Annot.BLAST <- Annot.NA %>% drop_na(BLAST_GO_Names) # 11443
Annot.SP <- Annot.NA %>% drop_na(SwissProt_GO_Names) # 30284
Annot.Trem <- Annot.NA %>% drop_na(TrEMBL_GO_Names) # 12178

Annot.all <- Annot.select %>% mutate_all(na_if,"") %>% drop_na() #43154 
```

# Load methylation data

``` r
# All Reduced Data
meth.data_5x <- list.files(path = "../output/WGBS/cov_to_cyto_reduced", pattern = ".bed$", full.names=TRUE) %>%
  set_names(.) %>% 
  map_dfr(read.csv,.id="Sample.ID", header=FALSE, sep="\t", na.string="NA", stringsAsFactors = FALSE) %>% 
  dplyr::select(-c(V3,V7:V14)) %>%
  group_by(Sample.ID)
colnames(meth.data_5x) <- c("Sample.ID", "scaffold", "position","per.meth","meth","unmeth","gene")
meth.data_5x$gene <- gsub(";.*","",meth.data_5x$gene) #remove extra characters
meth.data_5x$gene <- gsub("ID=","",meth.data_5x$gene) #remove extra characters
meth.data_5x$Sample.ID <- gsub("../output/WGBS/cov_to_cyto_reduced/","",meth.data_5x$Sample.ID) #remove extra characters
meth.data_5x$Sample.ID <- gsub("*._5x_sorted.tab_gene_CpG_5x_enrichment.bed","",meth.data_5x$Sample.ID) #remove extra characters 
meth.data_5x$Sample.ID <- gsub("_...","",meth.data_5x$Sample.ID) #remove extra characters

sample.info <- read.csv("../data/metadata/Thermal_Transplant_Metadata_Analysis.csv", stringsAsFactors = TRUE)
sample.info2 <- sample.info %>%
  dplyr::select(-Timepoint, -File_Name)

MD_5x <- merge(meth.data_5x, sample.info2, by="Sample.ID")
```
