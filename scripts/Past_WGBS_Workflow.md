# *Porites astreoides* Thermal Transplant WGBS Workflow

This script takes raw whole genome bisulfite sequencing files to bed files where they can then be analzyed in R. The laboratory preparation methods can be found [here](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

**Table of Contents:**
1. [**Methylation Quantification with Methylseq**](#methylseq)  
2. [**Merge strands with Biskmark**](#merge)  
3. [**Create bedgraphs post merge**](#bedgraph)   
4. [**CpGs for 5x coverage**](#5xcov)    
5. [**CpGs for 10x coverage**](#10xcov)    
6. [**Select samples**](#select)    
7. [**Intersect across samples and genes**](#intersect)    


## <a name="methylseq"></a> **1. Methylation Quantification with Methylseq**

A few testing parameters for the settings can be found in this [notebook post](https://kevinhwong1.github.io/KevinHWong_Notebook/Methylseq-trimming-test-to-remove-m-bias/). I decided to use the Trim_3 parameters as it keeps the most amount of data while still reducing the m-bias.

I also used the *Porites astreoides* genome from [this repository](https://github.com/hputnam/Past_Genome).

`nano methylseq_trim3.sh`

```bash
#!/bin/bash
#SBATCH --job-name="methylseq"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3
#SBATCH --exclusive

# load modules needed

module load Nextflow/21.03.0

# run nextflow methylseq

nextflow run nf-core/methylseq \
-profile singularity \
--aligner bismark \
--igenomes_ignore \
--fasta /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--save_reference \
--input '/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 15 \
--clip_r2 30 \
--three_prime_clip_r1 30 \
--three_prime_clip_r2 15 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir WGBS_methylseq
```

```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/MultiQC/multiqc_report.html MyProjects/Thermal_Transplant_Molecular/output/multiqc_report_trim3_full.html
```

## <a name="merge"></a> **2. Merge strands with Biskmark**

The Bismark coverage2cytosine command re-reads the genome-wide report and merges methylation evidence of both top and bottom strand.

`mkdir cov_to_cyto`

`nano cov_to_cyto.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --exclusive

# load modules needed

module load Bismark/0.20.1-foss-2018b

# run coverage2cytosine merge of strands

 find /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/methylseq_trim3/WGBS_methylseq/bismark_methylation_calls/methylation_coverage/*deduplicated.bismark.cov.gz \
 | xargs basename -s _L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz \
 | xargs -I{} coverage2cytosine \
 --genome_folder /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/methylseq_trim3/WGBS_methylseq/reference_genome/BismarkIndex \
 -o {} \
 --merge_CpG \
 --zero_based \
/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/methylseq_trim3/WGBS_methylseq/bismark_methylation_calls/methylation_coverage/{}_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
```

We now have to sort the merged files so the scaffolds are all in the same order and multiIntersectBed will run correctly. Run this for loop using bedtools to sort all .tab files

`nano bedtools.sort.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

module load BEDTools/2.27.1-foss-2018b

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  bedtools sort -i "${f}" \
  > "${STEM}"_sorted.cov
done
```

## <a name="bedgraph"></a> **3. Create bedgraphs post merge**

These bedgraph files are neccessary for the characterization of methylation in [genomic feature analysis](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Past_Genomic_Feature_Analysis_20221014.md).

**5x bedgraph**

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x_sorted.bedgraph
done
```

**10x bedgraph**

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x_sorted.bedgraph
done
```


## <a name="5xcov"></a> **4. CpGs for 5x coverage**

Run loop to filter CpGs for 5x coverage, creating tab files with raw count for glms.

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_5x_sorted.tab
done
```

`wc -l *5x_sorted.tab`

```
4527111 18-106_S163_5x_sorted.tab
5597804 18-118_S162_5x_sorted.tab
4104763 18-130_S172_5x_sorted.tab
4884343 18-142_S189_5x_sorted.tab
2463417 18-167_S166_5x_sorted.tab
5188393 18-178_S191_5x_sorted.tab
5235506 18-190_S186_5x_sorted.tab
3310274 18-202_S188_5x_sorted.tab
4670049 18-20_S202_5x_sorted.tab
4704858 18-227_S170_5x_sorted.tab
3954330 18-239_S185_5x_sorted.tab
6100768 18-250_S195_5x_sorted.tab
6854723 18-262_S179_5x_sorted.tab
5632004 18-311_S187_5x_sorted.tab
5462756 18-322_S180_5x_sorted.tab
5424011 18-32_S178_5x_sorted.tab
 377146 18-334_S164_5x_sorted.tab
7094939 18-346_S193_5x_sorted.tab
4928861 18-358_S201_5x_sorted.tab
5907039 18-370_S171_5x_sorted.tab
4421506 18-394_S192_5x_sorted.tab
5986501 18-406_S177_5x_sorted.tab
5640553 18-418_S196_5x_sorted.tab
7830166 18-442_S165_5x_sorted.tab
6501465 18-44_S198_5x_sorted.tab
7627601 18-454_S197_5x_sorted.tab
7787309 18-466_S199_5x_sorted.tab
5835415 18-55_S190_5x_sorted.tab
4083785 18-67_S176_5x_sorted.tab
7046855 18-79_S181_5x_sorted.tab
1035017 18-91_S160_5x_sorted.tab
 725233 18-9_S159_5x_sorted.tab
4867549 L-1029_S183_5x_sorted.tab
4555292 L-1038_S184_5x_sorted.tab
6352836 L-1053_S167_5x_sorted.tab
7551302 L-1059_S175_5x_sorted.tab
7108688 L-1093_S168_5x_sorted.tab
8781711 L-1257_S205_5x_sorted.tab
8085738 L-1263_S173_5x_sorted.tab
4471821 L-562_S174_5x_sorted.tab
7147304 L-571_S194_5x_sorted.tab
6156923 L-661_S182_5x_sorted.tab
6763146 L-704_S169_5x_sorted.tab
2165197 L-728_S161_5x_sorted.tab
4272213 L-862_S200_5x_sorted.tab
5074266 L-924_S204_5x_sorted.tab
6627315 L-933_S203_5x_sorted.tab
250925802 total
```

## <a name="10xcov"></a> **5. CpGs for 10x coverage**

Run loop to filter CpGs for 10x coverage, creating tab files with raw count for glms.

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_10x_sorted.tab
done
```

`wc -l *10x_sorted.tab`

```
676580 18-106_S163_10x_sorted.tab
1118168 18-118_S162_10x_sorted.tab
847700 18-130_S172_10x_sorted.tab
1142935 18-142_S189_10x_sorted.tab
210115 18-167_S166_10x_sorted.tab
1315924 18-178_S191_10x_sorted.tab
1055468 18-190_S186_10x_sorted.tab
411100 18-202_S188_10x_sorted.tab
825928 18-20_S202_10x_sorted.tab
1134939 18-227_S170_10x_sorted.tab
696483 18-239_S185_10x_sorted.tab
1787233 18-250_S195_10x_sorted.tab
1943403 18-262_S179_10x_sorted.tab
1543061 18-311_S187_10x_sorted.tab
1427802 18-322_S180_10x_sorted.tab
1068870 18-32_S178_10x_sorted.tab
  8494 18-334_S164_10x_sorted.tab
2079079 18-346_S193_10x_sorted.tab
1154658 18-358_S201_10x_sorted.tab
1390085 18-370_S171_10x_sorted.tab
865273 18-394_S192_10x_sorted.tab
1688362 18-406_S177_10x_sorted.tab
1464387 18-418_S196_10x_sorted.tab
2583003 18-442_S165_10x_sorted.tab
1917367 18-44_S198_10x_sorted.tab
2678858 18-454_S197_10x_sorted.tab
2749194 18-466_S199_10x_sorted.tab
1629133 18-55_S190_10x_sorted.tab
597189 18-67_S176_10x_sorted.tab
2694546 18-79_S181_10x_sorted.tab
 53547 18-91_S160_10x_sorted.tab
 27615 18-9_S159_10x_sorted.tab
976316 L-1029_S183_10x_sorted.tab
900254 L-1038_S184_10x_sorted.tab
1551870 L-1053_S167_10x_sorted.tab
3312304 L-1059_S175_10x_sorted.tab
2415616 L-1093_S168_10x_sorted.tab
3295825 L-1257_S205_10x_sorted.tab
3449982 L-1263_S173_10x_sorted.tab
1051976 L-562_S174_10x_sorted.tab
2590771 L-571_S194_10x_sorted.tab
1920044 L-661_S182_10x_sorted.tab
1790820 L-704_S169_10x_sorted.tab
205111 L-728_S161_10x_sorted.tab
642700 L-862_S200_10x_sorted.tab
1495979 L-924_S204_10x_sorted.tab
2332830 L-933_S203_10x_sorted.tab
68718897 total
```

## <a name="select"></a> **6. Select samples**

At this step, you have to make sure you are comparing all of the necessary samples in your analysis. I am only including a subset, therefore I have made a new directory with samples only pertaining to adult-larval pairs.

`mkdir cov_to_cyto_lifestage`

`cd cov_to_cyto_lifestage`

#### 5x coverage

```bash
cp ../cov_to_cyto_reduced/18-118_S162_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-202_S188_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-322_S180_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-358_S201_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-106_S163_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-190_S186_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-370_S171_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-454_S197_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-130_S172_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-142_S189_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-418_S196_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-55_S190_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-167_S166_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-227_S170_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-406_S177_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1029_S183_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-728_S161_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-924_S204_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-933_S203_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1053_S167_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1257_S205_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-704_S169_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-862_S200_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1038_S184_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1263_S173_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-562_S174_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-661_S182_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1059_S175_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1093_S168_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-571_S194_5x_sorted.tab ./
```

### 10x coverage

```bash
cp ../cov_to_cyto_reduced/18-118_S162_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-202_S188_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-322_S180_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-358_S201_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-106_S163_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-190_S186_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-370_S171_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-454_S197_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-130_S172_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-142_S189_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-418_S196_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-55_S190_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-167_S166_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-227_S170_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-406_S177_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1029_S183_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-728_S161_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-924_S204_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-933_S203_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1053_S167_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1257_S205_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-704_S169_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-862_S200_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1038_S184_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1263_S173_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-562_S174_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-661_S182_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1059_S175_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1093_S168_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-571_S194_10x_sorted.tab ./
```

## <a name="intersect"></a>**7. Intersect across samples and genes **

This script has multiple steps:
* Intersect positions found in all samples at a specific coverage (5x and 10x)
* Use `intersectBed` to find where loci and genes intersect, allowing loci to be mapped to annotated genes
* Intersect with file to subset only those positions found in all samples

`nano runall_meth_5x.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

#Intersect positions found in all samples at a specific coverage

multiIntersectBed -i *_5x_sorted.tab > CpG.lifestage.samps.5x_sorted.bed

cat CpG.lifestage.samps.5x_sorted.bed | awk '$4 ==30' > CpG.filt.lifestage.samps.5x_sorted.bed

# Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes
#wb: Print all lines in the second file
#a: file
#b: annotated gene list
#Save output in a new file that has the same base name

for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done

#Intersect with file to subset only those positions found in all samples

for i in *_5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.lifestage.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done
```

`wc -l *5x_enrichment.bed`

```
103636 18-106_S163_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-118_S162_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-130_S172_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-142_S189_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-167_S166_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-190_S186_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-202_S188_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-227_S170_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-322_S180_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-358_S201_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-370_S171_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-406_S177_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-418_S196_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-454_S197_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-55_S190_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1029_S183_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1038_S184_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1053_S167_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1059_S175_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1093_S168_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1257_S205_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1263_S173_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-562_S174_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-571_S194_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-661_S182_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-704_S169_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-728_S161_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-862_S200_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-924_S204_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-933_S203_5x_sorted.tab_gene_CpG_5x_enrichment.bed
 3109080 total
```

`nano runall_meth_10x.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

#Intersect positions found in all samples at a specific coverage

multiIntersectBed -i *_10x_sorted.tab > CpG.lifestage.samps.10x_sorted.bed

cat CpG.lifestage.samps.10x_sorted.bed | awk '$4 ==30' > CpG.filt.lifestage.samps.10x_sorted.bed

# Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes
#wb: Print all lines in the second file
#a: file
#b: annotated gene list
#Save output in a new file that has the same base name

for i in *10x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done

#Intersect with file to subset only those positions found in all samples

for i in *_10x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.lifestage.samps.10x_sorted.bed \
  > ${i}_CpG_10x_enrichment.bed
done
```

`wc -l *10x_enrichment.bed`

```
2687 18-106_S163_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-118_S162_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-130_S172_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-142_S189_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-167_S166_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-190_S186_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-202_S188_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-227_S170_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-322_S180_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-358_S201_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-370_S171_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-406_S177_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-418_S196_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-454_S197_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-55_S190_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1029_S183_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1038_S184_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1053_S167_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1059_S175_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1093_S168_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1257_S205_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1263_S173_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-562_S174_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-571_S194_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-661_S182_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-704_S169_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-728_S161_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-862_S200_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-924_S204_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-933_S203_10x_sorted.tab_gene_CpG_10x_enrichment.bed
80610 total
```

**Export bedfiles to local computer:**

```
scp 'kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage/*_5x_enrichment.bed' ~/MyProjects/Thermal_Transplant_Molecular/output/WGBS/meth_counts_5x

scp 'kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage/*_10x_enrichment.bed' ~/MyProjects/Thermal_Transplant_Molecular/output/WGBS/meth_counts_10x
```
