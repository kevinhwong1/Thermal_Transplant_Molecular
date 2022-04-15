## Genomic Feature Analysis

Following [this pipeline](https://github.com/hputnam/Geoduck_Meth/blob/master/code/Geoduck_Meth_Pipeline.md) created by H. Putnam.

### 1. Prepare Reference File

#### Creating a reference fasta file

```
module load SAMtools/1.9-foss-2018b

cd /data/putnamlab/kevin_wong1/Past_Genome/
samtools faidx past_filtered_assembly.fasta
```

This should create a fasta index file called `past_filtered_assembly.fasta.fai`

#### Making a putative promoter track

`nano put_promoter.sh`

```bash

#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

flankBed \
-i /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_noseq_v1.gff \
-g /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta.fai \
-l 1000 \
-r 0 \
-s |\
awk '{ gsub("gene","put_promoter",$3); print $0 }'|\
awk '{if($5-$4 > 3)print $0}'|\
tr ' ' '\t' \
> Pastreoides-v1.putprom.gff
```

#### Sort concatenated gff

```bash
interactive
module load BEDTools/2.27.1-foss-2018b

cat Pastreoides-v1.putprom.gff |\
awk -F"\t" '{print $1"\t"$4"\t"$5"\t"$3}'|\
sortBed \
-i - \
> Pastreoides-v1.sort.bed

exit
```

#### Bin features


```bash
interactive
module load BEDTools/2.27.1-foss-2018b

cat Pastreoides-v1.sort.bed |\
uniq |\
windowMaker \
-b - \
-w 2000 \
-i src | \
sortBed \
-i - \
> Pastreoides-v1.2Kbin.sort.bed

exit
```

#### Remove any duplicate lines from sorted binned concatenated gff

`cat Pastreoides-v1.2Kbin.sort.bed | uniq > Pastreoides-v1.2Kbin.uniq.bed `

note: some lines were removed, however under visual inspection, it seems like there are still duplicates. Revisit this later!

```
[kevin_wong1@n065 past_struc_annotations_v1]$ wc -l Pastreoides-v1.2Kbin.sort.bed
3684541 Pastreoides-v1.2Kbin.sort.bed
[kevin_wong1@n065 past_struc_annotations_v1]$ wc -l Pastreoides-v1.2Kbin.uniq.bed
3498993 Pastreoides-v1.2Kbin.uniq.bed
```

### 2. Determine Background Regions

#### Find all CpGs covered by at least 3/4 samples per group (use merge.tab files)

`nano sample_group_cov.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/genome_feature
#SBATCH --cpus-per-task=3

#first create variables for each file grouping
declare -a A_PAP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-118_S162_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-202_S188_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-322_S180_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-358_S201_5x_sorted.tab"

declare -a A_PAR="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-311_S187_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-346_S193_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-394_S192_5x_sorted.tab"

declare -a A_PHP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-106_S163_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-190_S186_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-370_S171_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-454_S197_5x_sorted.tab"

declare -a A_PHR="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-106_S163_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-190_S186_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-370_S171_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-454_S197_5x_sorted.tab"

declare -a A_PHR="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-262_S179_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-32_S178_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-67_S176_5x_sorted.tab"

declare -a A_RAP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-130_S172_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-142_S189_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-418_S196_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-55_S190_5x_sorted.tab"

declare -a A_RAR="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-178_S191_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-20_S202_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-239_S185_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-442_S165_5x_sorted.tab"

declare -a A_RHP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-167_S166_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-227_S170_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-406_S177_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-79_S181_5x_sorted.tab"

declare -a A_RHR="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-250_S195_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-44_S198_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/18-466_S199_5x_sorted.tab"

declare -a L_PAP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1029_S183_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-728_S161_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-924_S204_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-933_S203_5x_sorted.tab"

declare -a L_PHP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1053_S167_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1257_S205_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-704_S169_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-862_S200_5x_sorted.tab"

declare -a L_RAP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1038_S184_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1263_S173_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-562_S174_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-661_S182_5x_sorted.tab"

declare -a L_RHP="/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1059_S175_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-1093_S168_5x_sorted.tab /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced/L-571_S194_5x_sorted.tab"

#next create another list of variable names

groups=(A_PAP A_PAR A_PHP A_PHR A_RAP A_RAR A_RHP A_RHR L_PAP L_PHP L_RAP L_RHP)

#next loop through variable names to combine the files within each experimental group
#cat all group files together
#if the context is CG print the chromosome and position
#sort data
#unique and count the frequency of positions
#if the position is overlapping in 3/4 samples keep

for f in "${groups[@]}"
do
x="${f}"
cat ${!x} | \
awk '{print $1,$2,$3}' | \
sort | \
uniq -c | \
awk '{if($1>2)print $2"\t"$3"\t"$4}' \
> ${f}.3xCpG.bed
done
```

Checking the lengths of each BED group file
```
wc -l A*
  3387450 A_PAP.3xCpG.bed
  2733174 A_PAR.3xCpG.bed
  4540594 A_PHP.3xCpG.bed
  2479464 A_PHR.3xCpG.bed
  3829737 A_RAP.3xCpG.bed
  3822515 A_RAR.3xCpG.bed
  3716879 A_RHP.3xCpG.bed
  3765726 A_RHR.3xCpG.bed

wc -l L*
  3286169 L_PAP.3xCpG.bed
  5210126 L_PHP.3xCpG.bed
  4392384 L_RAP.3xCpG.bed
  4414681 L_RHP.3xCpG.bed
```

#### Sorting the merged files so scaffolds are all in the same order and multiIntersectBed will run correctly. Run for loop using bedtools to sort all .tab files

`nano bedtools.sort.cpg.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/genome_feature
#SBATCH --cpus-per-task=3

module load BEDTools/2.27.1-foss-2018b

for f in *.3xCpG.bed
do
  STEM=$(basename "${f}" .3xCpG.bed)
  bedtools sort -i "${f}" \
  > "${STEM}".3xCpG_sorted.bed
done
```

#### Find CpGs common to all groups in comparison (e.g. Life Stage (A vs L), Adult origin, treatment, transplant, Larval origin, treatment)

`nano cpg_group_comp.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/genome_feature
#SBATCH --cpus-per-task=3

# Adult and larval samples transplanted to the Patch site (Life stage) Comparisons

cat *P.3xCpG_sorted.bed |\
sort | \
uniq -c | \
awk '{if($1==3)print $2"\t"$3"\t"$4}' \
> lifestage.3xCpG.allgrps.bed

# Adult comparison

cat A*.3xCpG_sorted.bed |\
sort | \
uniq -c | \
awk '{if($1==3)print $2"\t"$3"\t"$4}' \
> adult.3xCpG.allgrps.bed

# Larval comparison

cat L*.3xCpG_sorted.bed |\
sort | \
uniq -c | \
awk '{if($1==3)print $2"\t"$3"\t"$4}' \
> larval.3xCpG.allgrps.bed

# Sort the bedfiles
module load BEDTools/2.27.1-foss-2018b

for f in *.3xCpG.allgrps.bed
do
  STEM=$(basename "${f}" .3xCpG.allgrps.bed)
  bedtools sort -i "${f}" \
  > "${STEM}".3xCpG.allgrps_sorted.bed
done

```

`wc -l lifestage.3xCpG.allgrps.bed`
825108 lifestage.3xCpG.allgrps.bed

`wc -l adult.3xCpG.allgrps.bed `
802968 adult.3xCpG.allgrps.bed

`wc -l larval.3xCpG.allgrps.bed`
1398867 larval.3xCpG.allgrps.bed


#### Determine binned features that common CpGs overlap with (bedtools intersect)

`nano bin_feat_cpg.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/genome_feature
#SBATCH --cpus-per-task=3

module load BEDTools/2.27.1-foss-2018b

FILES=*.3xCpG.allgrps_sorted.bed

for f in $FILES
do
intersectBed \
-a ${f} \
-b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides-v1.2Kbin.uniq.bed \
-wb \
> $(basename ${f%.3xCpG.allgrps_sorted.bed})_features.txt
done
```

#### Filter features for those with at least 3 CpGs

`nano filter_feat.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/genome_feature
#SBATCH --cpus-per-task=3

FILES=*_features.txt

for f in $FILES
do
cat ${f} |\
awk -F"\t" '{print $4,$5,$6,$7}' |\
sort |\
uniq -c |\
awk '{if($1>2)print $2"\t"$3"\t"$4"\t"$5}'\
> $(basename ${f%_features.txt})_features.3CpG.txt
done
```

```
[kevin_wong1@n063 genome_feature]$ wc -l adult_features.3CpG.txt
579330 adult_features.3CpG.txt
[kevin_wong1@n063 genome_feature]$ wc -l adult_features.txt
4182696 adult_features.txt

[kevin_wong1@n063 genome_feature]$ wc -l larval_features.3CpG.txt
808223 larval_features.3CpG.txt
[kevin_wong1@n063 genome_feature]$ wc -l larval_features.txt   
7322259 larval_features.txt

[kevin_wong1@n063 genome_feature]$ wc -l lifestage_features.3CpG.txt
592455 lifestage_features.3CpG.txt
[kevin_wong1@n063 genome_feature]$ wc -l lifestage_features.txt   
4293893 lifestage_features.txt
```

scp 'kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/genome_feature/*_features.3CpG.txt' ~/MyProjects/Thermal_Transplant_Molecular/output/WGBS/genome_feature
