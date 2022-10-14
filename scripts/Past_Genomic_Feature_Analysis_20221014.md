## Genomic Feature Analysis

Following [this pipeline](https://github.com/hputnam/Geoduck_Meth/blob/master/code/Geoduck_Meth_Pipeline.md) created by H. Putnam and [this pipeline](https://github.com/hputnam/Meth_Compare/blob/master/code/03.01-Generating-Genome-Feature-Tracks.ipynb) by Yaamini R. Venkataraman

### 1. Prepare Reference File

#### Creating a reference fasta file

```
module load SAMtools/1.9-foss-2018b

cd /data/putnamlab/kevin_wong1/Past_Genome/
samtools faidx past_filtered_assembly.fasta
```

This should create a fasta index file called that has the chromosome ID and length

`past_filtered_assembly.fasta.fai`

`mkdir 20221014_GF_Analysis`

Double checking lengths: obtain sequence lengths for each "chromosome"
```
awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' \
../../past_filtered_assembly.fasta \
> Past.genome_assembly-sequence-lengths.txt
```

Checking how many columns are in this file, seems like there are extra variables

`cat Past.genome_assembly-sequence-lengths.txt  | awk '{ print NF}'`

There are 6 columns, I only need the first and last (chromosome ID and sequence length)

`cat Past.genome_assembly-sequence-lengths.txt | awk '{ print $1, $6 }' > Past.genome_assembly-sequence-lengths_clean.txt`

Removing the first line as it is currently blank

`sed -i '1d' Past.genome_assembly-sequence-lengths_clean.txt`

`head Past.genome_assembly-sequence-lengths_clean.txt`

```
000000F 3369715
000001F 3221250
000002F 2231938
000003F 2177049
000004F 2151454
000005F 2047925
000006F 1978571
000007F 1916473
000008F 1871416
```

cat Past.genome_assembly-sequence-lengths_clean.txt  | awk '{ print NF}'

#### Adding introns to the gff

```
module load AGAT/0.8.1-foss-2020b
agat_sp_add_introns.pl --gff Pastreoides_noseq_v1.gff --out Pastreoides_noseq_v1_introns.gff
```

### 2. Filtering gff for specific features

#### Gene

`awk '{if ($3 == "gene") print $0;}' Pastreoides_noseq_v1_introns.gff > Past.GFFannotation.gene.gff`

`head Past.GFFannotation.gene.gff`

```
002995F	maker	gene	11	3258	.	-	.	ID=Pastreoides64635;Alias=maker-002995F-snap-gene-0.4;Name=Pastreoides64635
002995F	maker	gene	3264	4210	.	-	.	ID=Pastreoides64636;Alias=maker-002995F-snap-gene-0.5;Name=Pastreoides64636
001331F	maker	gene	182	10442	.	-	.	ID=Pastreoides54855;Alias=maker-001331F-augustus-gene-0.20;Name=Pastreoides54855
001331F	maker	gene	10546	15284	.	-	.	ID=Pastreoides54856;Alias=maker-001331F-augustus-gene-0.21;Name=Pastreoides54856
001331F	maker	gene	18161	27581	.	+	.	ID=Pastreoides54857;Alias=maker-001331F-augustus-gene-0.15;Name=Pastreoides54857
001331F	maker	gene	28749	33669	.	-	.	ID=Pastreoides54858;Alias=augustus_masked-001331F-processed-gene-0.32;Name=Pastreoides54858
001331F	maker	gene	46764	47017	.	+	.	ID=Pastreoides54859;Alias=maker-001331F-snap-gene-0.2;Name=Pastreoides54859
001331F	maker	gene	48363	49860	.	-	.	ID=Pastreoides54860;Alias=snap_masked-001331F-processed-gene-0.50;Name=Pastreoides54860
001331F	maker	gene	68142	74033	.	-	.	ID=Pastreoides54861;Alias=snap_masked-001331F-processed-gene-0.51;Name=Pastreoides54861
001331F	maker	gene	77242	78905	.	-	.	ID=Pastreoides54862;Alias=snap_masked-001331F-processed-gene-0.52;Name=Pastreoides54862
```

`wc -l Past.GFFannotation.gene.gff`

64636 Past.GFFannotation.gene.gff

#### CDS

`awk '{if ($3 == "CDS") print $0;}' Pastreoides_noseq_v1_introns.gff > Past.GFFannotation.cds.gff`

`head Past.GFFannotation.cds.gff`

```
002995F	maker	CDS	11	415	.	-	0	ID=Pastreoides64635-RA:cds;Parent=Pastreoides64635-RA
002995F	maker	CDS	1168	1411	.	-	1	ID=Pastreoides64635-RA:cds;Parent=Pastreoides64635-RA
002995F	maker	CDS	2486	2715	.	-	0	ID=Pastreoides64635-RA:cds;Parent=Pastreoides64635-RA
002995F	maker	CDS	2733	3019	.	-	2	ID=Pastreoides64635-RA:cds;Parent=Pastreoides64635-RA
002995F	maker	CDS	3054	3258	.	-	0	ID=Pastreoides64635-RA:cds;Parent=Pastreoides64635-RA
002995F	maker	CDS	3264	3667	.	-	2	ID=Pastreoides64636-RA:cds;Parent=Pastreoides64636-RA
002995F	maker	CDS	3676	3991	.	-	0	ID=Pastreoides64636-RA:cds;Parent=Pastreoides64636-RA
002995F	maker	CDS	4154	4210	.	-	0	ID=Pastreoides64636-RA:cds;Parent=Pastreoides64636-RA
001331F	maker	CDS	182	274	.	-	0	ID=Pastreoides54855-RA:cds;Parent=Pastreoides54855-RA
001331F	maker	CDS	740	808	.	-	0	ID=Pastreoides54855-RA:cds;Parent=Pastreoides54855-RA
```

`wc -l Past.GFFannotation.cds.gff`

315399 Past.GFFannotation.cds.gff

#### Intron

`awk '{if ($3 == "intron") print $0;}' Pastreoides_noseq_v1_introns.gff > Past.GFFannotation.intron.gff`

`head Past.GFFannotation.intron.gff `

```
002995F	maker	intron	416	1167	.	-	.	ID=intron_added-252663;Parent=Pastreoides64635-RA
002995F	maker	intron	1412	2485	.	-	.	ID=intron_added-252664;Parent=Pastreoides64635-RA
002995F	maker	intron	2716	2732	.	-	.	ID=intron_added-252665;Parent=Pastreoides64635-RA
002995F	maker	intron	3020	3053	.	-	.	ID=intron_added-252666;Parent=Pastreoides64635-RA
002995F	maker	intron	3668	3675	.	-	.	ID=intron_added-252667;Parent=Pastreoides64636-RA
002995F	maker	intron	3992	4153	.	-	.	ID=intron_added-252668;Parent=Pastreoides64636-RA
001331F	maker	intron	275	739	.	-	.	ID=intron_added-219834;Parent=Pastreoides54855-RA
001331F	maker	intron	809	1489	.	-	.	ID=intron_added-219835;Parent=Pastreoides54855-RA
001331F	maker	intron	1589	2463	.	-	.	ID=intron_added-219836;Parent=Pastreoides54855-RA
001331F	maker	intron	2584	3915	.	-	.	ID=intron_added-219837;Parent=Pastreoides54855-RA
```

`wc -l Past.GFFannotation.intron.gff`

252668 Past.GFFannotation.intron.gff


#### Flanking regions

##### All Regions

Create 1kb flanking regions. Subtract existing genes so flanks do not have any overlap.

`module load BEDTools/2.27.1-foss-2018b`


Sorting  gene file:

`sort -i Past.GFFannotation.gene.gff > Past.GFFannotation.gene_sorted.gff`

`head Past.GFFannotation.gene_sorted.gff`

```
000000F	maker	gene	1003060	1006443	.	+	.	ID=Pastreoides00083;Alias=augustus_masked-000000F-processed-gene-3.69;Name=Pastreoides00083
000000F	maker	gene	1009896	1010181	.	+	.	ID=Pastreoides00084;Alias=maker-000000F-snap-gene-3.4;Name=Pastreoides00084
000000F	maker	gene	1010221	1021971	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1025217	1033962	.	+	.	ID=Pastreoides00086;Alias=snap_masked-000000F-processed-gene-3.100;Name=Pastreoides00086
000000F	maker	gene	1035144	1042327	.	-	.	ID=Pastreoides00087;Alias=maker-000000F-snap-gene-3.22;Name=Pastreoides00087
000000F	maker	gene	1037615	1039890	.	+	.	ID=Pastreoides00088;Alias=maker-000000F-snap-gene-3.7;Name=Pastreoides00088
000000F	maker	gene	1046739	1081612	.	-	.	ID=Pastreoides00089;Alias=maker-000000F-snap-gene-3.26;Name=Pastreoides00089
000000F	maker	gene	1085981	1087852	.	+	.	ID=Pastreoides00090;Alias=maker-000000F-augustus-gene-3.44;Name=Pastreoides00090
000000F	maker	gene	1099357	1100275	.	-	.	ID=Pastreoides00091;Alias=maker-000000F-snap-gene-3.32;Name=Pastreoides00091
000000F	maker	gene	1100564	1101556	.	-	.	ID=Pastreoides00092;Alias=maker-000000F-snap-gene-3.31;Name=Pastreoides00092
```

```
flankBed \
-i Past.GFFannotation.gene_sorted.gff \
-g ../../past_filtered_assembly.fasta.fai \
-b 1000 \
| subtractBed \
-a - \
-b Past.GFFannotation.gene_sorted.gff \
> Past.GFFannotation.flanks.gff
```

`head Past.GFFannotation.flanks.gff`

```
000000F	maker	gene	1002060	1003059	.	+	.	ID=Pastreoides00083;Alias=augustus_masked-000000F-processed-gene-3.69;Name=Pastreoides00083
000000F	maker	gene	1006444	1007443	.	+	.	ID=Pastreoides00083;Alias=augustus_masked-000000F-processed-gene-3.69;Name=Pastreoides00083
000000F	maker	gene	1008896	1009895	.	+	.	ID=Pastreoides00084;Alias=maker-000000F-snap-gene-3.4;Name=Pastreoides00084
000000F	maker	gene	1010182	1010220	.	+	.	ID=Pastreoides00084;Alias=maker-000000F-snap-gene-3.4;Name=Pastreoides00084
000000F	maker	gene	1009221	1009895	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1010182	1010220	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1021972	1022971	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1024217	1025216	.	+	.	ID=Pastreoides00086;Alias=snap_masked-000000F-processed-gene-3.100;Name=Pastreoides00086
000000F	maker	gene	1033963	1034962	.	+	.	ID=Pastreoides00086;Alias=snap_masked-000000F-processed-gene-3.100;Name=Pastreoides00086
000000F	maker	gene	1034144	1035143	.	-	.	ID=Pastreoides00087;Alias=maker-000000F-snap-gene-3.22;Name=Pastreoides00087
```

`wc -l Past.GFFannotation.flanks.gff`

133449 Past.GFFannotation.flanks.gff

##### Upstream flanks

* Create flanking regions
* Create upstream flanks (-l) based on strand (-s)
* Subtract existing genes so flanks do not have any overlap

```
flankBed \
-i Past.GFFannotation.gene_sorted.gff \
-g ../../past_filtered_assembly.fasta.fai \
-l 1000 \
-r 0 \
-s \
| subtractBed \
-a - \
-b Past.GFFannotation.gene_sorted.gff \
> Past.GFFannotation.flanks.Upstream.gff
```

`head Past.GFFannotation.flanks.Upstream.gff`

```
000000F	maker	gene	1002060	1003059	.	+	.	ID=Pastreoides00083;Alias=augustus_masked-000000F-processed-gene-3.69;Name=Pastreoides00083
000000F	maker	gene	1008896	1009895	.	+	.	ID=Pastreoides00084;Alias=maker-000000F-snap-gene-3.4;Name=Pastreoides00084
000000F	maker	gene	1009221	1009895	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1010182	1010220	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1024217	1025216	.	+	.	ID=Pastreoides00086;Alias=snap_masked-000000F-processed-gene-3.100;Name=Pastreoides00086
000000F	maker	gene	1042328	1043327	.	-	.	ID=Pastreoides00087;Alias=maker-000000F-snap-gene-3.22;Name=Pastreoides00087
000000F	maker	gene	1081613	1082612	.	-	.	ID=Pastreoides00089;Alias=maker-000000F-snap-gene-3.26;Name=Pastreoides00089
000000F	maker	gene	1084981	1085980	.	+	.	ID=Pastreoides00090;Alias=maker-000000F-augustus-gene-3.44;Name=Pastreoides00090
000000F	maker	gene	1100276	1100563	.	-	.	ID=Pastreoides00091;Alias=maker-000000F-snap-gene-3.32;Name=Pastreoides00091
000000F	maker	gene	1101557	1101623	.	-	.	ID=Pastreoides00092;Alias=maker-000000F-snap-gene-3.31;Name=Pastreoides00092
```

`wc -l Past.GFFannotation.flanks.Upstream.gff`

66757 Past.GFFannotation.flanks.Upstream.gff

##### Downstream flanks

* Create flanking regions
* Create downstream flanks (-r) based on strand (-s)
* Subtract existing genes so flanks do not have any overlap

```
flankBed \
-i Past.GFFannotation.gene_sorted.gff \
-g ../../past_filtered_assembly.fasta.fai \
-l 0 \
-r 1000 \
-s \
| subtractBed \
-a - \
-b Past.GFFannotation.gene_sorted.gff \
> Past.GFFannotation.flanks.Downstream.gff
```

`head Past.GFFannotation.flanks.Downstream.gff`
```
000000F	maker	gene	1006444	1007443	.	+	.	ID=Pastreoides00083;Alias=augustus_masked-000000F-processed-gene-3.69;Name=Pastreoides00083
000000F	maker	gene	1010182	1010220	.	+	.	ID=Pastreoides00084;Alias=maker-000000F-snap-gene-3.4;Name=Pastreoides00084
000000F	maker	gene	1021972	1022971	.	+	.	ID=Pastreoides00085;Alias=snap_masked-000000F-processed-gene-3.99;Name=Pastreoides00085
000000F	maker	gene	1033963	1034962	.	+	.	ID=Pastreoides00086;Alias=snap_masked-000000F-processed-gene-3.100;Name=Pastreoides00086
000000F	maker	gene	1034144	1035143	.	-	.	ID=Pastreoides00087;Alias=maker-000000F-snap-gene-3.22;Name=Pastreoides00087
000000F	maker	gene	1045739	1046738	.	-	.	ID=Pastreoides00089;Alias=maker-000000F-snap-gene-3.26;Name=Pastreoides00089
000000F	maker	gene	1087853	1088852	.	+	.	ID=Pastreoides00090;Alias=maker-000000F-augustus-gene-3.44;Name=Pastreoides00090
000000F	maker	gene	1098357	1099356	.	-	.	ID=Pastreoides00091;Alias=maker-000000F-snap-gene-3.32;Name=Pastreoides00091
000000F	maker	gene	1100276	1100563	.	-	.	ID=Pastreoides00092;Alias=maker-000000F-snap-gene-3.31;Name=Pastreoides00092
000000F	maker	gene	1101557	1101623	.	-	.	ID=Pastreoides00093;Alias=maker-000000F-snap-gene-3.29;Name=Pastreoides00093
```

`wc -l Past.GFFannotation.flanks.Downstream.gff`

66696 Past.GFFannotation.flanks.Downstream.gff

#### Intergenic regions

* Find intergenic regions
* By definition, these are regions that are not in genes and do not include flanking regions
* Subtract defined flanks from intergenic regions


I first had to re-sort my already sorted "Past.GFFannotation.gene_sorted.gff" file. I have no idea why...
```
sort -i Past.GFFannotation.gene_sorted.gff > Past.GFFannotation.gene_sorted.gff
```
```
complementBed \
-i Past.GFFannotation.gene_sorted.gff \
-g ../../past_filtered_assembly_sorted.fasta.fai \
| subtractBed \
-a - \
-b Past.GFFannotation.flanks.gff \
| awk '{print $1"\t"$2"\t"$3}' \
> Past.GFFannotation.intergenic.bed
```

`head Past.GFFannotation.intergenic.bed`

```
000000F	0	5840
000000F	6840	7839
000000F	8841	16231
000000F	17231	23294
000000F	24294	26964
000000F	28766	29075
000000F	29932	31456
000000F	33346	45947
000000F	46947	47132
000000F	48132	61496
```

`wc -l Past.GFFannotation.intergenic.bed`

101146 Past.GFFannotation.intergenic.bed

### 3. Create CpG motifs

I am following [this script](https://robertslab.github.io/sams-notebook/2019/08/21/Data-Wrangling-Create-a-CpG-GFF-from-Pgenerosa_v074-using-EMBOSS-fuzznuc-on-Swoose.html) to characterize the CpGs in the Porites astreoides genome.

`module load EMBOSS/6.6.0-foss-2018b`

```
fuzznuc \
-sequence ../../past_filtered_assembly.fasta \
-pattern CG \
-outfile Past_CpG.gff \
-rformat gff
```

`head Past_CpG.gff`

```
##gff-version 3
##sequence-region 000000F 1 3369715
#!Date 2022-10-14
#!Type DNA
#!Source-version EMBOSS 6.6.0.0
000000F	fuzznuc	nucleotide_motif	1	2	2	+	.	ID=000000F.1;note=*pat pattern:CG
000000F	fuzznuc	nucleotide_motif	9	10	2	+	.	ID=000000F.2;note=*pat pattern:CG
000000F	fuzznuc	nucleotide_motif	22	23	2	+	.	ID=000000F.3;note=*pat pattern:CG
000000F	fuzznuc	nucleotide_motif	51	52	2	+	.	ID=000000F.4;note=*pat pattern:CG
000000F	fuzznuc	nucleotide_motif	56	57	2	+	.	ID=000000F.5;note=*pat pattern:CG
```

### 4. Characterizing CpG methylation (5x data)

* Using the 5x coverage data on the reduced dataset (i.e. removal of samples with low reads) that contains adult-larval pairs

Steps:

1. Characterize overlap between CG motifs and genome feature tracks
2. Download coverage files
3. Characterize methylation for each CpG dinucleotide
4. Characterize genomic locations of all sequenced data, methylated CpGs, sparsely methylated CpGs, and unmethylated CpGs for each sequencing type

#### Set up

`pwd`

/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS

`mkdir genome_feature_20221014`

`cd genome_feature_20221014`
