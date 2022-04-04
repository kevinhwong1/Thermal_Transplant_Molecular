# Symbiodiniaceae *S. linucheae* (A4) Thermal Transplant WGBS Workflow

## Methylation Quantification with Methylseq for Symbiodiniaceae *S. linucheae* (A4)

### methylseq_A4.sh

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
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_A4
#SBATCH --exclusive

# load modules needed

module load Nextflow/21.03.0

# run nextflow methylseq

nextflow run nf-core/methylseq \
-profile singularity \
--aligner bismark \
--igenomes_ignore \
--fasta /data/putnamlab/hputnam/Coral_Meth/A4_Past/S.linucheae_CCMP2456.genome.fasta \
--save_reference \
--input '/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir A4_Past_TT \
-name A4_Past_TT
```

```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_A4/A4_Past_TT/MultiQC/multiqc_report.html MyProjects/Thermal_Transplant_Molecular/output/multiqc_report_A4.html
```
