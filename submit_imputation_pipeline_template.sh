#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="ImputeGenotypes"

# These are needed modules in UT HPC to get Singularity and Nextflow running.
# Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# Define paths and arguments
nextflow_path=../../tools # folder where Nextflow executable is.
reference_path=../hg38 # folder where you unpacked the reference files.

cohort_name=[name of your cohort]
qc_input_folder=../../1_DataQC/output/postimpute/ # folder with QCd genotype and expression data, output of DataQC pipeline.
output_path=../output/ # Output path.

# Command
NXF_VER=21.10.6 ${nextflow_path}/nextflow run eQTLGenImpute.nf \
--qcdata ${qc_input_folder} \
--target_ref ${reference_path}/genome_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--ref_panel_hg38 ${reference_path}/harmonizing_reference/30x-GRCh38_NoSamplesSorted \
--eagle_genetic_map ${reference_path}/phasing_reference/genetic_map/genetic_map_hg38_withX.txt.gz \
--eagle_phasing_reference ${reference_path}/phasing_reference/phasing/ \
--minimac_imputation_reference ${reference_path}/imputation_reference/ \
--cohort_name ${cohort_name} \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume