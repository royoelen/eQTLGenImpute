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
module load java/11.0.2
module load singularity/3.5.3
module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

# Define paths and arguments
nextflow_path=../../tools # folder where Nextflow executable is.
reference_path=../hg38 # folder where you unpacked the reference files.

cohort_name=[name of your cohort]
qc_input_folder=../../1_DataQC/output/ # folder with QCd genotype and expression data, output of DataQC pipeline.
output_path=../output/ # Output path.
genome_build="GRCh37" # Adjust if your genotype data is in different genotype build. Options: hg18 or GRCh36, hg19 or GRCh37, hg38 or GRCh38 

# Command
NXF_VER=21.10.6 ${nextflow_path}/nextflow run eQTLGenImpute.nf \
--qcdata ${qc_input_folder} \
--target_ref ${reference_path}/genome_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--ref_panel_hg38 ${reference_path}/harmonizing_reference/30x-GRCh38_NoSamplesSorted \
--eagle_genetic_map ${reference_path}/phasing_reference/genetic_map/genetic_map_hg38_withX.txt.gz \
--eagle_phasing_reference ${reference_path}/phasing_reference/phasing/ \
--minimac_imputation_reference ${reference_path}/imputation_reference/ \
--cohort_name ${cohort_name} \
--genome_build ${genome_build} \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume
