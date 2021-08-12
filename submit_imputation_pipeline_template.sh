#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=[your e-mail address to send notification]
#SBATCH --job-name="ImputeGenotypes"

module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# Define paths
nextflow_path=[full path to your Nextflow executable]
reference_path=[full path to your folder with reference files]

input_path=[full path to your input genotype folder]
output_name=[name of the output files]
output_path=[name of the output path]

# Command
${nextflow_path}/nextflow run main.nf \
--bfile ${input_path} \
--source_ref ${reference_path}/hg19/ref_genome_QC/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
--target_ref ${reference_path}/hg38/ref_genome_QC/Homo_sapiens_assembly38.fasta \
--ref_panel_hg19 ${reference_path}/hg19/ref_panel_QC/1000G_GRCh37_variant_information \
--ref_panel_hg38 ${reference_path}/hg38/ref_panel_QC/30x-GRCh38_NoSamplesSorted \
--dbSNP_hg19 ${reference_path}/hg19/SNP_annotation/00-All \
--dbSNP_hg38 ${reference_path}/hg38/SNP_annotation/00-All \
--eagle_genetic_map ${reference_path}/hg38/phasing/genetic_map/genetic_map_hg38_withX.txt.gz \
--eagle_phasing_reference ${reference_path}/hg38/phasing/phasing_reference/ \
--minimac_imputation_reference ${reference_path}/hg38/imputation/ \
--output_name ${output_name} \
--outdir ${output_path}  \
-profile eqtlgen \
-resume
