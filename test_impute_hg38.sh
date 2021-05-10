#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="ImputeGenotypes"

module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=/gpfs/space/GI/GV/Projects/eQTLGenPhase2/tools/

# Using the hg38 reference fasta from here: https://www.biostars.org/p/492887/

${nextflow_path}/nextflow run main.nf \
--bfile /gpfs/space/GI/GV/Projects/eQTLGenPhase2/imputation/data/EGCUT_GSA_HT12v3/filtered/GSA_HT12v3_eQTL_samples \
--harmonise_genotypes true \
--ref_panel_hg38 /gpfs/space/GI/GV/Projects/eQTLGenPhase2/imputation/data/30x-GRCh38_no_samples/30x-GRCh38_NoSamples \
--ref_panel_hg19 /gpfs/hpc/projects/genomic_references/1000G/GRCh37/1000G_GRCh37_variant_information \
--target_ref2 /gpfs/space/GI/GV/Projects/eQTLGenPhase2/imputation/data/fasta/Homo_sapiens_assembly38.fasta \
--eagle_genetic_map /gpfs/space/GI/GV/Projects/eQTLGenPhase2/imputation/data/Eagle_genetic_map/genetic_map_hg38_withX.txt.gz \
--eagle_phasing_reference /gpfs/space/GI/GV/Projects/eQTLGenPhase2/imputation/data/1000Gp3v5_WGS_bcf/ \
--output_name EGCUT_GSA_HT12v3_1000Gp3HighCov_imputed \
--outdir /gpfs/space/GI/GV/Projects/eQTLGenPhase2/imputation/data/EGCUT_GSA_HT12v3_hg38_imputed \
-profile eqtlgen \
-resume