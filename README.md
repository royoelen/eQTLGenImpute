# eQTLGen/genimpute workflow
Genotype imputation and quality control workflow used by the eQTLGen phase II. This is modified from the genotype imputation workflow developed by eQTL Catalogue team (https://github.com/eQTL-Catalogue/genimpute).


Performs the following main steps:

**Pre-imputation QC:**
- Convert raw genotypes to GRCh38 coordinates with CrossMap.py v0.4.1.
- Align raw genotypes to the reference panel with [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer).
- Convert the genotypes to the VCF format with [PLINK](https://www.cog-genomics.org/plink/1.9/).
- Fixes alleles to match reference panel with [bcftools +fixref](https://samtools.github.io/bcftools/howtos/plugin.fixref.html).
- Exclude variants with Hardy-Weinberg p-value < 1e-6, missingness > 0.05 and minor allele frequency < 0.01 with [bcftools](https://samtools.github.io/bcftools/)
- Calculate individual-level missingness using [vcftools](https://vcftools.github.io/perl_module.html).

**Imputation:**
- Genotype pre-phasing with Eagle 2.4.1 
- Genotype imputation with Minimac4

## Usage information

### Input files

Pipeline expects as an input the folder with unimputed plink bed/bim/fam files which are in hg19.

### Reference files

Pipeline needs several reference files to do data processing, QC, and imputation:

- hg38 .vcf.gz reference for fixing the alleles and harmonizing the data after CrossMap to hg38
- hg38 reference genome .fasta file
- Phasing reference (hg38)
- Genetic map for phasing
- Imputation reference
- CrossMap hg19 --> hg38 chain file (comes with the pipeline)

These are organised to the on folder and all you need to do is to download the tar.gz file, unzip and change the path in the relevant script template.

### Running the imputation command

Replace the required paths in the script template.

´´´
#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="ImputeGenotypes"

# Load required modules
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
--target_ref ${reference_path}/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--ref_panel_hg38 ${reference_path}/hg38/ref_panel_QC/30x-GRCh38_NoSamplesSorted \
--eagle_genetic_map ${reference_path}/hg38/phasing/genetic_map/genetic_map_hg38_withX.txt.gz \
--eagle_phasing_reference ${reference_path}/hg38/phasing/phasing_reference/ \
--minimac_imputation_reference ${reference_path}/hg38/imputation/ \
--output_name ${output_name} \
--outdir ${output_path}  \
--profile slurm \
-resume

´´´

## Contributors

Original eQTL Catalogue pipeline was developed:

* Kaur Alasoo
* Liina Anette Pärtel
* Mark-Erik Kodar

Original pipeline was adjusted to work with 1000G p3 30X WGS reference panel:

* Ralf Tambets

Elements of those pipelines were adjusted to work with 1000G 30X WGS reference panel and accustomised for eQTLGen consortium analyses:

* Urmo Võsa
