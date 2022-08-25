# eQTLGen genotype imputation pipeline

Genotype imputation and quality control workflow used by the eQTLGen phase II. This is modified from the genotype imputation workflow developed by eQTL Catalogue team (https://github.com/eQTL-Catalogue/genimpute).

Performs the following main steps:

**Pre-imputation preprocessing and QC:**
- Convert raw genotypes to GRCh38 coordinates with CrossMap.py v0.4.1 (http://crossmap.sourceforge.net/).
- Align raw genotypes to the reference panel with [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer).
- Convert the genotypes to the VCF format with [PLINK](https://www.cog-genomics.org/plink/1.9/).
- Fixes alleles to match reference panel with [bcftools +fixref](https://samtools.github.io/bcftools/howtos/plugin.fixref.html).
- Exclude variants with Hardy-Weinberg p-value < 1e-6, missingness > 0.05 and minor allele frequency < 0.01 with [bcftools](https://samtools.github.io/bcftools/)
- Calculate individual-level missingness using [vcftools](https://vcftools.github.io/perl_module.html).

**Imputation:**
- Genotype pre-phasing with Eagle 2.4.1 
- Genotype imputation with Minimac4

## Usage information

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=8 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

`git clone https://github.com/eQTLGen/eQTLGenImpute.git`

Or just download this from the gitlab/github download link and unzip.

### Input files

- Input genotype files

Pipeline expects as an input the folder with unimputed plink `.bed/.bim/.fam` files which are in hg19. Full path to file name without extension.

- Reference files

Pipeline needs several reference files to do data processing, QC, and imputation:

- hg38 .vcf.gz reference for fixing the alleles and harmonizing the data after CrossMap to hg38
- hg38 reference genome .fasta file
- Phasing reference (hg38)
- Genetic map for phasing
- Imputation reference
- CrossMap hg19 --> hg38 chain file (comes with the pipeline)

These are organised to the on folder and all you need to do is to download the tar.gz file, unzip, and change the path in the relevant script template. Download the reference from [here](TODO).

### Running the imputation command

Go to folder `ConvertVcf2Hdf5` and modify the Slurm script template `submit_GenotypeConversion_pipeline_template.sh` with your input paths. Some of the paths are pre-filled, assuming that you follow [eQTLGen phase II cookbook](https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook) and its recommended folder structure, however you can also use custom paths.

```bash
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

cohortname=[name of your cohort]
qc_input_folder=../../1_DataQC/output # folder with QCd genotype and expression data, output of DataQC pipeline.
output_path=../output/ # Output path.

# Command
NXF_VER=21.10.6 ${nextflow_path}/nextflow run eQTLGenImpute.nf \
--qcdata ${qc_input_folder} \
--target_ref ${reference_path}/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--ref_panel_hg38 ${reference_path}/ref_panel_QC/30x-GRCh38_NoSamplesSorted \
--eagle_genetic_map ${reference_path}/phasing/genetic_map/genetic_map_hg38_withX.txt.gz \
--eagle_phasing_reference ${reference_path}/phasing/phasing_reference/ \
--minimac_imputation_reference ${reference_path}/imputation/ \
--cohort_name ${cohortname} \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume
```

You can save the modified script version to informative name, e.g. `submit_imputation_pipeline_[**CohortName_PlatformName**].sh`.

Then submit the job `sbatch submit_imputation_pipeline_[**CohortName_PlatformName**].sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate `work` directory is made to the folder and contains all interim files.

### Monitoring and debugging

- Monitoring:
  - Monitor the `slurm-***.out` log file and check if all the steps finish without error. Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.
  - Use `squeue -u [YourUserName]` to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file  is automatically written to your output folder `pipeline_info` subfolder, for potential errors or warnings. E.g. `output/pipeline_info/DataQcReport.html`.
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from `work` folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - `.command.sh`: script which was submitted
  - `.command.log`: log file for seeing the analysis outputs/errors.
  - `.command.err`: file which lists the errors, if any.
  
### Output

After successful completion of the pipeline, there should be three folders in your `output` path: `preimpute`, `postimpute` and `pipeline_info`. `postimpute` contains imputed genotype data in `.vcf.gz` format, filtered by MAF>0.01. `preimpute` data contains quality-controlled genotype data before imputation: this can be used for secondary purposes (e.g. calculate genotype MDSs) or just deleted. `pipeline_info` contains runlogs from the Nextflow run, useful for debugging.

## Acknowledgements

Original eQTL Catalogue pipeline was developed:

* Kaur Alasoo
* Liina Anette Pärtel
* Mark-Erik Kodar

Original pipeline was adjusted to work with 1000G 30X WGS reference panel:

* Ralf Tambets

Elements of those pipelines were adjusted to work with 1000G 30X WGS reference panel and accustomised for eQTLGen phase II consortium analyses:

* Urmo Võsa

### Citations

[Zhao, H., Sun, Z., Wang, J., Huang, H., Kocher, J. P., &#38; Wang, L. (2014). CrossMap: a versatile tool for coordinate conversion between genome assemblies. Bioinformatics, 30(7), 1006–1007. https://doi.org/10.1093/BIOINFORMATICS/BTT730](https://academic.oup.com/bioinformatics/article/30/7/1006/234947)

[Deelen, P., Bonder, M. J., van der Velde, K. J., Westra, H. J., Winder, E., Hendriksen, D., Franke, L., &#38; Swertz, M. A. (2014). Genotype harmonizer: Automatic strand alignment and format conversion for genotype data integration. BMC Research Notes, 7(1), 901. https://doi.org/10.1186/1756-0500-7-901](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-901)

[Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R. E., Lunter, G., Marth, G. T., Sherry, S. T., McVean, G., &#38; Durbin, R. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156–2158. https://doi.org/10.1093/BIOINFORMATICS/BTR330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296)

[Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., &#38; Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), 1–4. https://doi.org/10.1093/GIGASCIENCE/GIAB008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722)

[Loh, P.-R., Danecek, P., Palamara, P. F., Fuchsberger, C., A Reshef, Y., K Finucane, H., Schoenherr, S., Forer, L., McCarthy, S., Abecasis, G. R., Durbin, R., &#38; L Price, A. (2016). Reference-based phasing using the Haplotype Reference Consortium panel. Nature Genetics, 48(11), 1443–1448. https://doi.org/10.1038/ng.3679](https://www.nature.com/articles/ng.3679)

[Das, S., Forer, L., Schönherr, S., Sidore, C., Locke, A. E., Kwong, A., Vrieze, S. I., Chew, E. Y., Levy, S., McGue, M., Schlessinger, D., Stambolian, D., Loh, P.-R., Iacono, W. G., Swaroop, A., Scott, L. J., Cucca, F., Kronenberg, F., Boehnke, M., … Fuchsberger, C. (2016). Next-generation genotype imputation service and methods. Nature Genetics, 48(10), 1284–1287. https://doi.org/10.1038/ng.3656](https://www.nature.com/articles/ng.3656)


### Contacts

For this Nextflow pipeline: urmo.vosa at gmail.com
