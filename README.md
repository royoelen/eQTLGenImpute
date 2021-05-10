# eQTLGen/genimpute workflow
Genotype imputation and quality control workflow used by the eQTLGen phase II. This is modified from the genotype imputation workflow developed by eQTL Catalogue team (https://github.com/eQTL-Catalogue/genimpute).


Performs the following main steps:

**Pre-imputation QC:**
- Convert the genotypes to the VCF format with [PLINK](https://www.cog-genomics.org/plink/1.9/).
- Convert raw genotypes to GRCh38 coordinates with CrossMap.py v0.4.1
- Align raw genotypes to the reference panel with [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer).
- Exclude variants with Hardy-Weinberg p-value < 1e-6, missingness > 0.05 and minor allele frequency < 0.01 with [bcftools](https://samtools.github.io/bcftools/)
- Calculate individual-level missingness using [vcftools](https://vcftools.github.io/perl_module.html).

**Imputation:**
- Genotype pre-phasing with Eagle 2.4.1 
- Genotype imputation with Minimac4

## Contributors

Original pipeline was developed:

* Kaur Alasoo
* Liina Anette Pärtel
* Mark-Erik Kodar

Pipeline was adjusted to work with 1000G p3 30X WGS reference panel:

* Urmo Võsa
