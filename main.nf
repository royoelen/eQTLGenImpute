def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     eQTL-Catalogue/genimpute v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile eqtl_catalogue -resume\
        --bfile /gpfs/hpc/projects/genomic_references/CEDAR/genotypes/PLINK_100718_1018/CEDAR\
        --output_name CEDAR_GRCh37_genotyped\
        --outdir CEDAR

    Mandatory arguments:
      --bfile                       Path to the input unimputed plink bgen files (without extensions bed/bim/fam).
      --output_name                 Prefix for the output files.

    CrossMap arguments:
      --target_ref                  Reference genome fasta file for the target genome assembly (e.g. GRCh38).
      --chain_file                  Chain file to translate genomic cooridnates from the source assembly to target assembly.

    Genotype harmonisation & QC:
      --harmonise_genotypes         Run GenotypeHarmonizer on the raw genotypes to correct flipped/swapped alleles (default: true)
      --ref_panel_hg38              Reference panel used by GenotypeHarmonizer. Ideally should match the reference panel used for imputation (hg38).
      --ref_panel_hg19              Reference panel used for strand fixing. Ideally should match the reference your unimputed data is in (hg19).
      --ref_genome                  Reference genome fasta file for the raw genotypes (typically GRCh37).

    dbSNP annotation:
    --dbSNP_hg19                    dbSNP vcf for annotating with rs IDs (before LiftOver, in hg19).
    --dbSNP_hg38                    dbSNP vcf for annotating with rs IDs (after LiftOver, in hg38).

    Phasing & Imputation:
      --eagle_genetic_map           Eagle genetic map file
      --eagle_phasing_reference     Phasing reference panel for Eagle (1000 Genomes Phase 3 high coverage)
      --minimac_imputation_reference Imputation reference panel for Minimac4 in M3VCF format (1000 Genomes Phase 3 high coverage)
    
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    
    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Define input channels
Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome fasta file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

Channel
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

Channel
    .fromPath(params.ref_panel_hg19)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")] }
    .into { ref_panel_harmonise_genotypes_hg19; ref_panel_fixref_hg19 }

Channel
    .fromPath(params.ref_panel_hg38)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")] }
    .into { ref_panel_harmonise_genotypes_hg38; ref_panel_fixref_genotypes_hg38 }

Channel
    .fromPath( "${params.eagle_phasing_reference}*" )
    .ifEmpty { exit 1, "Eagle phasing reference not found: ${params.eagle_phasing_reference}" }
    .set { phasing_ref_ch }

Channel
    .fromPath( "${params.minimac_imputation_reference}*" )
    .ifEmpty { exit 1, "Minimac4 imputation reference not found: ${params.minimac_imputation_reference}" }
    .set { imputation_ref_ch }

Channel
    .fromPath(params.eagle_genetic_map)
    .ifEmpty { exit 1, "Eagle genetic map file not found: ${params.eagle_genetic_map}" } 
    .set { genetic_map_ch }

Channel
    .fromPath(params.chain_file)
    .ifEmpty { exit 1, "CrossMap.py chain file not found: ${params.chain_file}" } 
    .set { chain_file_ch }

Channel
    .fromPath(params.target_ref)
    .ifEmpty { exit 1, "CrossMap.py target reference genome file: ${params.target_ref}" } 
    .into { target_ref_ch; target_ref_ch2 }

Channel
    .fromPath(params.dbSNP_hg19)
    .map { dbSNP_hg19 -> [file("${dbSNP_hg19}.vcf.gz"), file("${dbSNP_hg19}.vcf.gz.tbi")]}
    .ifEmpty { exit 1, "dbSNP hg19 file: ${params.dbSNP_hg19}" }
    .into{dbSNP_hg19_ch; dbSNP_hg19_reffix_ch}

Channel
    .fromPath(params.dbSNP_hg38)
    .map { dbSNP_hg38 -> [file("${dbSNP_hg38}.vcf.gz"), file("${dbSNP_hg38}.vcf.gz.tbi")]}
    .ifEmpty { exit 1, "dbSNP hg38 file: ${params.dbSNP_hg38}" }
    .into{dbSNP_hg38_ch; dbSNP_hg38_fixref_ch}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
eQTL-Catalogue/genimpute v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'eQTL-Catalogue/genimpute'
summary['Pipeline Version']         = workflow.manifest.version
summary['Run Name']                 = custom_runName ?: workflow.runName
summary['PLINK bfile']              = params.bfile
summary['Reference genome']         = params.ref_genome
summary['Harmonise genotypes']      = params.harmonise_genotypes
summary['Harmonisation ref panel hg19']  = params.ref_panel_hg19
summary['Harmonisation ref panel hg38']  = params.ref_panel_hg38
summary['Eagle genetic map']        = params.eagle_genetic_map
summary['Eagle reference panel']    = params.eagle_phasing_reference
summary['Minimac4 reference panel'] = params.minimac_imputation_reference
summary['CrossMap reference genome'] = params.target_ref
summary['CrossMap chain file']      = params.chain_file
summary['dbSNP hg19 file']          = params.dbSNP_hg19
summary['dbSNP hg38 file']          = params.dbSNP_hg38
summary['Max Memory']               = params.max_memory
summary['Max CPUs']                 = params.max_cpus
summary['Max Time']                 = params.max_time
summary['Output name']              = params.output_name
summary['Output dir']               = params.outdir
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']            = params.awsregion
   summary['AWS Queue']             = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

process harmonise_genotypes_hg19{

    cpus 1
    memory '30 GB'
    time '24h'

    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch
    set file(hg19_ref_vcf_gz), file(hg19_ref_vcf_gz_index) from ref_panel_harmonise_genotypes_hg19.collect()

    output:
    set file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam") into harmonised_genotypes_hg19

    script:
    """
    java -jar /usr/bin/GenotypeHarmonizer.jar \
    --input ${study_name_bed.simpleName} \
    --inputType PLINK_BED \
    --ref ${hg19_ref_vcf_gz} \
    --refType VCF \
    --update-id \
    --output harmonised
    """
}

process plink_to_vcf{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from harmonised_genotypes_hg19

    output:
    file "harmonised_hg19.vcf" into harmonized_hg19_vcf_ch

    script:
    """
    plink2 --bfile ${study_name_bed.simpleName} --recode vcf-iid --chr 1-22 --out harmonised_hg19
    """
}

// Use genotypeharmonizer instead
/*
process annotate_with_hg19_rsID{
    input:
    file(vcf) from raw_vcf_ch
    set file(dbSnpRefVcf), file(dbSnpRefIndex) from dbSNP_hg19_ch

    output:
    file "vcf_annotated.vcf.gz" into annotated_vcf_ch

    script:
    """
    # Bgzip and tabix
    bgzip ${vcf}
    echo "Bgzipped!"
    tabix -p vcf ${vcf}.gz

    # Annotate with dbSNP rs IDs
    bcftools annotate \
    -a ${dbSnpRefVcf} \
    -c ID \
    -o vcf_temp.vcf \
    ${vcf}.gz

    echo "rs IDs added!"

    # remove anything before "rs" in the SNP ID column
    awk -F'\t' -vOFS='\t' '{ gsub("^.*rs", "rs", \$3) ; print }' vcf_temp.vcf > vcf_annotated.vcf
    rm vcf_temp.vcf

    bgzip vcf_annotated.vcf
    """
}
*/

// Alternative would be to use unsafe way of fixing reference alleles (https://samtools.github.io/bcftools/howtos/plugin.fixref.html)
process vcf_fixref_hg19{
    input:
    file input_vcf from harmonized_hg19_vcf_ch
    file fasta from ref_genome_ch.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_fixref_hg19

    output:
    file "fixref_hg19.vcf.gz" into crossmap_vcf_input

    script:
    """
    bgzip ${input_vcf}
    bcftools index ${input_vcf}.gz
    #bcftools +fixref ${input_vcf} -Oz -o fixref.vcf.gz -- -d -f ${fasta} -m flip

    bcftools +fixref ${input_vcf}.gz -- -f ${fasta} -i ${vcf_file} | \
    bcftools norm --check-ref x -f ${fasta} -Oz -o fixref_hg19.vcf.gz
    """
}

process crossmap_genotypes{
    input:
    file chain_file from chain_file_ch.collect()
    file target_ref from target_ref_ch.collect()
    file vcf from crossmap_vcf_input

    output:
    set file("${vcf.simpleName}_mapped_sorted.vcf.gz"), file("${vcf.simpleName}_mapped_sorted.vcf.gz.tbi") into lifted_vcf

    shell:
    """
    # Exclude structural variants, because they break latest version of CrossMap.py
    bcftools view --exclude-types other ${vcf} -Oz -o ${vcf.simpleName}_noSVs.vcf.gz
    
    # Run CrossMap.py
    CrossMap.py vcf ${chain_file} ${vcf.simpleName}_noSVs.vcf.gz ${target_ref} ${vcf.simpleName}_mapped.vcf
    
    # This did not work because of the contig order in the header changed after LiftOver
    #bcftools sort ${vcf.simpleName}_mapped.vcf > ${vcf.simpleName}_mapped_sorted.vcf
    grep "^#" ${vcf.simpleName}_mapped.vcf > ${vcf.simpleName}_mapped_sorted.vcf
    grep -v "^#" ${vcf.simpleName}_mapped.vcf | sort -k1,1V -k2,2g >> ${vcf.simpleName}_mapped_sorted.vcf

    bgzip ${vcf.simpleName}_mapped_sorted.vcf
    tabix -p vcf ${vcf.simpleName}_mapped_sorted.vcf.gz

    # TODO: check if the contig order in the header matters for any of the following steps!

    """
}

// Use VCF folder instead
process harmonise_genotypes_hg38{

    cpus 1
    memory '30 GB'
    time '24h'

    input:
    set file(study_name_vcf), file(study_name_tbi) from lifted_vcf
    set file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes_hg38.collect()

    output:
    set file("harmonised.vcf.gz"), file("harmonised.vcf.gz.tbi") into harmonised_genotypes

    script:
    """
    java -jar /usr/bin/GenotypeHarmonizer.jar \
    --input ${study_name_vcf.simpleName} \
    --inputType VCF \
    --ref ${vcf_file} \
    --refType VCF \
    --update-id \
    --output harmonised

    plink2 --bfile harmonised --recode vcf-iid --chr 1-22 --snps-only no-DI --out harmonised
    bgzip harmonised.vcf
    tabix -p vcf harmonised.vcf.gz
    """
}

// This should be removed when harmonization reference is already containing rs IDs!
/*
process annotate_with_hg38_rsID{
    input:
    set file(vcf), file(index) from harmonised_genotypes
    set file(dbSnpRefVcf), file(dbSnpRefIndex) from  dbSNP_hg38_ch

    output:
    set file("harmonized_annotated.vcf.gz"), file("harmonized_annotated.vcf.gz.tbi") into hg38annotated_vcf_ch

    script:
    """
    # Annotate with dbSNP rs IDs
    bcftools annotate \
    -a ${dbSnpRefVcf} \
    -c ID \
    -o vcf_annotated.vcf \
    ${vcf}

    echo "rs IDs added!"

    # Bgzip
    bgzip harmonized_annotated.vcf
    echo "Bgzipped!"
    tabix -p vcf harmonized_annotated.vcf.gz
    """
}
*/
process vcf_fixref_hg38{
    input:
    set file(input_vcf), file(input_vcf_tbi) from harmonised_genotypes
    file fasta from target_ref_ch2.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_fixref_genotypes_hg38

    output:
    file "fixref_hg38.vcf.gz" into fixed_to_filter

    script:
    """
    # Fixing 
    bcftools index ${input_vcf}
    bcftools +fixref ${input_vcf} -- -f ${fasta} -i ${vcf_file} | \
    bcftools norm --check-ref x -f ${fasta} -Oz -o fixref_hg38.vcf.gz
    """
}

process filter_preimpute_vcf{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}_preimpute.vcf.gz" else null }

    input:
    file input_vcf from fixed_to_filter

    output:
    set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into split_vcf_input, missingness_input

    script:
    """
    #Index
    bcftools index ${input_vcf}

    #Add tags
    bcftools +fill-tags ${input_vcf} -Oz -o tagged.vcf.gz

    #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
     bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
     bcftools filter -e "ALT='.'" |\
     bcftools norm -d all |\
     bcftools norm -m+any |\
     bcftools view -m2 -M2 -Oz -o filtered.vcf.gz

     #Index the output file
     bcftools index filtered.vcf.gz
    """
}

process calculate_missingness{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "genotypes.imiss") "${params.output_name}.imiss" else null }
    
    input:
    set file(input_vcf), file(input_vcf_index) from missingness_input 

    output:
    file "genotypes.imiss" into missing_individuals

    script:
    """
    vcftools --gzvcf ${input_vcf} --missing-indv --out genotypes
    """
}

process split_by_chr{
    publishDir "${params.outdir}/preimpute/split_chr", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".vcf.gz") > 0) filename else null }
    
    input:
    tuple file(input_vcf), file(input_vcf_index) from split_vcf_input
    each chr from Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    output:
    tuple val(chr), file("chr_${chr}.vcf.gz") into individual_chromosomes

    script:
    """
    bcftools view -r ${chr} ${input_vcf} -Oz -o chr_${chr}.vcf.gz
    """
}
 
process eagle_prephasing{
    input:
    tuple chromosome, file(vcf) from individual_chromosomes
    file genetic_map from genetic_map_ch.collect()
    file phasing_reference from phasing_ref_ch.collect()

    output:
    tuple chromosome, file("chr_${chromosome}.phased.vcf.gz") into phased_vcf_cf

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=chr${chromosome}.bcf \
    --geneticMapFile=${genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr_${chromosome}.phased \
    --numThreads=8
    """
}

process minimac_imputation{
    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.dose.vcf.gz"
 
    input:
    set chromosome, file(vcf) from phased_vcf_cf
    file imputation_reference from imputation_ref_ch.collect()

    output:
    tuple chromosome, file("chr_${chromosome}.dose.vcf.gz") into imputed_vcf_cf

    script:
    """
    minimac4 --refHaps chr${chromosome}.m3vcf.gz \
    --haps ${vcf} \
    --prefix chr_${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome
    """
}
