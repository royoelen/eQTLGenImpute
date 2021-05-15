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
      --bfile                       Path to the input unimputed plink bgen files (without extensions bed/bim/fam) 
      --output_name                 Prefix for the output files
      --chromosome_names            File with chromosome names for hg19 --> GRCh38.

    CrossMap arguments:
      --target_ref                  Reference genome fasta file for the target genome assembly (e.g. GRCh38).
      --chain_file                  Chain file to translate genomic cooridnates from the source assembly to target assembly.

    Genotype harmonisation & QC:
      --harmonise_genotypes         Run GenotypeHarmonizer on the raw genotypes to correct flipped/swapped alleles (default: true)
      --ref_panel_hg38              Reference panel used by GenotypeHarmonizer. Ideally should match the reference panel used for imputation (hg38).
      --ref_panel_hg19              Reference panel used for strand fixing. Ideally should match the reference your unimputed data is in (hg19).
      --ref_genome                  Reference genome fasta file for the raw genotypes (typically GRCh37).

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
    .from(params.ref_panel_hg38)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .into { ref_panel_harmonise_genotypes; ref_panel_fix_strands_hg38 }

Channel
    .from(params.ref_panel_hg19)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .into { ref_panel_fix_strands } 

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
    .fromPath(params.chromosome_names)
    .ifEmpty { exit 1, "Chromosome names file not found: ${params.chromosome_names}" } 
    .set { chr_names_file_ch }


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
summary['Harmonisation ref panel']  = params.ref_panel_hg38
summary['Strand fixing ref panel']  = params.ref_panel_hg19
summary['Eagle genetic map']        = params.eagle_genetic_map
summary['Eagle reference panel']    = params.eagle_phasing_reference
summary['Minimac4 reference panel'] = params.minimac_imputation_reference
summary['CrossMap reference genome'] = params.target_ref
summary['CrossMap chain file']      = params.chain_file
summary['Chromosomes names file']   = params.chromosome_names
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

process plink_to_vcf{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch

    output:
    file "raw.vcf" into raw_vcf_ch

    script:
    """
    # TODO: check if snps-only just-actg option is acceptable
    plink2 --bfile ${study_name_bed.simpleName} --recode vcf-iid --chr 1-22 --snps-only no-DI --out raw
    """
}

/* To consider, if the next step needs rs IDs, herin this step we could add rs IDs from hg19 reference panel 
process annotate_with_rsID{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch

    output:
    file "raw.vcf" into raw_vcf_ch

    script:
    """
    # Annotate with dbSNP rs IDs
    bcftools annotate \
    -a ${dbSnpRefVcf} \
    -c ID \
    -o chr\${chr}_FixedSnpNamesFiltered0005_rsIDs.vcf \
    ${FilteredVcfToAnnotation}

    echo "rs IDs added!"

    # Bgzip and tabix
    bgzip chr\${chr}_FixedSnpNamesFiltered0005_rsIDs.vcf
    echo "Bgzipped!"
    tabix -p vcf chr\${chr}_FixedSnpNamesFiltered0005_rsIDs.vcf.gz
    echo "Tabixed!"
    """
}

*/

// TODO, this is probably unsafe way of fixing reference alleles in hg19 (https://samtools.github.io/bcftools/howtos/plugin.fixref.html), later use Grch37 reference panel instead.
process vcf_fixref_hg19{
    input:
    file input_vcf from raw_vcf_ch
    file fasta from ref_genome_ch.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_fix_strands.collect()

    output:
    file "fixref.vcf.gz" into crossmap_vcf_input

    script:
    """
    bgzip ${input_vcf}
    bcftools index ${input_vcf}.gz
    bcftools +fixref ${input_vcf}.gz -Oz -o fixref.vcf.gz -- -d -f ${fasta} -m flip
    """
}

process crossmap_genotypes{
    input:
    file chain_file from chain_file_ch.collect()
    file target_ref from target_ref_ch.collect()
    file vcf from crossmap_vcf_input

    output:
    file("${vcf.simpleName}_mapped_sorted.vcf") into lifted_vcf

    shell:
    """
    #Exclude structural variants, beause they break latest version of CrossMap.py
    bcftools view --exclude-types other ${vcf} -Oz -o ${vcf.simpleName}_noSVs.vcf.gz
    
    #Run CrossMap.py
    CrossMap.py vcf ${chain_file} ${vcf.simpleName}_noSVs.vcf.gz ${target_ref} ${vcf.simpleName}_mapped.vcf
    
    # This did not work because of the contig order in the header changed after LiftOver
    #bcftools sort ${vcf.simpleName}_mapped.vcf > ${vcf.simpleName}_mapped_sorted.vcf
    grep "^#" ${vcf.simpleName}_mapped.vcf > ${vcf.simpleName}_mapped_sorted.vcf
    grep -v "^#" ${vcf.simpleName}_mapped.vcf | sort -k1,1V -k2,2g >> ${vcf.simpleName}_mapped_sorted.vcf

    # TODO: check if the contig order in the header matters for any of the following steps!

    """
}

process fix_chr_format{
    input:
    file lifted_vcf from lifted_vcf

    output:
    set file("raw_with_chr.vcf.gz"), file("raw_with_chr.vcf.gz.tbi") into mapped_vcf_ch

    shell:
    """   
    # Add "chr" to each chromosome name
    #awk '{ if(\$0 !~ /^#/ && \$0 !~ /^chr/) \
    #print "chr"\$0; else if(match(\$0,/(##contig=<ID=)(.*)/,m)) \
    #print m[1]"chr"m[2]; else print \$0 }'  ${lifted_vcf} > raw_with_chr.vcf

    awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' ${lifted_vcf} > raw_with_chr.vcf

    bgzip raw_with_chr.vcf
    #Index the output file
    tabix -p vcf raw_with_chr.vcf.gz
    """
}

process harmonise_genotypes{
    input:
    set file(study_name_vcf), file(study_name_tbi) from mapped_vcf_ch
    set file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes.collect()

    output:
    set file("harmonised.vcf.gz"), file("harmonised.vcf.gz.tbi") into harmonised_genotypes

    script:
    """
    java -jar /usr/bin/GenotypeHarmonizer.jar \
    --input ${study_name_vcf.simpleName} \
    --inputType VCF \
    --ref ${vcf_file.simpleName} \
    --refType VCF \
    --update-id \
    --output harmonised

    plink2 --bfile harmonised --recode vcf-iid --chr 1-22 --snps-only no-DI --out harmonised
    bgzip harmonised.vcf
    tabix -p vcf harmonised.vcf.gz
    """
}

process vcf_fixref_hg38{
    input:
    set file(input_vcf), file(input_vcf_tbi) from harmonised_genotypes
    file fasta from target_ref_ch2.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_fix_strands_hg38.collect()

    output:
    file "fixref_hg38.vcf.gz" into fixed_to_filter

    script:
    """
    gunzip -c ${input_vcf} > harmonised.vcf
    #cat harmonised.vcf | awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' | sed 's/contig=<ID=/contig=<ID=chr/g' > harmonised_with_chr.vcf
    awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' harmonised.vcf > harmonised_with_chr.vcf
    bgzip harmonised_with_chr.vcf

    # Add "chr" to fasta files
    sed 's/>/>chr/g' ${fasta} > fixed_fasta.fa

    # Fixing 
    bcftools index harmonised_with_chr.vcf.gz
    bcftools +fixref harmonised_with_chr.vcf.gz -- -f fixed_fasta.fa -i ${vcf_file} | \
    bcftools norm --check-ref x -f fixed_fasta.fa -Oz -o fixref_hg38.vcf.gz
    """
}

process filter_preimpute_vcf{
    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}_preimpute.vcf.gz" else null }

    input:
    set file(input_vcf), file(input_vcf_tbi) from fixed_to_filter

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
    bcftools view -r chr${chr} ${input_vcf} -Oz -o chr_${chr}.vcf.gz
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
    --vcfRef=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chromosome}.filtered.shapeit2-duohmm-phased.bcf \
    --geneticMapFile=${genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr_${chromosome}.phased \
    --numThreads=8
    """
}

adjust_chr_format = phased_vcf_cf.combine(chr_names_file_ch)  

process adjust_chr_format{
    input:
    set chromosome, file(vcf), file(chromosome_names) from adjust_chr_format

    output:
    tuple chromosome, file("chr_${chromosome}_fixed.phased.vcf.gz") into phased_vcf_fixed

    script:
    """
    bcftools annotate --rename-chrs ${chromosome_names} ${vcf} -Oz -o chr_${chromosome}_fixed.phased.vcf.gz
    """
}   

process minimac_imputation{
    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.dose.vcf.gz"
 
    input:
    set chromosome, file(vcf) from phased_vcf_fixed
    file imputation_reference from imputation_ref_ch.collect()

    output:
    tuple chromosome, file("chr_${chromosome}.dose.vcf.gz") into imputed_vcf_cf

    script:
    """
    minimac4 --refHaps CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chromosome}.filtered.shapeit2-duohmm-phased.m3vcf.gz \
    --haps ${vcf} \
    --prefix chr_${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome
    """
}

