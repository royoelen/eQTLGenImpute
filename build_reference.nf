//molecular trait data input data
vcf_file_ch = Channel.fromPath(params.vcf_list)
    .ifEmpty { error "Cannot find vcf list file in: ${params.vcf_list}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.chromosome, file(row.path) ]}

Channel
    .fromPath(params.chromosome_names)
    .ifEmpty { exit 1, "Chromosome names file not found: ${params.chromosome_names}" } 
    .set { chr_names_file_ch }

rename_chr_input = vcf_file_ch.combine(chr_names_file_ch)   

process rename_chromosomes{
    container = 'quay.io/biocontainers/bcftools:1.12--h45bccc9_1'

    input:
    tuple val(chr), file(vcf), file(chromosome_names) from rename_chr_input

    output:
    tuple val(chr), file("${chr}.renamed.vcf.gz") into renamed_vcf_ch

    script:
    """
    bcftools annotate --rename-chrs ${chromosome_names} ${vcf} -Oz -o ${chr}.renamed.vcf.gz
    """
}

process create_m3vcf{
    container = "quay.io/eqtlcatalogue/minimac3:v2.0.1"

    publishDir "${params.outdir}/", mode: 'copy', pattern: "*.m3vcf.gz"

    input:
    tuple val(chr), file(vcf) from renamed_vcf_ch

    output:
    tuple val(chr), file("${chr}.m3vcf.gz") into m3vcf_ch

    script:
    """
    minimac3 --refHaps ${vcf} --processReference --prefix ${chr}
    """
}