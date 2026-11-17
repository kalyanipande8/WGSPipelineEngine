// Authored-but-not-executed reference: the real Nextflow pipeline this
// project would run in production (FASTQ -> aligned BAM -> called VCF
// -> annotated VCF), containerized per-stage and deployed on AWS
// Batch. No Nextflow runtime or AWS account is available in this
// environment, so src/pipeline_stages.py implements the identical
// align -> call -> annotate stage sequence in Python as the executed
// mirror, and all measured results in this repo come from that
// mirror, not from this script.

nextflow.enable.dsl = 2

params.samples = "samples.csv"
params.reference = "GRCh38.fa"
params.outdir = "results"

process ALIGN {
    container 'biocontainers/bwa:v0.7.17'
    cpus 4
    input:
        tuple val(sample_id), path(fastq1), path(fastq2)
    output:
        tuple val(sample_id), path("${sample_id}.bam")
    script:
        "bwa mem -t ${task.cpus} ${params.reference} ${fastq1} ${fastq2} | samtools sort -o ${sample_id}.bam"
}

process CALL_VARIANTS {
    container 'broadinstitute/gatk:4.5.0.0'
    cpus 2
    input:
        tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path("${sample_id}.vcf")
    script:
        "gatk HaplotypeCaller -R ${params.reference} -I ${bam} -O ${sample_id}.vcf"
}

process ANNOTATE {
    container 'ensemblorg/ensembl-vep:release_110'
    input:
        tuple val(sample_id), path(vcf)
    output:
        path("${sample_id}.annotated.vcf")
    script:
        "vep -i ${vcf} -o ${sample_id}.annotated.vcf --vcf"
}

workflow {
    samples_ch = Channel.fromPath(params.samples).splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2)) }

    aligned_ch = ALIGN(samples_ch)
    called_ch = CALL_VARIANTS(aligned_ch)
    ANNOTATE(called_ch)
}

// Real concordance evaluation against a gold-standard call set (e.g.
// Genome in a Bottle) would compare CALL_VARIANTS' output VCFs per
// sample -- the same comparison src/variant_caller.py's
// compute_concordance() performs against the simulation's hidden true
// genotypes.
