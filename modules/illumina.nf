process trimPrimers {
    container "${params.wf.illumina_container}@${params.wf.illumina_container_sha}"
    publishDir "results/${sample_ID}/trimmed_reads_1", pattern: "*.trimmed.fq"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_files)
    path primer_list

    output:
    tuple val(sample_ID), path("*.trimmed.fq")

    script:
    """
    bbduk.sh in=${sample_ID_files[0]} in2=${sample_ID_files[1]} out=${sample_ID_files[0]}.trimmed.fq out2=${sample_ID_files[1]}.trimmed.fq ref=${primer_list} ktrimtips=33 k=21 mink=5 hdist=1 rcomp=f
    """
}

process mapReads {
    container "${params.wf.illumina_container}@${params.wf.illumina_container_sha}"
    publishDir "results/${sample_ID}/mapped_reads_2", pattern: "*.sorted.bam"
    
    debug true

    input:
    tuple val(sample_ID), path(sample_ID_trimmed)
    path ref_file

    output:
    tuple val(sample_ID), path("*.sorted.bam")

    script:
    """
    bwa index ${ref_file}
    bwa mem ${ref_file} ${sample_ID_trimmed} | samtools sort -o ${sample_ID}.sorted.bam
    """
}

process makeConsensus {
    container "${params.wf.illumina_container}@${params.wf.illumina_container_sha}"
    publishDir "results/consensus_sequences", pattern: "*consensus.fa"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_mapped)

    output:
    tuple val(sample_ID), path("*.consensus.fa")

    script:
    """
    samtools mpileup -aa -A -B -d 0 -Q 0 ${sample_ID_mapped} | ivar consensus -t 0.75 -m 10 -n N -p ${sample_ID}.consensus
    """
}