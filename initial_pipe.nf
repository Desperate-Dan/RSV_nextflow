inFiles_ch = Channel.fromFilePairs('./reads/*{R1,R2}*fastq.gz')
inPrimers_ch = Channel.value('/home/dmmalone/RSV_analysis/illumina_nf_testing/resources/RSV_primers_rev.fasta')
inRef_ch = Channel.value('/home/dmmalone/RSV_analysis/illumina_nf_testing/resources/RSVA.reference.fasta')

process trimPrimers {
    publishDir "results/${sample_ID}/trimmed_reads_1", pattern: "*.trimmed.fq" // Can use a pattern match with {} to match multiple things eg *.{bam,bai}

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
    publishDir "results/consensus_sequences", pattern: "*consensus.fa"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_mapped)

    output:
    tuple val(sample_ID), path("*.consensus.fa")

    script:
    """
    samtools mpileup -aa -A -B -d 0 -Q0 ${sample_ID_mapped} | ivar consensus -t 0.75 -m 10 -n N -p ${sample_ID}.consensus
    """
}

workflow {
    inFiles_ch.view()
    trimmed_ch = trimPrimers(inFiles_ch, inPrimers_ch)
    trimmed_ch.view()
    //Before mapping for the first time bwa needs some indexes of the ref, need to make sure the pipeline accounts for that, either always index or check "exists" and run.
    //Not quite so straightforward as they need to be ingested for the pipeline to work - for now ust index immediately before.
    mapped_ch = mapReads(trimmed_ch, inRef_ch)
    mapped_ch.view()
    makeConsensus(mapped_ch)
}