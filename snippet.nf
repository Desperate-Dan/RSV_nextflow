inFiles = Channel.fromFilePairs('./reads/*{R1,R2}*fastq.gz')


process LETS_DO_A_THING {
    publishDir "results/${sample_ID}", pattern: "*.file" // Can use a pattern match with {} to match multiple things eg *.{bam,bai}

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_paths)

    output:
    tuple val(sample_ID), path("*.file")

    script:
    """
    echo the sample is called $sample_ID and contains these files: $sample_ID_paths > ${sample_ID}.file
    """
}

process ONE_MORE_THING {
    debug true

    input:
    tuple val(sample_ID), path(new_file)

    script:
    """
    cat ${new_file}
    """
}

workflow {
    out_ch = LETS_DO_A_THING(inFiles)
    ONE_MORE_THING(out_ch)
}