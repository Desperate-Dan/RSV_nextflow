process ampliClean {
  container "${params.wf.container}@${params.wf.container_sha}"
  
  publishDir path: "${params.out_dir}/${barcode}/ampli_clean", mode: 'copy'

  input:
    tuple val(barcode), path(binned_reads)
    path refs
    path bed
    val min
    val max
    
    
  output:
    tuple val("${barcode}"), path("${barcode}.*.fastq.gz"), emit: reads, optional: true
    path "log.txt"

  script:
    """
    ampli_clean -f ${binned_reads} -r ${refs} -o ${barcode} -b ${bed} --min ${min} --max ${max} -s --fastq --log
    """
}

process articMinion {
  container "${params.wf.container}@${params.wf.container_sha}"

  publishDir path: "${params.out_dir}/${barcode}/artic", mode: 'copy'

  input:
    tuple val(barcode), path(input_reads)
    path schemes_dir
    val (medaka_model)

  output:
    path "${barcode}.${vir}.consensus.fasta"

  script:
    vir = input_reads.name.toString().tokenize('.').get(1)
    """
    artic minion --medaka --threads ${task.cpus} --scheme-directory ${schemes_dir} --read-file ${input_reads} --medaka-model ${medaka_model} --strict ${vir}/V1 ${barcode}.${vir}
    """
}