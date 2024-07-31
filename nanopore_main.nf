//Untested at the moment!
include { ampliClean; articMinion } from './modules/nanopore.nf'

//These lines for fastq dir parsing are taken from rmcolq's workflow https://github.com/rmcolq/pantheon
EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}

workflow {
//Define input channels  
  ref_ch = file("${params.refs}")
  bed_ch = file("${params.bed}")
  schemes_dir_ch = file("${params.schemes_dir}")
  min_ch = Channel.value("${params.min}")
  max_ch = Channel.value("${params.max}")
  med_mod_ch = Channel.value("${params.medaka_model}")
//These lines for fastq dir parsing are taken from rmcolq's workflow https://github.com/rmcolq/pantheon
  run_dir = file("${params.fastq}", type: "dir", checkIfExists:true)
  barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}

//Run the processes
  ampliClean(barcode_input, ref_ch, bed_ch, min_ch, max_ch)
  articMinion(ampliClean.out.reads, schemes_dir_ch, med_mod_ch)
}
