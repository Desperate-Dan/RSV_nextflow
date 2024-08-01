#!/usr/bin/env nextflow

//Get the modules for each pipeline
include { trimPrimers; mapReads; makeConsensus } from './modules/illumina.nf'
include { ampliClean; articMinion } from './modules/nanopore.nf'

//These lines for fastq dir parsing are taken from rmcolq's workflow https://github.com/rmcolq/pantheon
EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}

workflow nanopore_wf {
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

workflow illumina_wf {
    //Define the input channels
    inFiles_ch = Channel.fromFilePairs("${params.fastq}*{R1,R2}*fastq.gz")
    inPrimers_ch = Channel.value("${params.illumina_bed}")
    //For the moment, the easiest way to deal with the choice of RSV A or B ref is to require a specific flag. Could move this to conditional on reads themselves in the future.
    if (params.subtype == "RSVA") {
        inRef_ch = Channel.value("${params.rsva_ref}")
    } else if(params.subtype == "RSVB") {
        inRef_ch = Channel.value("${params.rsvb_ref}")
    } else {
        throw new IllegalArgumentException('Please choose --subtype="RSVA" or "RSVB"')
    }
    //Work on the input files
    trimmed_ch = trimPrimers(inFiles_ch, inPrimers_ch)
    //Before mapping for the first time bwa needs some indexes of the ref, need to make sure the pipeline accounts for that, either always index or check "exists" and run.
    //Not quite so straightforward as they need to be ingested for the pipeline to work - for now just index immediately before.
    mapped_ch = mapReads(trimmed_ch, inRef_ch)
    makeConsensus(mapped_ch)
}

workflow {
    if (params.illumina) 
        illumina_wf()
    else
        nanopore_wf()
}