#!/usr/bin/env nextflow

//Get the modules for each pipeline
include { trimPrimers; mapReads; makeConsensus } from './modules/illumina.nf'



workflow illumina_wf{
    //Define the input channels
    inFiles_ch = Channel.fromFilePairs("${params.fastq}*{R1,R2}*fastq.gz")
    inPrimers_ch = Channel.value('/home/dmmalone/RSV_analysis/illumina_nf_testing/resources/RSV_primers_rev.fasta')
    //For the moment, the easiest way to deal with the choice of RSV A or B ref is to require a specific flag. Could move this to conditional on reads themselves in the future.
    if (params.subtype == "RSVA") {
        inRef_ch = Channel.value('/home/dmmalone/RSV_analysis/illumina_nf_testing/resources/RSVA.reference.fasta')
    } else if(params.subtype == "RSVB") {
        inRef_ch = Channel.value('/home/dmmalone/RSV_analysis/illumina_nf_testing/resources/RSVB.reference.fasta')
    } else {
        throw new IllegalArgumentException('Please choose --subtype="RSVA" or "RSVB"')
    }
    //inFiles_ch.view()
    //Work on the input files
    trimmed_ch = trimPrimers(inFiles_ch, inPrimers_ch)
    //trimmed_ch.view()
    //Before mapping for the first time bwa needs some indexes of the ref, need to make sure the pipeline accounts for that, either always index or check "exists" and run.
    //Not quite so straightforward as they need to be ingested for the pipeline to work - for now ust index immediately before.
    mapped_ch = mapReads(trimmed_ch, inRef_ch)
    //mapped_ch.view()
    makeConsensus(mapped_ch)
}

workflow {
    illumina_wf()
}
