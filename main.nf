#!/usr/bin/env nextflow
/* Mapping reads in parallel, on the cloud. Currently very preliminary:
 * TODO:
 *  - Add support for multiple aligner softwares
 *  - Add pair end support
*/

// Should be the index
ref = file(params.reference, checkIfExists: true)
fastqs_ch = Channel.fromPath(params.input_dir + "/*").filter(~/.*(fastq|fq)(.gz)?$/)

// Split each input fastq into a predefined number of chunks
process splitFastq{
        tag "split_$sample"
	container 'cmdoret/seqkit:latest'

        input:
        val fastq from fastqs_ch
        
        output:
        file("fq_splits/${sample}*.fq.gz") into fq_splits
        
        script:
        sample = fastq.baseName.toString() - ~/.(fastq|fq)(.gz)?$/


        """
        seqkit split2 -p ${params.n_chunks} \
                      -w 0 \
                      -j ${task.cpus} \
                      -f \
                      -1 $fastq \
                      -O "fq_splits/"
        """
}
// Unnest fastq files in the channel
fq_splits = fq_splits.flatten()

// Align each fastq split independently
process mapSplit{

        tag "map $split"

        input:
        file(fq) from fq_splits

        output:
        path "bam_splits/${split}.bam" into bam_splits
        path "bam_splits/" into bam_dir_ch
        val sample into samples_ch

        script:
        split = fq.toString() - ~/\.f(ast)?q(\.gz)?$/
        sample = fq.baseName.toString() - ~/\.part_[0-9]+.*$/

        if ( params.mode == 'iterative' )
                """
                mkdir -p bam_splits
                hicstuff iteralign -t ${task.cpus} -g $ref $fq -o bam_splits/${split}.bam
                """
        else
                """
                mkdir -p bam_splits
                bowtie2 -p ${task.cpus} -x $ref -U $fq -S - \
                | samtools view -O BAM -o bam_splits/${split}.bam
                """
}

// Remove duplicated entries to get each individual sample 
// name and corresponding alignment directory
bam_dir_ch = bam_dir_ch.unique()
samples_ch = samples_ch.unique()

// Merge all split alignment files from a sample and publish the
// result to the output directory
process mergeBams{
        tag "merge $sample"
	publishDir "${params.output_dir}"
        input:
        path bam_dir from bam_dir_ch
        val sample from samples_ch

        output:
        file "${sample}.bam" into bam_merged

        script:
        """
        samtools merge -@ ${task.cpus} ${sample}.bam ${bam_dir}/*
        """
}

