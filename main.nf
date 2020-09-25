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
        tuple sample, file("fq_splits/${sample}*.fq.gz") into fq_splits
        
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
fq_splits = fq_splits.transpose().view()

// Align each fastq split independently
process mapSplit{
        tag "map $split"

        input:
        tuple sample, file(fq) from fq_splits

        output:
        tuple sample, file("${split}.bam")into bam_splits

        script:
        split = fq.toString() - ~/\.f(ast)?q(\.gz)?$/

        if ( params.mode == 'iterative' )
                """
                mkdir -p bam_splits
                hicstuff iteralign -t ${task.cpus} -g $ref $fq -o ${split}.bam
                """
        else
                """
                mkdir -p bam_splits
                bowtie2 -p ${task.cpus} -x $ref -U $fq -S - \
                | samtools sort -n -o ${split}.bam
                """
}

// Remove duplicated entries to get each individual sample 
// name and corresponding alignment directory
bam_groups = bam_splits.groupTuple().view()

// Merge all split alignment files from a sample and publish the
// result to the output directory
process mergeBams{
        tag "merge $sample"
        publishDir "${params.output_dir}"
        input:
        tuple sample, file(bam_splits) from bam_groups

        output:
        file "${sample}.bam" into bam_merged

        script:
        """
        samtools merge -n -@ ${task.cpus} ${sample}.bam ${bam_splits}
        """
}

