#!/usr/bin/env nextflow

// Should be the index
ch_genome = file(params.reference, checkIfExists: true)
ch_fastqs = Channel.fromPath(params.input_dir + "/*").filter(~/.*(fastq|fq)(.gz)?$/)
//split_names = [f'part_{s:03}' for s in range(1, params.n_chunks + 1)]

//log.info "Genome is: $ch_genome"
//log.info ch_fastqs.println { "Reads are: $it" }

process splitFastq{
        tag "split_$sample"
	container 'cmdoret/seqkit:latest'
	publishDir "${params.output_dir}/${sample}/"

        input:
        val fastq from ch_fastqs
        
        output:
        val [$sample]*${params.n_chunks} into ch_samples 
        file "splits/*.fq.gz" into fq_splits
        
        script:
        sample = fastq.baseName.toString() - ~/.(fastq|fq)(.gz)?$/


        """
        seqkit split2 -p ${params.n_chunks} \
                      -w 0 \
                      -f \
                      -1 $fastq \
                      -O "splits/"
        """
}

process mapSplit{

	publishDir "${params.output_dir}/${sample}"

        input:
        file index from ch_genome
        val sample from ch_samples
        file fq from fq_splits

        output:
        file "split_*.bam" into bam_splits

        script:
        if ( params.mode == 'iterative' )
                """
                hicstuff iteralign -g $index $fq -o $bam_splits
                """
        else
                """
                bowtie2 -x $index -U $fq -S $bam_splits
                """
}

process mergeBams{

	publishDir "${params.output_dir}"
        input:
        file bam_split from bam_splits

        output:
        file "sample.bam" into bam_merged

        """
        samtools merge $bam_merged $sample_*
        """
}

