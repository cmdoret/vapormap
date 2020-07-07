#!/usr/bin/env nextflow

// Should be the index
ref = file(params.reference, checkIfExists: true)
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
        file("${sample}_splits/*.fq.gz") into fq_splits
        
        script:
        sample = fastq.baseName.toString() - ~/.(fastq|fq)(.gz)?$/


        """
        seqkit split2 -p ${params.n_chunks} \
                      -w 0 \
                      -f \
                      -1 $fastq \
                      -O "${sample}_splits/"
        """
}

fq_splits = fq_splits.flatten()

process mapSplit{

	publishDir "${params.output_dir}/${sample}"

        input:
        file(fq) from fq_splits

        output:
        file "${prefix}.bam" into bam_splits

        script:
        prefix = fq.toString() - ~/.f(ast)?q(\.gz)?$/

        if ( params.mode == 'iterative' )
                """
            hicstuff iteralign -g $ref $fq -o ${prefix}.bam
                """
        else
                """
                bowtie2 -x $ref -U $fq -S ${prefix}.split.bam
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

