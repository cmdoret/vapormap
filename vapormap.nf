#!/usr/bin/env nextflow

params.str = 'Hello world!'
fastqs = Channel.fromPath(params.input_dir + "/*.fq*")
samples = params.input_dir.listfiles('*.fq').collect({it =~ /\(.*\)\.fq.*/})
split_names = [f'part_{s:03}' for s in range(1, params.n_chunks + 1)]
conda: 'envs/hic_processing.yaml'

samples, = glob_wildcards(GS.remote(str(bucket / IN / "{sample}.fq.gz")))

process indexGenome{
        container 'koszullab/hicstuff:latest'

        input:
        file genome from params.reference

        output:
        file params.genome + "*bt2" into index

        """
        bowtie2-build $genome genome_idx
        """
}

process splitFastq{
        
        input:
        val fastq from fastqs
        
        output:
        file "split_*.fq.gz" into fq_splits

        """
        seqkit split2 -p {params.n_splits} \
                      -w 0 \
                      -f \
                      -1 $sample \
                      -O {params.split_dir}
        printf ${params.str} | split -b 6 - chunk_
        """
}

process mapSplit{
        input:
        val sample from samples
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
        input:
        file bam_split from bam_splits

        output:
        file "sample.bam" into bam_merged

        """
        samtools merge $bam_merged $sample_*
        """
}

result.view { it.trim() }
