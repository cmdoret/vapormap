from snakemake.utils import validate
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import pathlib
shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

GS = GSRemoteProvider()
configfile: "config.yaml"
bucket = pathlib.Path(config['remote']['bucket_name'])
validate(config, schema="schemas/config.schema.yaml")
IN = config['input_dir']
OUT = config['output_dir']
TMP =  pathlib.Path('tmp')
GENOME = config['reference']
N_SPLITS = config['n_chunks']
# Make splits from fastq files to speed up mapping. [between 1 and 999 per fastq]
split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)] # base names of split files
conda: 'envs/hic_processing.yaml'

samples, = glob_wildcards(GS.remote(str(bucket / IN / "{sample}.fq.gz")))

rule all:
  input: GS.remote(expand(str(bucket / OUT / "{sample}.bam"), sample=samples))

rule bt2_index:
  output: GS.remote(str(bucket / TMP / 'genome.rev.2.bt2'))
  params:
    genome = GS.remote(str(bucket / GENOME)),
    idx = GS.remote(str(bucket / TMP / 'genome'))
  shell: "bowtie2-build --threads 1 {params.genome} {params.idx}"


rule split_hic_fastq:
  input: GS.remote(str(bucket / IN / '{sample}.fq.gz'))
  output: GS.remote(expand(str(bucket / TMP / 'split_reads' / '{{sample}}.{split}.fq.gz'), split=split_names))
  params:
    n_splits = N_SPLITS,
    split_dir = GS.remote(str(bucket / TMP / 'split_reads'))
  message: "Splitting {wildcards.fastq} into {params.n_splits} split fastq"
  shell:
     """
     mkdir -p {params.split_dir}
     # 100 split fastqs will be created with name pattern 00000.fq - 000100.fq
     seqkit split2 -p {params.n_splits} \
                   -w 0 \
                   -f \
                   -1 {input} \
                   -O {params.split_dir}
     """

# Iterative alignment of a single fastq split from a Hi-C sample
rule split_iter_align_hic:
  input:
    index_flag = GS.remote(str(bucket / TMP / 'genome.rev.2.bt2')),
    fq = GS.remote(str(bucket / TMP / 'split_reads' / '{sample}.{split}.fq.gz')),
  output: GS.remote(str(TMP / 'split_reads' / '{sample}.{split}.bam'))
  params:
    tmp_dir = lambda w: GS.remote(TMP / "split_reads" / f"{w.sample}_{w.split}"),
    index = GS.remote(str(TMP / 'genome')),
    iteralign_presets = config['iteralign_params']
  shell:
    """
    hicstuff iteralign {params.iteralign_presets} \
                       -T {params.tmp_dir} \
                       -g {params.index} \
                       -o {output} \
                       {input.fq}
    """

# Merge splits from individual mapping jobs into one bam file per library end
rule merge_split_alignments:
  input:
    GS.remote(expand(
      str(TMP / 'split_reads' / '{{sample}}.{split}.bam'),
      split=split_names
    ))
  output: GS.remote(str(bucket / OUT / '{sample}.bam'))
  threads: 1
  singularity: "docker://biocontainers/samtools:latest"
  shell: "samtools merge -O BAM -@ {threads} {output} {input}"
