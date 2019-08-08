from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

bucket = 'test_hbv_cancer'
GS = GSRemoteProvider()
rule cat_file:
  input: GS.remote(bucket + "/data/input/centro_hg19.bed")
  output: GS.remote(bucket + "/data/output/centro_copied.bed")
  singularity: "docker://ubuntu:18.04"
  shell: "cat {input} > {output}"
