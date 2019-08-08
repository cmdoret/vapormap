from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

bucket = 'smk-k8s-demo-bucket'
GS = GSRemoteProvider()
rule cat_file:
  input: GS.remote(bucket + "/data/input/centro_hg19.bed")
  output: GS.remote(bucket + "/data/output/centro_copied.bed")
  shell: "cat {input} > {output}"
