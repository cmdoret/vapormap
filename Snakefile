from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

bucket = 'smk-k8s-demo-bucket'
GS = GSRemoteProvider()
rule cat_file:
  input: GS.remote(bucket + "/test.txt")
  output: GS.remote(bucket + "//test.txt")
  shell: "cat {input} > {output}"
