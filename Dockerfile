FROM continuumio/miniconda3:4.8.2
LABEL authors="Cyril Matthey-Doret" \
      description="Docker image containing all software requirements for cmdoret/vapormap read alignment pipeline"

RUN conda config --add channels bioconda
RUN conda install -c conda-forge -y pip bowtie2 bwa samtools seqkit \
  && conda clean -afy

RUN pip install hicstuff

