#!/bin/bash
# Assuming google-cloud-sdk is installed and setup, and a project is set up
# Also assuming you have a file named test.txt inside your google storage  bucket
NODES=1
CLUSTER_NAME=smk-k8s-demo-cluster
BUCKET_NAME=smk-k8s-demo-bucket

# 1: setup kubernetes cluster
gcloud container clusters create "$CLUSTER_NAME" --num-nodes="$NODES" --scopes storage-rw --region eu-west1-d
gcloud container clusters get-credentials "$CLUSTER_NAME"

# 2: log into google storage
gcloud auth application-default login

# 3: Run snakemake through kubernetes and use data on the GS bucket
snakemake --default-remote-provider GS \
          --default-remote-prefix "$BUCKET_NAME" \
          --keep-remote \
          --use-singularity \
          --kubernetes
