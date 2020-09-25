#!/usr/bin/env bash
# Wrapper script to:
# 1. Get nextflow.
# 2. Create a kubernetes cluster on google cloud.
# 3. Run the pipeline there.
# 4. Shutdown the cluster.

# Define environment variables to install experimental
# nextflow with google cloud compatibility
export NXF_VER=20.04.1
export NXF_MODE=google

# Get nextflow executable if not already present
if [ ! -f './nextflow' ]; then
        wget -qO- https://get.nextflow.io | bash
fi

# Fetch cluster name in config file
BUCKET_NAME=$(grep 'bucket_name *= *' nextflow.config | tr -d " \"\'")
BUCKET_NAME=${BUCKET_NAME//*bucket_name:/}

if [ -z ${BUCKET_NAME} ]; then
  echo "Error: bucket_name not defined in config.yaml"
  exit 1
fi
# Run through kubernetes and use data on the GS bucket
./nextflow cloud create nf-k8s-vapormap-cluster -c $CLUSTER_SIZE
./nextflow kuberun vapormap.nf -v gs://bucket_name/mount/path
./nextflow cloud shutdown nf-k8s-vapormap-cluster
