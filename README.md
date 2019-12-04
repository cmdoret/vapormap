# Vapormap

Vapormap is a cloud-based alignment pipeline for Hi-C. Each fastq file is split into a number of chunks which are all aligned in parallel on the cloud.

Hi-C reads can be aligned using iterative mapping. All configuration is done through config.yaml.
The `gcloud_setup` script takes care of creating the kubernetes cluster, but input fastq files should have been uploaded on a google-cloud bucket whose name must be defined in the config.yaml file. All paths defined in the config are defined inside the gcloud bucket.

### Installation

You should have installed and configured the gcloud sdk for your account. You also need python 3.6 You can install python dependencies using:

```bash
pip3 install -Ur requirements.txt
```

### Usage

Setup the cluster using:

```bash
./gcloud_setup "cluster-name"
```

Run the alignment using:

```bash
./vapormap 'bucket-name'
```
Where `bucket-name` matches an existing bucket on your gcloud account where files have already been uploaded.

Don't forget to delete the cluster once you're finished, using:

```bash
gcloud container clusters delete "cluster-name"
```
