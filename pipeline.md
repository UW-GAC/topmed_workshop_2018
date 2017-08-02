# Analysis Pipeline

The DCC's analysis pipeline is hosted on github:
[https://github.com/smgogarten/analysis_pipeline](https://github.com/smgogarten/analysis_pipeline)

## Running on a local cluster

To run a burden test on our local SGE cluster, first we create a config file and call it `assoc_window_burden.config`:

```
out_prefix "test"
gds_file "testdata/1KG_phase3_subset_chr .gds"
phenotype_file "testdata/1KG_phase3_subset_annot.RData"
pcrelate_file "testdata/round2_pcrelate.gds"
pca_file "testdata/round2_pcair.RData"
sample_include_file "testdata/sample_include.RData"
variant_include_file "testdata/variant_include_chr .RData"
outcome outcome
covars "sex"
n_pcs 4
alt_freq_range "0 0.1"
test "burden"
test_type "score"
```

We will use the python script `assoc.py` to submit all jobs. First we look at the available options:

```
setenv PIPELINE /projects/topmed/working_code/analysis_pipeline
$PIPELINE/assoc.py --help
```

Let's run a sliding window test on chromosomes 1-10. We will also specify the cluster type, although UW_Cluster is actually the default. The last argument is our config file.

First, we print the commands that will be be run without actually submitting jobs:

```
$PIPELINE/assoc.py --chromosomes 1-10 --cluster_type UW_Cluster --print_only window testdata/assoc_window_burden.config
```

The default segment length is 10,000 kb, but we can change that to 50,000 kb when we submit:

```
$PIPELINE/assoc.py --chromosomes 1-10 --cluster_type UW_Cluster --segment_length 50000 window testdata/assoc_window_burden.config
```

We can use the `qstat` command to check the status of our jobs.


## Running on AWS Batch

To run a burden test on AWS Batch, we do the following general steps:
1. Log into the the docker AMI instance
2. cd to a working directory on our EFS volume
3. Create the configuration file `assoc_window_burden.config`
4. Optionally execute the association pipeline specifying the AWS Batch service to print out the commands (not running the pipeline)
5. Execute the association pipeline specifying the AWS Batch service to run the pipeline
6. Monitor the pipeline via the AWS Batch console

### Log into AWS docker image
ssh into our image which is running docker.  Various docker commands can be executed including running TOPMed version of R (note: TOPMed data is not part of the docker image).
```
ssh -i ~/.ssh/<some private key> kuraisa@54.244.25.94
[kuraisa@ip-172-255-46-100]~
_4816$ docker images
...
[kuraisa@ip-172-255-46-100]~
_4817$ docker run -it uwgac/r-topmed:dev /bin/bash
/# which R
...
/# R
...
> .libPaths()
...
> library(SeqArray)
...
> q()
...
/# exit
[kuraisa@ip-172-255-46-100]~
_4818$  
```
### cd to working directory and create config file
```
cd /projects/topmed/analysts/kuraisa/tm-workshop/
vi assoc_window_burden.config
...
```
### Print out AWS commands if executing the pipeline
```
python /projects/topmed/dev_code/analysis_pipeline/assoc.py \
 single ./assoc_window_burden.config \
  --cluster_type AWS_Batch --verbose \
  --cluster_file \
   /projects/topmed/dev_code/analysis_pipeline/testdata_batch.json \
  --print > single_print.log  2>&1
```
### Execute the pipeline
```
python /projects/topmed/dev_code/analysis_pipeline/assoc.py \
 single ./assoc_window_burden.config \
  --cluster_type AWS_Batch --verbose \
  --cluster_file \
   /projects/topmed/dev_code/analysis_pipeline/testdata_batch.json \
  > burden_print.log  2>&1

```
### Monitor the jobs
From the web browser, log into the AWS account and select the **Batch Services** to monitor:
- Summary via **Dashboard**
- Job queue **Optimal_topmed_testdata**
- View high-level job logs

You can switch to **ec2 services** to monitor instances being created or running to support the various jobs.
