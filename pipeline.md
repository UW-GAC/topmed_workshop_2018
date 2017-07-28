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
