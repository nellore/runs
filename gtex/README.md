If anything below doesn't make sense, ask us questions: [![Join the chat at https://gitter.im/nellore/runs](https://badges.gitter.im/nellore/runs.svg)](https://gitter.im/nellore/runs?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) .

## Redoing GTEx Rail-RNA runs on [Amazon Elastic MapReduce](https://aws.amazon.com/elasticmapreduce/)

1. Install Rail-RNA v0.2.1, which is available for download [here](https://github.com/nellore/rail/raw/master/releases/install_rail-rna-0.2.1).
2. Follow the instructions [here](http://docs.rail.bio/dbgap/) to set up an [Amazon Web Services (AWS)](https://aws.amazon.com/) account with an [Identity and Access Management (IAM)](https://aws.amazon.com/iam/) user configured to analyze dbGaP data securely. Name the  CloudFormation stack `dbgap-1` rather than `dbgap`, as those instructions recommend.  The secure bucket name created with the CloudFormation template is referenced as `s3://gtex-bucket` here.
3. To faciliate submitting job flows in multiple availability zones of the US Standard region (i.e., `us-east-1`), create three more CloudFormation stacks as described [here](http://docs.rail.bio/dbgap/#create-a-secure-cloudformation-stack-administrator), except with [this](https://raw.githubusercontent.com/nellore/rail/master/src/cloudformation/dbgap_minus_cloudtrail.template) CloudFormation template, which requires specification of an availability zone for the public subnet into which the [Rail-RNA](http://rail.bio) Elastic MapReduce cluster will be launched. Choose from among `us-east-1a`, `us-east-1b`, `us-east-1c`, `us-east-1d`, and `us-east-1e`, and name the stacks `dbgap-2`, `dbgap-3`, and `dbgap-4`.
4. Download the [dbGaP repository key](http://www.ncbi.nlm.nih.gov/books/NBK63512/#Download.are_downloaded_files_encrypted) granting access to GTEx data. It should have the extension `.ngc` and is referenced as `/path/to/dbgap/key.ngc` here.
5. Run

        python gen.py
          --s3-bucket s3://gtex-bucket
          --prep-stack-names <one or more of the dbgap-* stack names above separated by spaces>
          --align-stack-names <one or more of the dbgap-* stack name above separated by spaces>
          --dbgap-key /path/to/dbgap/key.ngc
  
  to generate scripts for preprocessing and aligning GTEx data. 60 scripts representing a partitioning of GTEx RNA-seq data into 30 batches are generated: 30 for preprocessing and 30 for aligning.
6. Run the scripts generated in the previous step to submit job flows to Elastic MapReduce. Each `prep_gtex_batch_k.sh` file for `k` between `0` and `29` inclusive should be run and its job flow completed before the corresponding `align_gtex_batch_k.sh` is run to align data preprocessed and uploaded to S3. It is recommended that only three preprocessing job flows are submitted at a time. Tweak shell scripts to change the argument of `--stack-name` in the `rail-rna` command as necessary if Elastic MapReduce complains that there aren't enough IPs in the subnet of a VPC in a given availability zone to launch more job flows.
7. Use [this script](https://github.com/nellore/runs/blob/master/gtex/download.sh) to download all results from S3 to local storage. Command-line parameters are described in its comments.
8. Compute total number of reads across samples with the script [`total.sh`](total.sh). Its only command-line parameter is the local GTEx output directory specified in the previous step, which is where all analysis results should have been dumped. The figure we obtained was `896,466,227,499`. Here, "read" refers to a mate for paired-end samples.

## Reproducing figures from the Rail-dbGaP paper

### Figure 1: security architecture
This figure was generated using Keynote; see [`security_figure.key`](security_figure.key).

### Figures 2 and 3: core activity and costs
Run the Mathematica 10 notebook [`rail_dbgap_plots.nb`](rail_dbgap_plots.nb). It uses [`costs.csv`](costs.csv), costs downloaded from the AWS Cost Explorer, as well as [`activity.tsv`](activity.tsv), which has start and end times of all GTEx preprocess and align job flows. The file `activity.tsv` was generated with [`reconstruct_activity.py`](reconstruct_activity.py) from the saved Elastic MapReduce web interface HTML files in [`logs/`](logs/). If you don't have Mathematica, check [`rail_dbgap_plots.pdf`](rail_dbgap_plots.pdf) for its output.
