# Reproducing TCGA runs

We securely reanalyzed TCGA in the cloud with [Amazon Elastic MapReduce](https://aws.amazon.com/emr/). Reproducing our TCGA runs requires [dbGaP authorized access](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v1.p1), a [CGC account](https://cgc.sbgenomics.com/), and an AWS account set up for runs of [Rail-RNA](http://rail.bio) on dbGaP-protected data. If you already have an AWS account, see [our documentation](http://docs.rail.bio/dbgap/) for instructions on preparing the account for analysis of dbGaP-protected data.

1. Clone the repo and change to this directory at the command line. Then run `python tcga_file_list.py >tcga_file_list.tsv` to obtain a list of file paths from a SPARQL query of CGC. See the docstring of `tcga_file_list.py` for its requirements. The user's list may be different from the list we obtained when we performed the query (9/29/2016). Our list is `tcga_file_list.tsv`, and the user may skip this step and simply use our file, assuming all file paths on the CGC are the same.
2. [Download](http://rail.bio) and install Rail-RNA v0.2.4a. Set Rail-RNA up for analyzing dbGaP-protected data by following the instructions at http://docs.rail.bio/dbgap/.
3. Use the Python script `gen.py` to regenerate all the Rail-RNA manifest files (`*.manifest`) in this directory as well as scripts that run Rail-RNA to preprocess (`prep_tcga_batch_*.sh`) and align (`align_tcga_batch_*.sh`) TCGA data on [Amazon Elastic MapReduce](https://aws.amazon.com/emr/). Refer to `gen.py`'s docstring for the precise command to execute; be sure to change the output bucket on S3. The script divides TCGA into 30 batches, each with about 380 randomly selected samples. A given batch is associated with a different Rail-RNA manifest file, preprocess script, and alignment script. Note that `gen.py` requires an authorization token provided by CGC. To obtain one, sign up for a [CGC account](https://cgc.sbgenomics.com/), confirm dbGaP authorized access to TCGA, generate a token using the CGC web interface, and store it in a text file locally.
4. For each batch b (a number between 0 and 29 inclusive), run

        sh prep_tcga_batch_b.sh
wait for for the Rail-RNA preprocess job on Elastic MapReduce to finish successfully, and next run

        sh align_tcga_batch_b.sh
5. Download all results from the output bucket on S3 you chose in step 3 to a dbGaP-compliant local cluster using either the [AWS CLI](https://aws.amazon.com/cli/) or the console.

# Reproducing metadata tables

Run `sh tcga_query.sh` to ultimately obtain `all_cgc_metadata.tsv.gz`.
