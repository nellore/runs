# Reproducing TCGA runs

We securely reanalyzed TCGA in the cloud with [Amazon Elastic MapReduce](https://aws.amazon.com/emr/). Reproducing our TCGA runs requires [dbGaP authorized access](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v1.p1), a [CGC account](https://cgc.sbgenomics.com/), and an AWS account set up for runs of [Rail-RNA](http://rail.bio) on dbGaP-protected data. If you already have an AWS account, see [our documentation](http://docs.rail.bio/dbgap/) for instructions on preparing the account for analysis of dbGaP-protected data.

1. Clone the repo and change to this directory at the command line. Then run `python tcga_file_list.py >tcga_file_list.tsv` to obtain a list of file paths from a SPARQL query of CGC. See the docstring of `tcga_file_list.py` for its requirements. The user's list may be different from the list we obtained when we performed the query (9/29/2016). Our list is `tcga_file_list.tsv`, and the user may skip this step and simply use our file, assuming all file paths are the same.
2. 
