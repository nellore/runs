#!/usr/bin/env python
"""
mismatch_mean.py

Obtains mean of mean mismatch bigwigs. Considers only those samples for which
coverage AUC > 0 when weighting mean mismatch bigwigs for a given batch. For
example, if there were only two batches, and the first batch had 5 bigwigs with
nonzero AUC while the second batch had 6 bigwigs with nonzero AUC, the mean
bigwig would be 5 / 11 * (first batch bigwig) + 6 / 11 * (second batch bigwig).
"""
import os

if __name__ == '__main__':
	import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--batch-dir', type=str, required=True,
            help=('root directory of all-of-SRA results; this is the '
            	  'directory with the batch_* subdirectories')
        )
    for path