#!/usr/bin/env python
"""
reconstruct_activity.py
Abhi Nellore / December 14, 2015

Reads saved EMR console HTML files in log/ to create TSV with cluster activity.
We ran

cat <(echo -e "batch number\tflow type\tcreation date\tend date\tnumber of instances\tinstance type\tnumber of active worker cores") <(python reconstruct_activity.py | sort -k1,1n -k2,2r) >activity.tsv

and added two columns, "exclude" and "notes", to activity.tsv .
"""
import glob
import os

if __name__ == '__main__':
    for log in glob.glob(os.path.join(
                os.path.dirname(os.path.realpath(__file__)), 'logs', '*.html'
            )):
        with open(log) as log_stream:
            f = log_stream.read()
            dates = []
            flow = None
            for k in xrange(len(f)):
                if f[k:k+13] == 'Creation date':
                    for l in xrange(k+13, k+500):
                        if f[l:l+5] == '2015-':
                            dates.append(f[l:l+24])
                            break
                elif f[k:k+8] == 'End date':
                    for l in xrange(k+13, k+500):
                        if f[l:l+5] == '2015-':
                            dates.append(f[l:l+24])
                            break
                elif f[k:k+30] == 's3://dbgap-stack-361204003210/':
                    flow = f[k+30:k+300].partition('.')[0]
                    flow_type = flow.split('_')[1]
                    batch_number = flow.split('_')[-1]
            assert flow is not None
            print '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        batch_number, flow_type, dates[0], dates[1],
                        81 if flow_type == 'align' else 21,
                        'c3.8xlarge' if flow_type == 'align' else 'm3.xlarge',
                        80*32 if flow_type == 'align' else 20*4
                    )