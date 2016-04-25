#!/usr/bin/env python
"""
figshare.py

Command-line tool that uploads file(s) to figshare. Based heavily on
https://docs.figshare.com/api/upload_example/. Requires requests library.
"""
import json
import os
import requests
import sys

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--token', type=str, required=True,
            help='access token'
        )
    parser.add_argument('--paths', type=str, nargs='+', required=True,
            help='space-separated list of paths to files to upload'
        )
    parser.add_argument('--article-id', type=str, required=True,
            help='article id under which files are to be uploaded'
        )
    args = parser.parse_args()

    BASE_URL = 'https://api.figshare.com/v2/{endpoint}'
    HEADERS = {'Authorization': 'token ' + args.token}

    for file_to_upload in args.paths:
        file_name = os.path.basename(file_to_upload)
        # Get file info
        with open(file_to_upload, 'rb') as fin:
            fin.seek(0, 2)  # Go to end of file
            size = fin.tell()
        data = json.dumps({'name': file_name, 'size': size})

        # Initiate upload
        endpoint = 'account/articles/{}/files'.format(args.article_id)
        resp = requests.post(
                    BASE_URL.format(endpoint=endpoint), headers=HEADERS, 
                    data=data
                )

        file_id = json.loads(resp.content)['location'].rsplit('/', 1)[1]

        # Get upload/parts info
        endpoint = 'account/articles/{}/files/{}'.format(args.article_id,
                                                            file_id)
        resp = requests.get(
                        BASE_URL.format(endpoint=endpoint), headers=HEADERS
                    )

        url = '{upload_url}'.format(**json.loads(resp.content))
        parts = json.loads(requests.get(url).content)['parts']

        # Upload parts
        with open(file_to_upload, 'rb') as fin:
            for part in parts:
                size = part['endOffset'] - part['startOffset'] + 1
                address = '{}/{}'.format(url, part['partNo'])
                requests.put(address, data=fin.read(size))

        # Mark file upload as completed
        requests.post(BASE_URL.format(endpoint=endpoint), headers=HEADERS)

    # List files
    endpoint = 'account/articles/{}/files'.format(args.article_id)
    resp = requests.get(BASE_URL.format(endpoint=endpoint), headers=HEADERS)
    print >>sys.stderr, 'Files uploaded follow.'
    print >>sys.stderr, resp.content
