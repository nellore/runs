#!/usr/bin/env python
"""
true_manifest.py

Replaces storage paths in dummy manifest files with temporary S3 URLs obtained
from the CGC API. This script takes a single command-line argument: the path to
a text file with a single line, the CGC authorization token, which is
generated at https://cgc.sbgenomics.com/account/#developer .

Based on https://github.com/sbg/docs/blob/master/cgc/
SPARQL/SPARQL_download_notebook.ipynb
"""
import uuid
import json
import requests
import sys

def api(api_url, path, auth_token, method='GET', query=None, data=None):
    """ Performs request on CGC API

        api_url: base URI of API
        path: query path
        auth_token: CGC authorization token
        method: HTTP method (GET or POST here)
        query: query string parameters
        data: POST data

        Return value: requests response dictionary
    """
    data = (
        json.dumps(data) if isinstance(data, dict)
        or isinstance(data,list) else None
    )
    base_url = api_url
    headers = { 
        'X-SBG-Auth-Token': auth_token, 
        'Accept': 'application/json', 
        'Content-type': 'application/json', 
    } 

    response = requests.request(
                  method, base_url + path, params=query,
                  data=data, headers=headers
              ) 
    print >>sys.stderr, "URL: ",  response.url
    print >>sys.stderr, "RESPONSE CODE: ", response.status_code
    response_dict = json.loads(response.content) if response.content else {} 
    response_headers = dict(response.headers)
    return response_dict

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--cgc-auth-token', type=str, required=True,
            help='path to text file containing CGC authorization token; '
        )
    args = parser.parse_args()
    base = 'https://cgc-api.sbgenomics.com/v2/'
    manifest = []
    for line in sys.stdin:
        manifest.append(line.strip().split('\t'))
    download_urls = api(
            api_url=base, auth_token=open(args.cgc_auth_token).read().strip(),
            path='action/files/get_download_url', method='POST', query=None,
            data=[tokens[0] for tokens in manifest]
        )
    for i in xrange(len(download_urls)):
        print '\t'.join([download_urls[i], manifest[i][1], manifest[i][2]])
