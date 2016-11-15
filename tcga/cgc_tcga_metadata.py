#!/usr/bin/env python
"""
cgc_tcga_metadata.py

Gets all TCGA metadata using SPARQL queries and writes results to stdout in
TSV format.

This file is based on the Seven Bridges example IPython notebook at 
https://github.com/sbg/docs/blob/master/cgc/SPARQL/
SPARQL_download_notebook.ipynb.

Requires SPARQLWrapper, which can be installed with `pip install SPARQLWrapper`

We ran this script on Tuesday, November 15 at 1:36 PM CDT as

python cgc_tcga_metadata.py | gzip -9 >cgc_tcga_metadata.tsv.gz

to obtain cgc_tcga_metadata.tsv.gz in the same directory as this script.
"""
from SPARQLWrapper import SPARQLWrapper, JSON
import json
import csv
import sys

if __name__ == '__main__':
    # Get query type
    query_type = sys.argv[1]
    query_type_lower = sys.argv[1].lower()
    # Use the public endpoint
    sparql_endpoint = (
        'https://opensparql.sbgenomics.com/blazegraph/namespace/'
        'tcga_metadata_kb/sparql'
    )

    # Initialize the SPARQL wrapper with the endpoint
    sparql = SPARQLWrapper(sparql_endpoint)

    props = (
"""PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>

SELECT distinct ?p
WHERE
{{
  ?{query_type_lower} a tcga:{query_type} .
  ?{query_type_lower} ?p ?o
}}
"""
    ).format(query_type=query_type, query_type_lower=query_type_lower)

sparql.setQuery(props)

sparql.setReturnFormat(JSON)
results = sparql.query().convert()
hases = [result['p']['value'].rpartition('#')[2]
          for result in results['results']['bindings']
          if 'www.w3.org' not in result['p']['value']]

data = {}
fields = set()

for has in hases:
    after_has = has[3:]
    underscored = ''.join([('_' + letter.lower() if letter.isupper() 
                                and (i > 0
                                     and after_has[i-1].islower())
                                and i != 0
                                else letter.lower())
                            for i, letter in enumerate(after_has)])
    fields.add(underscored)
    current_query = (
"""
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>
 
SELECT ?{query_type_lower}_label ?{underscored}
WHERE
{{
  ?{query_type_lower} a tcga:{query_type} .
  ?{query_type_lower} rdfs:label ?{query_type_lower}_label .
  ?{query_type_lower} tcga:{has} ?{underscored}
}}
"""
    ).format(
        query_type_lower=query_type_lower,
        query_type=query_type,
        underscored=underscored,
        has=has
    )
    sparql.setQuery(current_query)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    '''if not results['results']['bindings']:
        # No RDFS label; try different query
        current_query = (
"""
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>
 
SELECT ?{query_type_lower}_label ?{underscored}
WHERE
{{
  ?case a tcga:Case .
  ?case rdfs:label ?{query_type_lower}_label .
  ?case tcga:{has} ?{underscored}
}}
LIMIT 10
"""
        ).format(
            query_type_lower=query_type_lower,
            underscored=underscored,
            has=has
        )'''
    label = '{query_type_lower}_label'.format(
                                query_type_lower=query_type_lower
                            )
    other_label = [el for el in results['head']['vars'] if el != label][0]
    for result in results['results']['bindings']:
        if result[label]['value'] not in data:
            data[result[label]['value']] = {}
        data[result[label]['value']][underscored] = result[
                                                        other_label
                                                    ]['value']
print data
"""for result in results['results']['bindings']:
    print '\t'.join(
            [result[label]['value'] for label in
                ['uuid', 'path', 'sample_details', 'disease']]
        )

print '\t'.join(
        ['GDC_UUID', 'CGC_path', 'sample_type', 'disease']
    )

query_prototype = (
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>
 
SELECT ?case_label ?{}
WHERE
{
 ?case a tcga:Case .
 ?case rdfs:label ?case_label .
 ?case tcga:hasAgeAtDiagnosis ?age
}

        )

for result in results['results']['bindings']:
    print '\t'.join(
            [result[label]['value'] for label in
                ['uuid', 'path', 'sample_details', 'disease']]
        )"""
