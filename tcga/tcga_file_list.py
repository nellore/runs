#!/usr/bin/env python
"""
tcga_file_list.py

Gets paths to all TCGA RNA-seq samples on CGC and writes results to stdout
in TSV format.

This file is based on the Seven Bridges example IPython notebook at 
https://github.com/sbg/docs/blob/master/cgc/SPARQL/
SPARQL_download_notebook.ipynb.

Requires SPARQLWrapper, which can be installed with `pip install SPARQLWrapper`

We ran this script on Thursday, September 29 at 1:36 PM CDT as

python tcga_file_list.py >tcga_file_list.tsv

to obtain tcga_file_list.tsv in the same directory as this script.
"""

from SPARQLWrapper import SPARQLWrapper, JSON
import json
import csv

if __name__ == '__main__':
    # Use the public endpoint
    sparql_endpoint = (
        'https://opensparql.sbgenomics.com/blazegraph/namespace/'
        'tcga_metadata_kb/sparql'
    )

    # Initialize the SPARQL wrapper with the endpoint
    sparql = SPARQLWrapper(sparql_endpoint)

    query = (
"""PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>

SELECT distinct ?uuid ?path ?sample_details ?disease ?file_name ?xs_label ?subtype_label
WHERE
{
  ?file a tcga:File .
  ?file rdfs:label ?file_name .
     
  ?file tcga:hasDiseaseType ?dt.
  ?dt rdfs:label ?disease .
  
  ?file tcga:hasExperimentalStrategy ?xs .
  ?xs rdfs:label ?xs_label .
  filter(?xs_label="RNA-Seq") .

  ?file tcga:hasDataSubtype ?subtype .
  ?subtype rdfs:label ?subtype_label .
  filter(?subtype_label="Unaligned reads") .

  ?file tcga:hasSample ?sample .
  ?sample tcga:hasSampleType ?sample_type .
  ?sample_type rdfs:label ?sample_details .
  
  ?file tcga:hasStoragePath ?path .
  
  ?file tcga:hasGDCFileUUID ?uuid .
}
"""
        )

sparql.setQuery(query)

sparql.setReturnFormat(JSON)
results = sparql.query().convert()
print '\t'.join(
        ['GDC_UUID', 'CGC_path', 'sample_type', 'disease']
    )

for result in results['results']['bindings']:
    print '\t'.join(
            [result[label]['value'] for label in
                ['uuid', 'path', 'sample_details', 'disease']]
        )
