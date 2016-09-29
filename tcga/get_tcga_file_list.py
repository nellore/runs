#!/usr/bin/env python
"""
tcga_dummy_manifest.py

Creates dummy Rail-RNA manifest file with pointers to TCGA samples on CGC. This
dummy manifest is interpreted by temp_link_manifest.py to assign download
links that expire in 48 hours to each sample and write a new temporarily valid
manifest file. The new manifest file should be passed to Rail-RNA's prep job
flow in secure mode. The dummy manifest file can be passed to Rail-RNA's align
mode for alignment of preprocessed samples.

This file is based on the Seven Bridges example IPython notebook at 
https://github.com/sbg/docs/blob/master/cgc/SPARQL/
SPARQL_download_notebook.ipynb.
"""

from SPARQLWrapper import SPARQLWrapper, JSON
import json

if __name__ == '__main__':
    # Use the public endpoint
    sparql_endpoint = (
        'https://opensparql.sbgenomics.com/blazegraph/namespace/'
        'tcga_metadata_kb/sparql'
    )

    # Initialize the SPARQL wrapper with the endpoint
    sparql = SPARQLWrapper(sparql_endpoint)

    query = """
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>

    SELECT distinct ?uuid ?path ?sa_ty ?di_ty ?file_name ?xs_label ?subtype_label
    WHERE
    {
      ?file a tcga:File .
      ?file rdfs:label ?file_name .
         
      ?file tcga:hasDiseaseType ?dt.
      ?dt rdfs:label ?di_ty .
      
      ?file tcga:hasExperimentalStrategy ?xs .
      ?xs rdfs:label ?xs_label .
      filter(?xs_label="RNA-Seq") .

      ?file tcga:hasDataSubtype ?subtype .
      ?subtype rdfs:label ?subtype_label.
      filter(?subtype_label="Unaligned reads") .

      ?file tcga:hasSample ?sample .
      ?sample tcga:hasSampleType ?sample_type .
      ?sample_type rdfs:label ?sa_ty .
      
      ?file tcga:hasStoragePath ?path .
      
      ?file tcga:hasGDCFileUUID ?uuid .
    }
    """

sparql.setQuery(query)

sparql.setReturnFormat(JSON)
results = sparql.query().convert()

# From results, we grab a list of files. TCGA metadata database returns a list of filepaths. 
filelist = [result['path']['value'] for result in results['results']['bindings']]

