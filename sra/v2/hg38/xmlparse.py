#!/usr/bin/env python
"""
xmlparse.py

Parses XML in table created by get_biosample_data.sh so it's more readable.
Requires https://github.com/martinblech/xmltodict .
To install that, run "pip install xmltodict".
"""
import sys
import xmltodict

def silent_value(*args):
    """ Returns value from dictionary or NA if there's a KeyError

        args: dictionary followed by chain of keys; ex: 'a', 'hello',
            'howareyou' translates to a['hello']['howareyou']

        Return value: NA or value
    """
    try:
        return eval(args[0] + ''.join('["%s"]' % arg for arg in args[1:]))
    except KeyError:
        return 'NA'

print '\t'.join(
        ['sample', 'submission date', 'publication date', 'last update',
         'title', 'ids', 'attributes']
    )
for line in sys.stdin:
    sample, _, xml = line.strip().split('\t')
    xml = xmltodict.parse(xml)
    alt_sample_names = ';'.join([item.rpartition(';')[-1]
                        for item in [
                            ';'.join(i for el in db.items() for i in el)
                                for db
                                in xml['SampleData']['BioSample']['Ids']['Id']]
                                if 'sample name' in item.lower()])
    if not alt_sample_names: alt_sample_names = 'NA'
    try:
        attributes = ';'.join(
                [':'.join([el['@attribute_name'], el['#text']])
                    for el
                    in xml['SampleData'][
                                    'BioSample'
                                ]['Attributes']['Attribute']])
    except TypeError:
        attributes = 'NA'
    print ('\t'.join(
        [sample,
            silent_value('xml', 'SampleData',
                         'BioSample', '@submission_date'),
            silent_value('xml', 'SampleData',
                         'BioSample', '@publication_date'),
            silent_value('xml', 'SampleData',
                         'BioSample', '@last_update'),
            silent_value('xml', 'SampleData',
                              'BioSample', 'Description', 'Title'),
            alt_sample_names, attributes
        ]).encode('ascii', 'ignore')
    )