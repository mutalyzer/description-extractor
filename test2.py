#!/usr/bin/env python

from __future__ import unicode_literals

import monoseq
from suds.client import Client

from extractor import describe


URL = 'https://mutalyzer.nl/services/?wsdl'


def add_annotation(dest, start, end):
    if len(dest) == 2:
        if len(dest[0]) == len(dest[1]):
            dest[0].append((start, end))
        else:
            dest[1].append((start, end))
    else:
        dest.append([(start, end)])


client = Client(URL, cache=None)
service = client.service

response = service.runMutalyzer(
    #'NM_002001.2:c.[10del;37_38del]')
    #'NM_002001.2:c.[10del;37_38del;38_39insTGATGATGATGATGATGATGATGATGA]')
    'NM_002001.2:c.[10del;37_38delinsTGATGATGATGATGATGATGATGATG;51_52del]')

description = describe.describe_protein(
    response.origProtein, response.newProtein)

ref_annotation = []
alt_annotation = []
for rv in description:
    for i in rv.annotated_inserted:
        #ref_annotation.append([(i.start - 1, i.end)])
        #alt_annotation.append([(i.sample_start - 1, i.sample_end)])
        add_annotation(ref_annotation, i.start - 1, i.end)
        add_annotation(alt_annotation, i.sample_start - 1, i.sample_end)
        #print i.start, i.end, i.sample_start, i.sample_end
        #print i.start, i.end, rv.start, i.sample_start, i.sample_end

print '    p.{}'.format(description)
print '    q.{}'.format(description.nhgvs())

print ref_annotation
print alt_annotation

print
print monoseq.pprint_sequence(response.origProtein, ref_annotation,
    format=monoseq.AnsiFormat)
print
print monoseq.pprint_sequence(response.newProtein, alt_annotation,
    format=monoseq.AnsiFormat)
print
