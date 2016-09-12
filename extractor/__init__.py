"""
extractor: Extract a list of differences between two sequences.


Copyright (c) 2013 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2016 Jonathan K. Vis <j.k.vis@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import describe, extractor


__version_info__ = tuple(extractor.VERSION.split('.'))


__version__ = extractor.VERSION
__author__ = 'LUMC, Jonathan K. Vis'
__contact__ = 'j.k.vis@lumc.nl'
__homepage__ = 'https://github.com/mutalyzer/description-extractor'


describe_dna = describe.describe_dna
describe_protein = describe.describe_protein
describe_repeats = describe.describe_repeats
extract = extractor.extract
