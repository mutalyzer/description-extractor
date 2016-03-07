"""
General utility definitions.
"""


from __future__ import (absolute_import, division, print_function,
    unicode_literals)

import sys

from Bio.Data import CodonTable
from Bio.Data.IUPACData import (protein_letters_1to3,
    protein_letters_1to3_extended)
from Bio.SeqUtils import seq3

PY2 = sys.version_info[0] == 2


def codon_table_string(table_id):
    """
    Return the codon table referenced by {table_id} in compresed from. The
    result consists of a string of amino acids sorted by the codon that
    translates to them. For example, the codon 'AAG' has position 3 in the
    sorted list of codons, so its translation 'K' occurs in the third position
    of the output.

    :arg table_id: ID of a codon table.
    :type table_id: int

    :returns: String representation of code table referenced by {table_id}.
    :rtype: str
    """
    codons = CodonTable.unambiguous_dna_by_id[table_id].forward_table.items()

    codons += map(lambda x: (x, '*'),
        CodonTable.unambiguous_dna_by_id[table_id].stop_codons)

    return ''.join(map(lambda x: x[1], sorted(codons)))


def swig_str(s, ascii_only=True):
    """
    Given a unicode string, returns the representation expected by SWIG and
    its (UTF-8 encoded) length. Unless `ascii_only=False`, the string must
    contain only characters in the ASCII range.

    Unfortunately, SWIG encodes unicode strings on Python 3 (the `str` type)
    automatically as UTF-8, while it doesn't do so on Python 3 (the `unicode`)
    type. So we have to encode ourselves on Python 2. Hence this function.

    The SWIG documentation doesn't really discuss this, so this was a real
    pain to debug.

    Note that to correctly calculate the length of the resulting *char value,
    we also have to encode to UTF-8 on Python 3. Since this means decoding is
    done twice (once here for the length calculation and once by SWIG), we
    instead assume all characters are in the ACII range so we can use the
    length of the unicode string directly. This assumption can be removed by
    specifying `ascii_only=False`.

    http://www.swig.org/Doc2.0/SWIGDocumentation.html#Python_nn49
    https://github.com/swig/swig/blob/master/Lib/python/pystrings.swg
    http://comments.gmane.org/gmane.comp.programming.swig.devel/23268
    """
    if PY2 or not ascii_only:
        s_encoded = s.encode('utf-8')
        if PY2:
            return s_encoded, len(s_encoded)
        return s, len(s_encoded)

    return s, len(s)


#: Python 3 behaviour for `str` on both Python 2 and 3.
str = unicode if PY2 else str


def python_2_unicode_compatible(cls):
    """
    A decorator that defines `__unicode__` and `__str__` methods under Python
    2. Under Python 3 it does nothing.

    To support Python 2 and 3 with a single code base, define a `__str__`
    method returning unicode text and apply this decorator to the class.

    The implementation comes from django.utils.encoding.
    """
    if PY2:
        cls.__unicode__ = cls.__str__
        cls.__str__ = lambda self: self.__unicode__().encode('utf-8')

    return cls
