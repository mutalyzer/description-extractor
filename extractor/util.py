"""
General utility definitions.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys


PY2 = sys.version_info[0] == 2


# From BioPython.
protein_letters_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp',
    'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met',
    'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp',
    'Y': 'Tyr',
}
protein_letters_1to3_extended = dict(list(protein_letters_1to3.items()) + list({
    'B': 'Asx', 'X': 'Xaa', 'Z': 'Glx', 'J': 'Xle',
    'U': 'Sel', 'O': 'Pyl',
}.items()))


# From BioPython.
def seq3(seq, custom_map={'*': 'Ter'}, undef_code='Xaa'):
    """Turn a one letter code protein sequence into one with three letter codes.

    The single input argument 'seq' should be a protein sequence using single
    letter codes, either as a python string or as a Seq or MutableSeq object.

    This function returns the amino acid sequence as a string using the three
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters B for "Asx", J for "Xle" and X for "Xaa", and also U
    for "Sel" and O for "Pyl") plus "Ter" for a terminator given as an asterisk.
    Any unknown character (including possible gap characters), is changed into
    'Xaa'.

    e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq3("MAIVMGRWKGAR*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer'

    You can set a custom translation of the codon termination code using the
    "custom_map" argument, e.g.

    >>> seq3("MAIVMGRWKGAR*", custom_map={"*": "***"})
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArg***'

    You can also set a custom translation for non-amino acid characters, such
    as '-', using the "undef_code" argument, e.g.

    >>> seq3("MAIVMGRWKGA--R*", undef_code='---')
    'MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer'

    If not given, "undef_code" defaults to "Xaa", e.g.

    >>> seq3("MAIVMGRWKGA--R*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaXaaXaaArgTer'

    This function was inspired by BioPerl's seq3.
    """
    # not doing .update() on IUPACData dict with custom_map dict
    # to preserve its initial state (may be imported in other modules)
    threecode = dict(list(protein_letters_1to3_extended.items()) +
                     list(custom_map.items()))
    #We use a default of 'Xaa' for undefined letters
    #Note this will map '-' to 'Xaa' which may be undesirable!
    return ''.join(threecode.get(aa, undef_code) for aa in seq)


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
