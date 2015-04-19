"""
Models for the description extractor.
"""

from __future__ import unicode_literals

from . import extractor


# From BioPython.
protein_letters_1to3  = {
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


WEIGHTS = {
    'subst': extractor.WEIGHT_SUBSTITUTION,
    'del': extractor.WEIGHT_DELETION,
    'ins': extractor.WEIGHT_INSERTION,
    'dup': extractor.WEIGHT_INSERTION,
    'inv': extractor.WEIGHT_INVERSION,
    'delins': extractor.WEIGHT_DELETION_INSERTION
}


class HGVSList(object):
    """
    Container for a list of sequences or variants.
    """
    def __init__(self, items=[]):
        self.items = list(items)


    def __getitem__(self, index):
        return self.items[index]


    def __bool__(self):
        return bool(len(self.items) > 0)


    def __nonzero__(self): # Python 2.x compatibility.
        return self.__bool__()


    def __unicode__(self):
        if len(self.items) > 1:
            return '[{}]'.format(';'.join(map(unicode, self.items)))
        return unicode(self.items[0])


    def append(self, item):
        self.items.append(item)


    def weight(self):
        weight = sum(map(lambda x: x.weight(), self.items))

        if len(self.items) > 1:
            return weight + (len(self.items) + 1) * extractor.WEIGHT_SEPARATOR
        return weight


class Allele(HGVSList):
    pass


class ISeqList(HGVSList):
    pass


class ISeq(object):
    """
    Container for an inserted sequence.
    """
    def __init__(self, sequence='', start=0, end=0, reverse=False,
            weight_position=1):
        """
        Initialise the class with the appropriate values.

        :arg unicode sequence: Literal inserted sequence.
        :arg int start: Start position for a transposed sequence.
        :arg int end: End position for a transposed sequence.
        :arg bool reverse: Inverted transposed sequence.
        """
        self.sequence = sequence
        self.start = start
        self.end = end
        self.reverse = reverse
        self.weight_position = weight_position

        self.type = 'trans'
        if self.sequence:
            self.type = 'ins'


    def __unicode__(self):
        if self.type == 'ins':
            return self.sequence

        if not (self.start or self.end):
            return ''

        inverted = 'inv' if self.reverse else ''
        return '{}_{}{}'.format(self.start, self.end, inverted)


    def __bool__(self):
         return bool(self.sequence)


    def __nonzero__(self): # Python 2.x compatibility.
        return self.__bool__()


    def weight(self):
        if self.type == 'ins':
            return len(self.sequence) * extractor.WEIGHT_BASE

        inverse_weight = WEIGHTS['inv'] if self.reverse else 0
        return (self.weight_position * 2 + extractor.WEIGHT_SEPARATOR +
            inverse_weight)


class DNAVar(object):
    """
    Container for a DNA variant.
    """
    def __init__(self, start=0, start_offset=0, end=0, end_offset=0,
            sample_start=0, sample_start_offset=0, sample_end=0,
            sample_end_offset=0, type='none', deleted=ISeqList([ISeq()]),
            inserted=ISeqList([ISeq()]), shift=0, weight_position=1):
        """
        Initialise the class with the appropriate values.

        :arg int start: Start position.
        :arg int start_offset:
        :arg int end: End position.
        :arg int end_offset:
        :arg int sample_start: Start position.
        :arg int sample_start_offset:
        :arg int sample_end: End position.
        :arg int sample_end_offset:
        :arg unicode type: Variant type.
        :arg unicode deleted: Deleted part of the reference sequence.
        :arg ISeqList inserted: Inserted part.
        :arg int shift: Amount of freedom.
        """
        # TODO: Will this container be used for all variants, or only genomic?
        #       start_offset and end_offset may be never used.
        self.start = start
        self.start_offset = start_offset
        self.end = end
        self.end_offset = end_offset
        self.sample_start = sample_start
        self.sample_start_offset = sample_start_offset
        self.sample_end = sample_end
        self.sample_end_offset = sample_end_offset
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.weight_position = weight_position
        self.shift = shift


    def __unicode__(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        :returns unicode: The HGVS description of the raw variant stored in
            this class.
        """
        if self.type == 'unknown':
            return '?'
        if self.type == 'none':
            return '='

        description = '{}'.format(self.start)

        if self.start != self.end:
            description += '_{}'.format(self.end)

        if self.type != 'subst':
            description += '{}'.format(self.type)

            if self.type in ('ins', 'delins'):
                return description + '{}'.format(self.inserted)
            return description

        return description + '{}>{}'.format(self.deleted, self.inserted)


    def weight(self):
        if self.type == 'unknown':
            return -1
        if self.type == 'none':
            return 0

        weight = self.weight_position
        if self.start != self.end:
            weight += self.weight_position + extractor.WEIGHT_SEPARATOR

        return weight + WEIGHTS[self.type] + self.inserted.weight()


class ProteinVar(object):
    """
    Container for a protein variant.

    """
    #NOTE: This is experimental code. It is not used at the moment.
    def __init__(self, start=0, end=0, sample_start=0, sample_end=0,
            type='none', deleted=ISeqList([ISeq()]),
            inserted=ISeqList([ISeq()]), shift=0, term=0):
        """
        Initialise the class with the appropriate values.

        :arg int start: Start position.
        :arg int end: End position.
        :arg int sample_start: Start position.
        :arg int sample_end: End position.
        :arg unicode type: Variant type.
        :arg unicode deleted: Deleted part of the reference sequence.
        :arg ISeqList inserted: Inserted part.
        :arg int shift: Amount of freedom.
        :arg int term: Number of positions until stop codon.
        """
        self.start = start
        self.end = end
        self.sample_start = sample_start
        self.sample_end = sample_end
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.shift = shift
        self.term = term


    def __unicode__(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        Note that this function relies on the absence of values to make the
        correct description. The method used in the DNAVar is better.

        :returns unicode: The HGVS description of the raw variant stored in
            this class.
        """
        if self.type == 'unknown':
            return '?'
        if self.type == 'none':
            return '='

        description = ''
        if not self.deleted:
            if self.type == 'ext':
                description += '*'
            else:
                description += '{}'.format(seq3(self.start_aa))
        else:
            description += '{}'.format(seq3(self.deleted))
        description += '{}'.format(self.start)
        if self.end:
            description += '_{}{}'.format(seq3(self.end_aa), self.end)
        if self.type not in ('subst', 'stop', 'ext', 'fs'): # fs is not a type
            description += self.type
        if self.inserted:
            description += '{}'.format(seq3(self.inserted))

        if self.type == 'stop':
            return description + '*'
        if self.term:
            return description + 'fs*{}'.format(self.term)
        return description
