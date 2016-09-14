"""
Models for the description extractor.
"""


from __future__ import (
    absolute_import, division, print_function, unicode_literals)

from . import extractor
from extractor.util import python_2_unicode_compatible, seq3, str


WEIGHTS = {
    'subst': extractor.WEIGHT_SUBSTITUTION,
    'del': extractor.WEIGHT_DELETION,
    'ins': extractor.WEIGHT_INSERTION,
    'dup': extractor.WEIGHT_INSERTION,
    'inv': extractor.WEIGHT_INVERSION,
    'delins': extractor.WEIGHT_DELETION_INSERTION
}
FS = {
    '1': extractor.FRAME_SHIFT_1,
    '2': extractor.FRAME_SHIFT_2,
    'inv': extractor.FRAME_SHIFT_REVERSE,
    '1inv': extractor.FRAME_SHIFT_REVERSE_1,
    '2inv': extractor.FRAME_SHIFT_REVERSE_2
}


@python_2_unicode_compatible
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


    def __str__(self):
        if len(self.items) > 1:
            return '[{0}]'.format(';'.join(map(str, self.items)))
        return str(self.items[0])


    def append(self, item):
        self.items.append(item)


    def weight(self):
        weight = sum(map(lambda x: x.weight(), self.items))

        if len(self.items) > 1:
            return weight + (len(self.items) + 1) * extractor.WEIGHT_SEPARATOR
        return weight


class Allele(HGVSList):
    pass


class ProteinAllele(HGVSList):
    def nhgvs(self):
        if len(self.items) > 1:
            return '[{0}]'.format(';'.join(map(lambda x: x.nhgvs(),
                self.items)))
        return self.items[0].nhgvs()


class ISeqList(HGVSList):
    pass


class AISeqList(ISeqList):
    def get_sequence(self):
        return ''.join(map(lambda x: x.sequence, self.items))


@python_2_unicode_compatible
class ISeq(object):
    """
    Container for an inserted sequence.
    """
    def __init__(
            self, sequence='', start=0, end=0, reverse=False,
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


    def __str__(self):
        if self.type == 'ins':
            return self.sequence

        if not (self.start or self.end):
            return ''

        if self.start == self.end:
             return '{}'.format(self.start)

        inverted = 'inv' if self.reverse else ''
        return '{0}_{1}{2}'.format(self.start, self.end, inverted)


    # TODO: Is this still used?
    def __bool__(self):
         return bool(self.sequence)


    # TODO: Is this still used?
    def __nonzero__(self): # Python 2.x compatibility.
        return self.__bool__()


    def weight(self):
        if self.type == 'ins':
            return len(self.sequence) * extractor.WEIGHT_BASE

        inverse_weight = WEIGHTS['inv'] if self.reverse else 0
        return (
            self.weight_position * 2 + extractor.WEIGHT_SEPARATOR +
            inverse_weight)


@python_2_unicode_compatible
class AISeq(object):
    """
    Container for an annotated inserted sequence.
    """
    def __init__(
            self, sequence='', start=0, end=0, sample_start=0, sample_end=0,
            frames=[], weight_position=1):
        """
        Initialise the class with the appropriate values.

        :arg unicode sequence: Literal inserted sequence.
        :arg int start: Start position for a transposed sequence.
        :arg int end: End position for a transposed sequence.
        """
        self.sequence = sequence
        self.start = start
        self.end = end
        self.sample_start = sample_start
        self.sample_end = sample_end
        self.weight_position = weight_position
        self.frames = frames

        self.type = 'trans'
        if self.sequence:
            self.type = 'ins'
        if self.frames:
            self.type = 'fs'


    def __str__(self):
        if self.type == 'ins':
            return self.sequence

        if self.type == 'trans':
            if self.start != self.end:
                return '{}_{}'.format(self.start, self.end)
            return '{}'.format(self.start)

        return '{}_{}{}|{}'.format(
            self.start, self.end, self.sequence, '|'.join(self.frames))


@python_2_unicode_compatible
class DNAVar(object):
    """
    Container for a DNA variant.
    """
    def __init__(
            self, start=0, start_offset=0, end=0, end_offset=0, sample_start=0,
            sample_start_offset=0, sample_end=0, sample_end_offset=0, count=0,
            type='none', deleted=ISeqList([ISeq()]),
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
        :arg int count
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
        self.count = count
        self.type = type
        self.deleted = deleted
        self.inserted = inserted
        self.weight_position = weight_position
        self.shift = shift


    def __str__(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        :returns unicode: The HGVS description of the raw variant stored in
            this class.
        """
        if self.type == 'unknown':
            return '?'
        if self.type == 'none':
            return '='

        if self.count > 0:
            return '{0}({1})'.format(self.inserted, self.count)

        description = str(self.start)

        if self.start != self.end:
            description += '_{0}'.format(self.end)

        if self.type != 'subst':
            description += str(self.type)

            if self.type in ('ins', 'delins'):
                return description + str(self.inserted)
            elif self.type == 'struct':
                return str(self.inserted)
            return description

        return description + '{0}>{1}'.format(self.deleted, self.inserted)


    def weight(self):
        if self.type == 'unknown':
            return -1
        if self.type == 'none':
            return 0

        weight = self.weight_position
        if self.start != self.end:
            weight += self.weight_position + extractor.WEIGHT_SEPARATOR

        return weight + WEIGHTS[self.type] + self.inserted.weight()


@python_2_unicode_compatible
class ProteinVar(object):
    """
    Container for a protein variant.

    """
    def __init__(
            self, s1='', s2='', start=0, end=0, sample_start=0, sample_end=0,
            type='none', deleted=ISeqList([ISeq()]),
            inserted=AISeqList([AISeq()]), shift=0, term=0, weight_position=1):
        """
        Initialise the class with the appropriate values.

        :arg unicode s1: Reference sequence.
        :arg unicode s2: Sample sequence.
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
        self.start_aa = s1[start - 1]
        self.end_aa = s1[end - 1]
        self.sample_start_aa = s2[sample_start - 1]
        self.sample_end_aa = s2[sample_end - 1]
        self.type = type
        self.deleted = deleted
        self.inserted = ISeqList([ISeq(inserted.get_sequence())])
        self.annotated_inserted = inserted
        self.shift = shift
        self.term = term


    def __str__(self):
        """
        Give the HGVS description of the raw variant stored in this class.

        :returns unicode: The HGVS description of the raw variant stored in
            this class.
        """
        # TODO: ext*
        if self.type == 'unknown':
            return '?'
        if self.type == 'none':
            return '='

        description = '{}{}'.format(seq3(self.start_aa), self.start)
        if self.term:
            return description + '{}fs*{}'.format(
                seq3(self.inserted[0].sequence[0]), self.term)
        if self.start != self.end:
            description += '_{}{}'.format(seq3(self.end_aa), self.end)

        if self.type != 'subst':
            description += self.type

            if self.type in ('ins', 'delins'):
                return description + seq3(str(self.inserted)) # FIXME: str
            return description
        return description + seq3(self.inserted)


    def nhgvs(self):
        """
        """
        if self.type == 'unknown':
            return '?'
        if self.type == 'none':
            return '='

        description = str(self.start)
        if self.start != self.end:
            description += '_{}'.format(self.end)

        if self.type != 'subst':
            description += self.type

            if self.type in ('ins', 'delins'):
                return description + str(self.annotated_inserted)
            return description
        return description + '{}>{}'.format(
            self.deleted, self.annotated_inserted)


@python_2_unicode_compatible
class FrameShiftAnnotation(object):
    """
    Container for frame shift annotation.
    """
    def __init__(self, s2='', start=0, end=0, sample_start=0, sample_end=0,
            type='none'):
        """
        Initialise the class with the appropriate values.

        :arg unicode s2: Sample sequence.
        :arg int start: Start position.
        :arg int end: End position.
        :arg int sample_start: Start position.
        :arg int sample_end: End position.
        :arg unicode type: Frame shift type.
        """
        self.start = start
        self.end = end
        self.sample_start = sample_start
        self.sample_end = sample_end
        self.seq = s2[sample_start - 1:sample_end]
        for fs_type in FS:
            if FS[fs_type] & type:
                self.type = fs_type


    def __str__(self):
        """
        """
        return '{}_{}{}|{}'.format(self.start, self.end, self.seq, self.type)
