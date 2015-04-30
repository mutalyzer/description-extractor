"""
Generate a HGVS description of the variant(s) leading from one sequence to an
other.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math

from .variant import ISeq, ISeqList, DNAVar, ProteinVar, Allele
from . import extractor, util


# Taken from BioPython.
AMBIGUOUS_DNA_COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'M': 'K',
    'R': 'Y',
    'W': 'W',
    'S': 'S',
    'Y': 'R',
    'K': 'M',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'B': 'V',
    'X': 'X',
    'N': 'N'}
AMBIGUOUS_RNA_COMPLEMENT = {
    'A': 'U',
    'C': 'G',
    'G': 'C',
    'U': 'A',
    'M': 'K',
    'R': 'Y',
    'W': 'W',
    'S': 'S',
    'Y': 'R',
    'K': 'M',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'B': 'V',
    'X': 'X',
    'N': 'N'}


def _make_translation_table(complement_mapping):
    before = list(complement_mapping.keys())
    before += [b.lower() for b in before]
    after = list(complement_mapping.values())
    after += [b.lower() for b in after]
    return dict((ord(k), v) for k, v in zip(before, after))


_dna_complement_table = _make_translation_table(AMBIGUOUS_DNA_COMPLEMENT)
_rna_complement_table = _make_translation_table(AMBIGUOUS_RNA_COMPLEMENT)


def reverse_complement(sequence):
    """
    Reverse complement of a sequence represented as unicode string.
    """
    if 'U' in sequence or 'u' in sequence:
        table = _rna_complement_table
    else:
        table = _dna_complement_table

    return ''.join(reversed(sequence.translate(table)))


def roll(s, first, last):
    """
    Determine the variability of a variant by looking at cyclic
    permutations. Not all cyclic permutations are tested at each time, it
    is sufficient to check ``aW'' if ``Wa'' matches (with ``a'' a letter,
    ``W'' a word) when rolling to the left for example.

        >>> roll('abbabbabbabb', 4, 6)
        (3, 6)
        >>> roll('abbabbabbabb', 5, 5)
        (0, 1)
        >>> roll('abcccccde', 4, 4)
        (1, 3)

    @arg s: A reference sequence.
    @type s: any sequence type
    @arg first: First position of the pattern in the reference sequence.
    @type first: int
    @arg last: Last position of the pattern in the reference sequence.
    @type last: int

    @return: tuple:
        - left  ; Amount of positions that the pattern can be shifted to
                  the left.
        - right ; Amount of positions that the pattern can be shifted to
                  the right.
    @rtype: tuple(int, int)
    """
    pattern = s[first - 1:last]   # Extract the pattern
    pattern_length = len(pattern)

    # Keep rolling to the left as long as a cyclic permutation matches.
    minimum = first - 2
    j = pattern_length - 1
    while minimum > -1 and s[minimum] == pattern[j % pattern_length]:
        j -= 1
        minimum -= 1

    # Keep rolling to the right as long as a cyclic permutation matches.
    maximum = last
    j = 0
    while maximum < len(s) and s[maximum] == pattern[j % pattern_length]:
        j += 1
        maximum += 1

    return first - minimum - 2, maximum - last


def palinsnoop(s):
    """
    Check a sequence for a reverse-complement-palindromic prefix (and
    suffix). If one is detected, return the length of this prefix. If the
    string equals its reverse complement, return -1.

        >>> palinsnoop('TACGCTA')
        2
        >>> palinsnoop('TACGTA')
        -1
        >>> palinsnoop('TACGCTT')
        0

    @arg s: A nucleotide sequence.
    @type s: unicode

    @return: The number of elements that are palindromic or -1 if the string
             is a 'palindrome'.
    @rtype: int
    """
    s_revcomp = reverse_complement(s)

    for i in range(int(math.ceil(len(s) / 2.0))):
        if s[i] != s_revcomp[i]:
            # The first i elements are 'palindromic'.
            return i

    # Perfect 'palindrome'.
    return -1


def var_to_rawvar(s1, s2, var, seq_list=[], container=DNAVar,
        weight_position=1):
    """
    Convert a variant from the extractor module to one of the RawVar
    subclasses.

    :arg unicode s1: Reference sequence.
    :arg unicode s2: Sample sequence.
    :arg str var: Variant from the extractor module.
    :arg str seq_list: Container for an inserted sequence.
    :arg str container: Destination container.
    :arg str weight_position: Weight of a position.
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [container(type='unknown', weight_position=weight_position)]

    # Insertion / Duplication.
    if var.reference_start == var.reference_end:
        ins_length = var.sample_end - var.sample_start
        shift5, shift3 = roll(s2, var.sample_start + 1, var.sample_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3
        var.sample_start += shift3
        var.sample_end += shift3

        if (var.sample_start - ins_length >= 0 and
                s1[var.reference_start - ins_length:var.reference_start] ==
                s2[var.sample_start:var.sample_end]):
            # NOTE: We may want to omit the inserted / deleted sequence and
            # use the ranges instead.
            return container(start=var.reference_start - ins_length + 1,
                end=var.reference_end, type='dup', shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.sample_end],
                    weight_position=weight_position)]),
                weight_position=weight_position)

        return container(start=var.reference_start,
            end=var.reference_start + 1,
            inserted=seq_list or
            ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
                weight_position=weight_position)]),
            type='ins', shift=shift, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, weight_position=weight_position)

    # Deletion.
    if var.sample_start == var.sample_end:
        shift5, shift3 = roll(s1, var.reference_start + 1, var.reference_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3

        return container(start=var.reference_start + 1,
            end=var.reference_end, type='del', shift=shift,
            sample_start=var.sample_start, sample_end=var.sample_end + 1,
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
            var.sample_start + 1 == var.sample_end):
        return container(start=var.reference_start + 1,
            end=var.reference_end, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, type='subst',
            deleted=ISeqList([ISeq(sequence=s1[var.reference_start],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[var.sample_start],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # Inversion.
    if var.type & extractor.REVERSE_COMPLEMENT:
        trim = palinsnoop(s1[var.reference_start:var.reference_end])

        if trim > 0: # Partial palindrome.
            var.reference_end -= trim
            var.sample_end -= trim

        return container(start=var.reference_start + 1,
            end=var.reference_end, type='inv',
            sample_start=var.sample_start + 1, sample_end=var.sample_end,
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # InDel.
    return container(start=var.reference_start + 1,
        end=var.reference_end, deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
        inserted=seq_list or
        ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
            weight_position=weight_position)]),
        type='delins', sample_start=var.sample_start + 1,
        sample_end=var.sample_end, weight_position=weight_position)


def describe_dna(s1, s2):
    """
    Give an allele description of the change from {s1} to {s2}.

    :arg unicode s1: Sequence 1.
    :arg unicode s2: Sequence 2.

    :returns list(RawVar): A list of RawVar objects, representing the allele.
    """
    description = Allele()
    in_transposition = 0

    s1_swig = util.swig_str(s1)
    s2_swig = util.swig_str(s2)
    extracted = extractor.extract(s1_swig[0], s1_swig[1],
                                  s2_swig[0], s2_swig[1], extractor.TYPE_DNA)

    for variant in extracted.variants:
       # print(variant.type, variant.reference_start,
       #     variant.reference_end, variant.sample_start,
       #     variant.sample_end, variant.transposition_start,
       #     variant.transposition_end)
       # print(variant.type & extractor.TRANSPOSITION_OPEN, variant.type &
       #     extractor.TRANSPOSITION_CLOSE)

        if variant.type & extractor.TRANSPOSITION_OPEN:
            if not in_transposition:
                seq_list = ISeqList()
            in_transposition += 1

        if in_transposition:
            if variant.type & extractor.IDENTITY:
                seq_list.append(ISeq(start=variant.transposition_start + 1,
                    end=variant.transposition_end, reverse=False,
                    weight_position=extracted.weight_position))
            elif variant.type & extractor.REVERSE_COMPLEMENT:
                seq_list.append(ISeq(start=variant.transposition_start + 1,
                    end=variant.transposition_end, reverse=True,
                    weight_position=extracted.weight_position))
            else:
                seq_list.append(ISeq(
                    sequence=s2[variant.sample_start:variant.sample_end],
                    weight_position=extracted.weight_position))
        elif not (variant.type & extractor.IDENTITY):
            description.append(var_to_rawvar(s1, s2, variant,
                weight_position=extracted.weight_position))

        if variant.type & extractor.TRANSPOSITION_CLOSE:
            in_transposition -= 1

            if not in_transposition:
                description.append(var_to_rawvar(s1, s2, variant, seq_list,
                    weight_position=extracted.weight_position))

    if not description:
        return Allele([DNAVar()])
    return description


def describe_protein(s1, s2):
    # Placeholder for the interface of the protein extraction.
    pass
