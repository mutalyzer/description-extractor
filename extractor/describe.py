"""
Generate a HGVS description of the variant(s) leading from one sequence to an
other.
"""


from __future__ import (absolute_import, division, print_function,
    unicode_literals)

import math

from Bio.Seq import reverse_complement

from .variant import (ISeq, AISeq, ISeqList, AISeqList, DNAVar, ProteinVar,
    Allele, ProteinAllele, FS)
from . import extractor, util

from crossmapper import Crossmap


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
    s_revcomp = reverse_complement(str(s)) # FIXME str inserted.

    for i in range(int(math.ceil(len(s) / 2.0))):
        if s[i] != s_revcomp[i]:
            # The first i elements are 'palindromic'.
            return i

    # Perfect 'palindrome'.
    return -1


def var_to_dna_var(s1, s2, var, seq_list=[], weight_position=1):
    """
    Convert a variant from the extractor module to a DNAVar.

    :arg unicode s1: Reference sequence.
    :arg unicode s2: Sample sequence.
    :arg str var: Variant from the extractor module.
    :arg str seq_list: Container for an inserted sequence.
    :arg int weight_position: Weight of a position.
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [DNAVar(type='unknown', weight_position=weight_position)]

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
                '$' not in s1[var.reference_start - ins_length:var.reference_start] and
                '$' not in s2[var.sample_start:var.sample_end] and
                s1[var.reference_start - ins_length:var.reference_start] ==
                s2[var.sample_start:var.sample_end] and
                s2[var.sample_start - ins_length:var.sample_start] ==
                s2[var.sample_start:var.sample_end]):
            # NOTE: We may want to omit the inserted / deleted sequence and
            # use the ranges instead.

            return DNAVar(start=var.reference_start - ins_length + 1,
                end=var.reference_end, type='dup', shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.sample_end],
                    weight_position=weight_position)]),
                weight_position=weight_position)

        return DNAVar(start=var.reference_start,
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

        return DNAVar(start=var.reference_start + 1,
            end=var.reference_end, type='del', shift=shift,
            sample_start=var.sample_start, sample_end=var.sample_end + 1,
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
            var.sample_start + 1 == var.sample_end):
        return DNAVar(start=var.reference_start + 1,
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

        return DNAVar(start=var.reference_start + 1,
            end=var.reference_end, type='inv',
            sample_start=var.sample_start + 1, sample_end=var.sample_end,
            deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            inserted=ISeqList([ISeq(sequence=s2[
                var.sample_start:var.sample_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # InDel.
    return DNAVar(start=var.reference_start + 1,
        end=var.reference_end, deleted=ISeqList([ISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
        inserted=seq_list or
        ISeqList([ISeq(sequence=s2[var.sample_start:var.sample_end],
            weight_position=weight_position)]),
        type='delins', sample_start=var.sample_start + 1,
        sample_end=var.sample_end, weight_position=weight_position)


def var_to_protein_var(s1, s2, var, seq_list=[], weight_position=1):
    """
    Convert a variant from the extractor module to a ProteinVar.

    :arg unicode s1: Reference sequence.
    :arg unicode s2: Sample sequence.
    :arg str var: Variant from the extractor module.
    :arg str seq_list: Container for an inserted sequence.
    :arg str weight_position: Weight of a position.
    """
    # Unknown.
    if s1 == '?' or s2 == '?':
        return [ProteinVar(type='unknown', weight_position=weight_position)]

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
            return ProteinVar(s1=s1, s2=s2,
                start=var.reference_start - ins_length + 1,
                end=var.reference_end, type='dup', shift=shift,
                sample_start=var.sample_start + 1, sample_end=var.sample_end,
                inserted=AISeqList([AISeq(sequence=s2[
                var.sample_start:var.sample_end],
                    weight_position=weight_position)]),
                weight_position=weight_position)

        return ProteinVar(s1=s1, s2=s2, start=var.reference_start,
            end=var.reference_start + 1,
            inserted=seq_list or
            AISeqList([AISeq(sequence=s2[var.sample_start:var.sample_end],
                weight_position=weight_position)]),
            type='ins', shift=shift, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, weight_position=weight_position)

    # Deletion.
    if var.sample_start == var.sample_end:
        shift5, shift3 = roll(s1, var.reference_start + 1, var.reference_end)
        shift = shift5 + shift3

        var.reference_start += shift3
        var.reference_end += shift3

        return ProteinVar(s1=s1, s2=s2, start=var.reference_start + 1,
            end=var.reference_end, type='del', shift=shift,
            sample_start=var.sample_start, sample_end=var.sample_end + 1,
            deleted=AISeqList([AISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # Substitution.
    if (var.reference_start + 1 == var.reference_end and
            var.sample_start + 1 == var.sample_end):
        return ProteinVar(s1=s1, s2=s2, start=var.reference_start + 1,
            end=var.reference_end, sample_start=var.sample_start + 1,
            sample_end=var.sample_end, type='subst',
            deleted=AISeqList([AISeq(sequence=s1[var.reference_start],
                weight_position=weight_position)]),
            inserted=AISeqList([AISeq(sequence=s2[var.sample_start],
                weight_position=weight_position)]),
            weight_position=weight_position)

    # InDel.
    return ProteinVar(s1=s1, s2=s2, start=var.reference_start + 1,
        end=var.reference_end, deleted=AISeqList([AISeq(sequence=s1[
                var.reference_start:var.reference_end],
                weight_position=weight_position)]),
        inserted=seq_list or
        AISeqList([AISeq(sequence=s2[var.sample_start:var.sample_end],
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
            description.append(var_to_dna_var(s1, s2, variant,
                weight_position=extracted.weight_position))

        if variant.type & extractor.TRANSPOSITION_CLOSE:
            in_transposition -= 1

            if not in_transposition:
                description.append(var_to_dna_var(s1, s2, variant, seq_list,
                    weight_position=extracted.weight_position))

    if not description:
        return Allele([DNAVar()])
    return description



def mask_string(string, units):
    MASK = '$'
    repeats = []
    for unit in units:
        found = string.find(unit)
        while found != -1:
            last = len(repeats) - 1
            if last > -1 and repeats[last]['start'] + repeats[last]['count'] * len(unit) == found:
                repeats[last]['count'] += 1
            else:
                repeats.append({'start': found, 'count': 1, 'unit': unit})
            string = string[:found] + MASK * len(unit) + string[found + len(unit):]
            found = string.find(unit, found + len(unit))
    return string, sorted(repeats, key=lambda repeat: repeat['start'])

def describe_repeats(reference, sample, units):
    MASK = '$'
    masked_ref, _ = mask_string(reference, units)
    masked_alt, repeats = mask_string(sample, units)

    reference_start = max(0, masked_ref.find(MASK))
    reference_end = min(len(masked_ref), masked_ref.rfind(MASK)) + 1
    sample_start = max(0, masked_alt.find(MASK))
    sample_end = min(len(masked_alt), masked_alt.rfind(MASK)) + 1

    description = Allele()
    prefix = describe_dna(reference[:reference_start], sample[:sample_start])
    for variant in prefix:
        if variant.type != 'none':
            description.append(variant)

    ref_swig = util.swig_str(masked_ref[reference_start:reference_end])
    alt_swig = util.swig_str(masked_alt[sample_start:sample_end])
    extracted = extractor.extract(ref_swig[0], ref_swig[1],
                                  alt_swig[0], alt_swig[1], extractor.TYPE_DNA)

    variant_list = []
    for variant in extracted.variants:
        replaced = []
        i = variant.sample_start
        start = i
        while i < variant.sample_end:
            while i < variant.sample_end and masked_alt[i + sample_start] != MASK:
                i += 1
            if i < variant.sample_end:
                split = extractor.Variant()
                split.reference_start = variant.reference_start + reference_start
                split.reference_end = variant.reference_end + reference_start
                split.sample_start = start
                split.sample_end = i
                split.type = variant.type
                split.transposition_start = variant.transposition_start + reference_start
                split.transposition_end = variant.transposition_end + reference_start
                replaced.append(split)
                while i < variant.sample_end and masked_alt[i + sample_start] == MASK:
                    i += 1
                start = i

        if len(replaced) > 0:
            split = extractor.Variant()
            split.reference_start = variant.reference_start + reference_start
            split.reference_end = variant.reference_end + reference_start
            split.sample_start = start
            split.sample_end = variant.sample_end
            split.type = variant.type
            split.transposition_start = variant.transposition_start + reference_start
            split.transposition_end = variant.transposition_end + reference_start
            replaced.append(split)
            variant_list += replaced
        else:
            variant_list.append(variant)

    in_transposition = 0
    index = 0
    repeat = 0
    seq_list = ISeqList()
    for variant in variant_list:
        while variant.sample_start > index and repeat < len(repeats):
            seq_list.append(DNAVar(type='repeat',inserted=repeats[repeat]['unit'],count=repeats[repeat]['count']))
            index = repeats[repeat]['start'] + repeats[repeat]['count'] * len(repeats[repeat]['unit'])
            repeat += 1

        if variant.type & extractor.TRANSPOSITION_OPEN:
            in_transposition += 1

        if in_transposition:
            if variant.type & extractor.IDENTITY:
                seq_list.append(ISeq(start=variant.transposition_start + 1 + reference_start,
                    end=variant.transposition_end + reference_start, reverse=False,
                    weight_position=extracted.weight_position))
            elif variant.type & extractor.REVERSE_COMPLEMENT:
                seq_list.append(ISeq(start=variant.transposition_start + 1 + reference_start,
                    end=variant.transposition_end + reference_start, reverse=True,
                    weight_position=extracted.weight_position))
            else: #bases insertion
                seq_list.append(ISeq(
                    sequence=sample[variant.sample_start + sample_start:variant.sample_end + sample_start],
                    weight_position=extracted.weight_position))
        elif variant.type & extractor.IDENTITY:
            seq_list.append(ISeq(start=variant.reference_start + 1 + reference_start,
                end=variant.reference_end + reference_start, reverse=False,weight_position=extracted.weight_position))
        elif variant.type & extractor.REVERSE_COMPLEMENT:
            seq_list.append(ISeq(start=variant.reference_start + 1 + reference_start,
                end=variant.reference_end + reference_start, reverse=True,weight_position=extracted.weight_position))
        elif variant.sample_end != variant.sample_start:
            seq_list.append(ISeq(sequence=sample[variant.sample_start + sample_start:variant.sample_end + sample_start],
                weight_position=extracted.weight_position))

        if variant.type & extractor.TRANSPOSITION_CLOSE:
            in_transposition -= 1

        index = variant.sample_end

    while repeat < len(repeats):
        seq_list.append(DNAVar(type='repeat',inserted=repeats[repeat]['unit'],count=repeats[repeat]['count']))
        repeat += 1

    if len(variant_list) > 0 or len(repeats) > 0:
        description.append(DNAVar(start=reference_start + 1,end=reference_end,sample_start=sample_start,sample_end=sample_end,type='struct',inserted=seq_list))

        suffix = describe_dna(reference[reference_end:], sample[sample_end:])
        for variant in suffix:
            if variant.type != 'none':
                variant.start += reference_end
                variant.end += reference_end
                description.append(variant)
    else:
        description = prefix

    cm = Crossmap([reference_start + 1, reference_end], [], 1)
    for variant in description:
        for inserted in variant.inserted:
            inserted.start = cm.tuple2string(cm.g2x(int(inserted.start)))
            inserted.end = cm.tuple2string(cm.g2x(int(inserted.end)))
        variant.start = cm.tuple2string(cm.g2x(int(variant.start)))
        variant.end = cm.tuple2string(cm.g2x(int(variant.end)))

    return description, reference_start, reference_end


def print_var(variant):
    print('({:3}, {:3}), ({:3}, {:3}), {:08b}, {}, {}'.format(variant.reference_start,
        variant.reference_end, variant.sample_start, variant.sample_end,
        variant.type, variant.type, variant.sample_end - variant.sample_start))


def get_frames(flags):
    result = []

    for fs in FS:
        if flags & FS[fs]:
            result.append(fs)

    return result


def describe_protein(s1, s2, codon_table=1):
    """
    """
    codons = util.codon_table_string(codon_table)

    description = ProteinAllele()

    s1_swig = util.swig_str(s1)
    s2_swig = util.swig_str(s2)
    codons_swig = util.swig_str(codons)

    extracted = extractor.extract(s1_swig[0], s1_swig[1],
        s2_swig[0], s2_swig[1], extractor.TYPE_PROTEIN, codons_swig[0])
    variants = extracted.variants

    #for variant in variants:
    #    print_var(variant)
    #print()

    index = 0
    while index < len(variants):
        if variants[index].type != extractor.IDENTITY:
            variant = variants[index]
            index += 1
            seq_list = AISeqList()

            # NOTE: This is for filling.
            last_end = variants[index].reference_start

            while (index < len(variants) and
                    variants[index].type & extractor.FRAME_SHIFT):

                if last_end < variants[index].sample_start:
                    seq_list.append(AISeq(
                        s2[last_end:variants[index].sample_start]))

                last_end = variants[index].sample_end

                seq_list.append(AISeq(
                    s2[variants[index].sample_start:
                        variants[index].sample_end],
                    start=variants[index].reference_start + 1,
                    end=variants[index].reference_end,
                    sample_start=variants[index].sample_start + 1,
                    sample_end=variants[index].sample_end,
                    frames=get_frames(variants[index].type)))

                # NOTE: Perhaps use trans_open, trans_close to ...
                index += 1

            if last_end < variant.sample_end:
                seq_list.append(AISeq(s2[last_end:variant.sample_end]))

            var = var_to_protein_var(s1, s2, variant, seq_list,
                weight_position=extracted.weight_position)
            description.append(var)
        else:
            index += 1

    if not description:
        return ProteinAllele([ProteinVar()])
    return description
