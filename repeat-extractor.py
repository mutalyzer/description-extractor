#!/usr/bin/env python

from __future__ import unicode_literals

from extractor import *

#ref = 'AGCTGTGGGAGGGAGCCAGTGGATTTGGAAACAGAAATGGCTTGGCCTTGCCTGCCTGCCTGCCTGCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCTGCAATCCTTTAACTTACTGAATAACTCATGATTATGGGCCACCTGCAGGTACCATGCTAG'
#alt = 'AGCTGTGGGAGGGAGCCAGTGGATTTGGAAACAGAAATGGCTTCGCCTTGCCTGCCTGCCTGCCTGCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCGTCCTTCCTTCCCTCCTGCAATCCTATAACTTACTGAATAACTCATGATTATGGGCCACCTGCAGGTACCATGCTAG'
#units = ['TCCT', 'GCCT']

#print 'version:', extractor.VERSION

#describe_repeats(ref, alt, units)

THRESHOLD = 10000


class Repeat(object):
    """
    Simple repeat structure.
    """
    __slots__ = ['start', 'end', 'count']

    def __init__(self, start, end, count=0):
        """
        Initialize the repeat structure.

        :arg integer start: Start position of the repeat.
        :arg integer end: End position of the repeat.
        :arg integer count: Number of repeats.
        """
        self.start = start
        self.end = end
        self.count = count


def short_sequence_repeat_extractor(string, min_length=1):
    """
    Extract the short tandem repeat structure from a string.

    :arg string string: The string.
    :arg integer min_length: Minimum length of the repeat structure.
    """
    length = len(string)

    k_max = length // 2 + 1
    if k_max > THRESHOLD:
        k_max = THRESHOLD // 2

    repeats = []

    i = 0
    last_repeat = i
    while i < length:
        max_count = 0
        max_k = 1
        for k in range(min_length, k_max):
            count = 0
            for j in range(i + k, length - k + 1, k):
                if string[i:i + k] != string[j:j + k]:
                    break
                count += 1

            if count > 0 and count >= max_count:
                max_count = count
                max_k = k

        if max_count > 0:
            if last_repeat < i:
                repeats.append(Repeat(last_repeat, i))
            repeats.append(Repeat(i, i + max_k, max_count))
            last_repeat = i + max_k * (max_count + 1)

        i += max_k * (max_count + 1)

    if last_repeat < i:
        repeats.append(Repeat(last_repeat, i))

    return repeats


min_count = 3
min_length = 3

with open('strlist.txt', 'r') as infile:
    lines = infile.readlines()

sequences = {}

for line in lines:
    label, string = line.split('\t')
    if label in sequences:
        sequences[label].append(string.strip())
    else:
        sequences[label] = [string.strip()]

#select = 'D13S317'
#unit_list = ['TATC']
#reference = sequences[select][0]
#sample = sequences[select][0]
#description, _, _ = describe_repeats(reference, sample, unit_list)
#print 'l.{}'.format(description)


for sequence in sequences:
    best = 0
    for string in sequences[sequence]:
        repeats = short_sequence_repeat_extractor(string, min_length)
        score = 0
        for repeat in repeats:
            if repeat.count + 1 >= min_count:
                score += (repeat.end - repeat.start) * (repeat.end - repeat.start) * (repeat.count + 1)
        if score > best:
            reference = string
            best = score

    repeats = short_sequence_repeat_extractor(reference, min_length)
    units = {}
    for repeat in repeats:
        if repeat.count + 1 >= min_count:
            units[reference[repeat.start:repeat.end]] = repeat.count + 1
    unit_list = []
    for unit in units:
        unit_list.append(unit)

    reference = sequences[sequence][0]
    print sequence + ':',
    print reference
    if best > 0:
        print 'repeat units:', unit_list
    else:
        print 'repeat units: []'
    for string in sequences[sequence]:
        rep_start = 1
        rep_end = len(reference)
        if best > 0:
            description, rep_start, rep_end = describe_repeats(reference, string, unit_list)
        else:
            description = describe_dna(reference, string)
        print '{}({}_{}):l.{}'.format(sequence, rep_start, rep_end, description)
    print


