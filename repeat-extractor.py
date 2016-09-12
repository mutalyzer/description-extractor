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
min_length = 2

with open('strlist.txt', 'r') as infile:
    lines = infile.readlines()

sequences = {}

for line in lines:
    label, string = line.split('\t')
    if label in sequences:
        sequences[label].append(string.strip())
    else:
        sequences[label] = [string.strip()]

for sequence in sequences:
    reference = sequences[sequence][0]
    repeats = short_sequence_repeat_extractor(reference, min_length)
    units = {}
    for repeat in repeats:
        if repeat.count + 1 >= min_count:
            units[reference[repeat.start:repeat.end]] = repeat.count + 1
    unit_list = []
    for unit in units:
        unit_list.append(unit)

    print sequence, unit_list
    for string in sequences[sequence]:
        describe_repeats(reference, string, unit_list)
    print


