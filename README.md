# HGVS variant description extractor

Unambiguous sequence variant descriptions are important in reporting
the outcome of clinical diagnostic DNA tests. The standard
nomenclature of the Human Genome Variation Society (HGVS) describes
the observed variant sequence relative to a given reference sequence.
We propose an efficient algorithm for the extraction of
HGVS descriptions from two sequences with three main requirements in
mind: minimizing the length of the resulting descriptions, minimizing
the computation time, and keeping the unambiguous descriptions
biologically meaningful.

This algorithm is able to compute the HGVS descriptions of complete
chromosomes or other large DNA strings in a reasonable amount of
computation time and its resulting descriptions are relatively small.
Additional applications include updating of gene variant database
contents and reference sequence liftovers.

    >>> from extractor import describe_dna
    >>> print describe_dna('TAACAATGGAAC', 'TAAACAATTGAA')
    [3dup;8G>T;12del]


## Implementation

The core algorithm is implemented in C++ with a Python wrapper providing a
developer friendly interface.


## Installation

### Python package

You need [SWIG](http://www.swig.org/) installed. Then:

    pip install description-extractor


### C++ library only

Run `make`.

Optionally set the `__debug__` flag to trace the algorithm.

For direct use within a C/C++ environment just
`#include "extractor.h"` and add `extractor.cc` to your project's
source files.


## Testing

There are some unit tests for the Python interface. After installing the
Python package, run them using [pytest](http://pytest.org/):

    pip install pytest
    python setup.py develop
    py.test

Alternatively, use [tox](https://tox.readthedocs.org/) to automatically run
the tests on all supported versions of Python:

    pip install tox
    tox
