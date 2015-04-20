# Extractor

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

## Installation

### C++ Library only

Run `make`.

Optionally set the `__debug__` flag to trace the algorithm.

For direct use within a C/C++ environment just
`#include "extractor.h"` and add `extractor.cc` to your project's
source files.

### Within Mutalyzer

Todo

## Examples

```
Variant_List extract(char_t const* const reference,
                     size_t const        reference_length,
                     char_t const* const sample,
                     size_t const        sample_length,
                     int const           type = TYPE_DNA,
                     char_t const* const codon_string = 0)
```

## Testing

There are some unit tests for the Python interface. After installing the
Python package, run them using [http://pytest.org/](pytest):

    pip install pytest
    python setup.py develop
    py.test

Alternatively, use [https://tox.readthedocs.org/](tox) to automatically run
the tests on all supported versions of Python:

    pip install tox
    tox

