import pytest

import extractor


REFERENCE = 'ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT'


TESTS = [
    # No variants
    (REFERENCE,
     REFERENCE,
     [{'location': {'end':   {'position': 44, 'type': 'point'},
                    'start': {'position': 0, 'type': 'point'},
                    'type': 'range'},
                    'type': 'equal'}]),

    # Single substitution: 7A>G
    (REFERENCE,
     'ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single deletion: 7del
    (REFERENCE,
     'ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single deletion multiple bases: 7_8del
    (REFERENCE,
     'ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single insertion: 6_7insC
    (REFERENCE,
     'ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single insertion multiple bases: 6_7insCC
    (REFERENCE,
     'ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single inversion: 7_11inv
    (REFERENCE,
     'ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single delins: 7delinsCC
    (REFERENCE,
     'ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single delins: 21_23delinsTTTT
    (REFERENCE,
     'ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single inversion: 6_7inv
    (REFERENCE,
     'ACGTCTCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

    # Single inversion: 7_8delinsTC
    (REFERENCE,
     'ACGTCGTCTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{}]),

]


@pytest.mark.parametrize('reference, observed, variants', TESTS)
def test_variants(reference, observed, variants):
    assert variants == extractor.describe_dna(reference, observed)
