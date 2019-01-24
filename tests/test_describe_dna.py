import pytest

import extractor


TESTS = [
    # No variants
    ('ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     'ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'location': {'end':   {'position': 44, 'type': 'point'},
                    'start': {'position': 0, 'type': 'point'},
                    'type': 'range'},
                    'type': 'equal'}])
]


@pytest.mark.parametrize('reference,observed,model', TESTS)
def test_variants(reference, observed, model):
    assert model == extractor.describe_dna(reference, observed)
