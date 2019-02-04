# description-extractor
HGVS variant description extractor  
Extract a list of differences between two sequences.


## Installation

`git clone git@github.com:mutalyzer/description-extractor.git`

`git checkout mutalyzer3`

`python setup.py install` or `pip install .`

## Testing

`pip install -r requirements-dev` or `pip install pytest`

`py.test` or `py.test tests/`

## Use

`import extractor`

`variants = extractor.describe_dna('TAACAATGGAAC', 'TAAACAATTGAA')`
