# description-extractor
HGVS variant description extractor  
Extract a list of differences between two sequences.


## Installation

`git clone git@github.com:mutalyzer/description-extractor.git`
`git checkout mutalyzer3`

`python setup.py install` or `pip install .`

## Use

`import extractor`

`variants = extractor.describe_dna('AAATAA', 'AAAGAAA')`
