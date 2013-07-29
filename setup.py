import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('extractor requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = []

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import extractor as distmeta

setup(
    name='extractor',
    version=distmeta.__version__,
    description='FASTA/FASTQ analysis and manipulation toolkit.',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['extractor'],
    install_requires=requires,
    entry_points = {
        'console_scripts': [
        ]
    },
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
