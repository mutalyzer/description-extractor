# Patch for swig.
from distutils.command.build import build as _build
_build.sub_commands = [
    ('build_ext', _build.has_ext_modules),
    ('build_py', _build.has_pure_modules),
    ('build_clib', _build.has_c_libraries),
    ('build_scripts', _build.has_scripts)
] 

#from setuptools.command.bdist_egg import bdist_egg
#old_run = bdist_egg.run
#
#def run(self):
#    old_run(self)
#    self.run_command("install_lib")
#
#bdist_egg.run = run

import sys
from setuptools import setup
from distutils.core import Extension

if sys.version_info < (2, 6):
    raise Exception('extractor requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = []

import extractor as distmeta

setup(
    name='extractor',
    ext_modules=[Extension('_extractor', ['extractor/extractor.i',
        'extractor/extractor.cc'], swig_opts=['-c++'])],
    py_modules=['extractor.extractor'],
    version=distmeta.__version__,
    description=distmeta.usage[0],
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
