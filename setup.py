import os
import sys
from setuptools import setup
from distutils.core import Extension

# We do some trickery to assure SWIG is always run before installing the
# generated files.
# http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module
from setuptools.command.install import install
from distutils.command.build import build
class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)
class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        self.do_egg_install()
custom_cmdclass = {'build': CustomBuild, 'install': CustomInstall}

if sys.version_info < (2, 6):
    raise Exception('extractor requires Python 2.6 or higher.')

# This is quite the hack, but we don't want to import our package from here
# since that's recipe for disaster (it might have some uninstalled
# dependencies, or we might import another already installed version).
distmeta = {}
for line in open(os.path.join('extractor', '__init__.py')):
    try:
        field, value = (x.strip() for x in line.split('='))
    except ValueError:
        continue
    value = value.strip('\'"')
    distmeta[field] = value

# The __version__ value is actually defined in extractor.h.
for line in open(os.path.join('extractor', 'extractor.h')):
    if ' VERSION = ' in line:
        version = line.split('=')[-1].replace(';', '').replace('"', '') \
                                     .replace("'", '').strip()
        distmeta['__version__'] = version
        distmeta['__version_info__'] = tuple(version.split('.'))
        break

try:
    with open('README.md') as readme:
        long_description = readme.read()
except IOError:
    long_description = 'See ' + distmeta['__homepage__']

setup(
    name='description-extractor',
    cmdclass=custom_cmdclass,
    ext_modules=[Extension('_extractor', ['extractor/extractor.i',
        'extractor/extractor.cc'], swig_opts=['-c++'])],
    version=distmeta['__version__'],
    description='HGVS variant description extractor',
    long_description=long_description,
    author=distmeta['__author__'],
    author_email=distmeta['__contact__'],
    url=distmeta['__homepage__'],
    license='MIT License',
    platforms=['any'],
    packages=['extractor'],
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics',
    install_requires=['biopython']
)
