from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess

class git_clone_external(build_ext):
    def run(self):
        subprocess.check_call(['rm', '-rf', 'extractor-core'])
        subprocess.check_call(['git', 'clone', 'git@github.com:mutalyzer/extractor-core.git'])
        build_ext.run(self)

extractor = Extension('extractor', sources = ['extractor-module.cc',
                                              'extractor-core/src/extractor.cc'])
setup(name = 'description-extractor',
      version = '3.0.0',
      cmdclass = {'build_ext': git_clone_external},
      description = 'HGVS variant description extractor',
      ext_modules = [extractor]
)
