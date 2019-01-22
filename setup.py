from setuptools import setup, Extension

from setuptools.command.build_ext import build_ext
import subprocess

class git_clone_external(build_ext):
    def run(self):
        subprocess.check_call(['git', 'clone', 'https://github.com/mutalyzer/extractor-core.git'])
        build_ext.run(self)


extractor = Extension('description-extractor', sources = ['extractor-module.cc',
                                                          'extractor-core/src/extractor.cc'])

setup(name = 'description-extractor',
      version = '3.0',
      cmdclass = {'build_ext': git_clone_external},
      description = 'This is the extractor package',
      ext_modules = [extractor]
)
