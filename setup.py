# -*- coding: utf-8 -*-
import sys, os, re
from setuptools import setup
from setuptools.command.test import test as TestCommand

def find_version():
    filepath = os.path.join('gseapy', '__main__.py')
    with open(os.path.join(os.path.dirname(__file__), filepath), encoding="utf8") as fp:
        content = fp.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", content, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

__author__ = 'Zhuoqing Fang'
__version__ = find_version()


if sys.argv[-1] == 'publish':
    os.system("python setup.py sdist bdist_wheel register upload")
    print("You probably want to also tag the version now:")
    print("  git tag -a %s -m 'version %s'" % (__version__,__version__))
    print("  git push --tags")
    sys.exit()

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='gseapy',
      version=__version__,
      description='Gene Set Enrichment Analysis in Python',
      long_description=readme(),
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries'],
      keywords= ['Gene Ontology', 'GO','Biology', 'Enrichment',
          'Bioinformatics', 'Computational Biology',],
      url='https://github.com/zqfang/gseapy',
      author='Zhuoqing Fang',
      author_email='fzq518@gmail.com',
      license='MIT',
      packages=['gseapy'],
      package_data={'gseapy': ["data/*.txt"],},
      include_package_data=True,
      install_requires=[
                        'numpy>=1.13.0',
                        'scipy',
                        'pandas',
                        'matplotlib',
                        'bioservices',
                        'requests',
                        'joblib'],
      entry_points={'console_scripts': ['gseapy = gseapy.__main__:main'],},
      tests_require=['pytest'],
      cmdclass = {'test': PyTest},
      zip_safe=False,
      download_url='https://github.com/zqfang/gseapy',)


