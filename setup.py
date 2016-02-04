# -*- coding: utf-8 -*-

from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='gseapy',
      version='0.3.1',
      description='Gene Set Enrichment Analysis in Python',
      long_description=readme(),
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 2.7',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries'],
      keywords= ['GO', 'Gene Ontology', 'Biology', 'Enrichment',
          'Bioinformatics', 'Computational Biology',],
      url='https://github.com/BioNinja/gseapy',
      author='Zhuoqing Fang',
      author_email='fangzhuoqing@sibs.ac.cn',
      license='MIT',
      packages=['gseapy'],
      
      install_requires=[
          'numpy>=1.8.0',
          'pandas>=0.16',
          'matplotlib>=1.4.3',
          'beautifulsoup4>=4.4.1',],
      entry_points={'console_scripts': ['gseapy = gseapy.__main__:main'],},
      
      zip_safe=False,
      download_url='https://github.com/BioNinja/gseapy',)
      
__author__ = 'Zhuoqing Fang'
