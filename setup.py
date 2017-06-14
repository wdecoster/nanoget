# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

exec(open('nanoget/version.py').read())


setup(
    name='nanoget',
    version=__version__,
    description='Functions to extract information from Oxford Nanopore sequencing data and alignments.',
    long_description='Functions to extract information from Oxford Nanopore sequencing data and alignments.',
    url='https://github.com/wdecoster/nanoget',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore sequencing plotting quality control',
    packages=find_packages(),
    install_requires=['pandas', 'numpy', 'biopython', 'pysam>0.10.0.0', 'nanomath'],
    package_dir={'nanoget': 'nanoget'})
