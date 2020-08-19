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
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    url='https://github.com/wdecoster/nanoget',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore sequencing plotting quality control',
    python_requires='>=3',
    packages=find_packages() + ['scripts'],
    install_requires=['pandas>=0.22.0',
                      'numpy',
                      'biopython',
                      'pysam>0.10.0.0'],
    package_dir={'nanoget': 'nanoget'},
    data_files=[("", ["LICENSE"])])
