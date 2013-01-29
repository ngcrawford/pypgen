import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup

setup(
    name='pypgen',
    version='0.2.0',
    author='Nicholas G. Crawford',
    author_email='ngcrawford@gmail.com',
    url='http://pypi.python.org/pypi/pypgen/',
    license='LICENSE.txt',
    description='Genetic diversity metrics from popoulation genomic datasets.',
    long_description=open('README.rst').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=[
        "numpy >= 1.6.1",
        "pysam >= 0.7.2"
    ],

    packages=[
            'pypgen.fstats',
            'pypgen.parser',
            'pypgen.misc'
            ],

    package_data={
            '': ['*.txt',
                 '*.rst'],     # READMEs, etc
            'pypgen/data': [
                'example.vcf.gz',
                'example.vcf.gz.tbi',
            ]
        },

    include_package_data=True,

    scripts=[
             'scripts/vcf2Dadi.py',
             'scripts/vcf2phylip.py',
             'scripts/vcf_sliding_window.py',
             'scripts/vcf_snpwise_fstats.py',
            ],
)