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
    test_suite='src/tests',
    # install_requires=[
    #     "numpy == 1.6.1",
    #     "pysam == 0.7.2"
    # ],

    packages=[
            'src.fstats',
            'src.parser',
            'src.misc'
            ],

    package_data={
            '': ['*.txt',
                 '*.rst'],     # READMEs, etc
            'src/data': [
                'butterfly.vcf.gz',
                'butterfly.vcf.gz.tbi',
            ]
        },

    include_package_data=True,

    scripts=[
             'src/vcf2Dadi.py',
             'src/vcf2phylip.py',
             'src/vcf_sliding_window.py',
             'src/vcf_snpwise_fstats.py',
            ],
)
