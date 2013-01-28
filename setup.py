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
    test_suite='vcf_tests.py',
    install_requires=[
        "numpy >= 1.6.1",
        "pysam >= 0.6"
    ],

    packages=[
            'fstats',
            ],

    package_data = {
            '': ['*.txt'],     # READMEs, etc
            'test_data': [
                'butterfly.vcf.gz',
                'butterfly.vcf.gz.tbi',
            ]
        },
    include_package_data = True,
    scripts = [
              'cloudforest/nexus2oneliner.py',
              'cloudforest/phylip2oneliner.py',
              'cloudforest/process.py'
              ],
)
