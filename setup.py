import setuptools
from setuptools import setup, find_packages
from codecs import open
from os import path

import evilmc

__version__ = evilmc.__version__

here = path.abspath(path.dirname(__file__))

# get the dependencies and installs
with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '') for x in all_reqs if x.startswith('git+')]

setuptools.setup(
    name="evilmc",
    version="0.1.4",
    url="https://github.com/BoiseStatePlanetary/evilmc",
    download_url='https://github.com/BoiseStatePlanetary/evilmc/archive/'+__version__+'.tar.gz',
    license='BSD',

    author="Brian Jackson",
    author_email="bjackson@boisestate.edu",

    description="A python version of the EVIL-MC code",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages(),

    include_package_data = True,
    package_dir={'evilmc': 'evilmc'}, 
    package_data={'evilmc': ['../data/kepler_response_hires1.txt']},
    install_requires=install_requires,
    dependency_links=dependency_links,

    classifiers=[
        'Development Status :: 2 - Pre-Alpha'
    ],
)
