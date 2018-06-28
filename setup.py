import setuptools
from setuptools import setup, find_packages
from codecs import open
from os import path

import evilmc

__version__ = evilmc.__version__

setuptools.setup(
    name="evilmc",
    version=__version__,
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
    package_data={'evilmc': ['kepler_response_hires1.txt']},
    install_requires=['numpy', 'PyAstronomy', 'six', 'astropy'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha'
    ],
)
