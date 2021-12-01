#! /usr/bin/env python

import os
from setuptools import setup
import sys

PACKAGE = "seqcol"

# Additional keyword arguments for setup().
extra = {}

# Ordinary dependencies
DEPENDENCIES = []
with open("requirements/requirements-all.txt", "r") as reqs_file:
    for line in reqs_file:
        if not line.strip():
            continue
        DEPENDENCIES.append(line)

extra["install_requires"] = DEPENDENCIES

with open("{}/_version.py".format(PACKAGE), 'r') as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")

with open('README.md') as f:
    long_description = f.read()

setup(
    name=PACKAGE,
    packages=[PACKAGE],
    version=version,
    description="Python implementation of seqcol protocol",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: System :: Distributed Computing"
    ],
    keywords="genome, assembly, bioinformatics, reference, sequence",
    url="https://github.com/refgenie/seqcol",
    author=u"Nathan Sheffield, Michal Stolarczyk",
    author_email=u"nathan@code.databio.org", 
    license="BSD2",
    entry_points={
        "console_scripts": [
            'seqcol = seqcol.seqcol:main'
        ],
    },    
    include_package_data=True,
    test_suite="tests",
    tests_require=(["mock", "pytest"]),
    setup_requires=(["pytest-runner"] if {"test", "pytest", "ptr"} & set(sys.argv) else []),
    **extra
)
