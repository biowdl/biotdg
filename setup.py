# Copyright (c) 2020 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from pathlib import Path

from setuptools import find_packages, setup

LONG_DESCRIPTION = Path("README.rst").read_text()

setup(
    name="biotdg",
    version="0.2.0-dev",
    description="Bioinformatics Test Data Generator",
    author="Leiden University Medical Center",
    author_email="sasc@lumc.nl",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/x-rst",
    license="MIT",
    keywords="bioinformatics test data generator creator",
    zip_safe=False,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url="https://github.com/biowdl/biotdg",
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
    ],
    # >=3.6 Because of NamedTuple class constructors
    python_requires=">=3.6",
    install_requires=[
        "biopython",
        "cyvcf2",
        "setuptools"  # To get the version during runtime.
    ],
    entry_points={"console_scripts": [
        "biotdg = biotdg:main"]}
)
