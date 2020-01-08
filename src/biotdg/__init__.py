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

import argparse
import subprocess
from typing import NamedTuple

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaIterator


class Mutation(NamedTuple):
    contig: str
    start: int
    end: int
    sequence: str


def argument_parser() -> argparse.ArgumentParser:
    """
    Creates the argument parser for biotdg
    :return: a parser object
    """
    parser = argparse.ArgumentParser(
        description="Bioinformatics Test Data Generator")
    parser.add_argument("--vcf", required=True,
                        help="VCF file with mutations.")
    parser.add_argument("-p", "--ploidy-table", required=True,
                        help="Tab-delimited file with two columns specifying "
                             "the chromosome name and its ploidity. By "
                             "default all chromosomes have a ploidity of 2")
    parser.add_argument("-s", "--sample-name", required=True,
                        help="name of the sample to generate. The sample must "
                             "be in the VCF file")
    return parser


def dwgsim(*args):
    """
    Runs `dwgsim` locally with specified args.
    :param args:
    :return:
    """
    subprocess.run(["dwgsim"] + list(args))


def main():
    args = argument_parser().parse_args()


if __name__ == "__main__":
    main()