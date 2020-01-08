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
from pathlib import Path
from typing import Dict, Generator, Iterable, List, NamedTuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaIterator

import cyvcf2

class Mutation(NamedTuple):
    start: int
    end: int
    sequence: str


def vcf_to_mutations(vcf_path: str, sample: str) -> Dict[str, Dict[int, List[Mutation]]]:
    mutations_dict = {}
    vcf = cyvcf2.VCFReader(vcf_path)
    try:
        vcf.set_samples([sample])
        for record in vcf:  # type: cyvcf2.Variant
            contig = record.CHROM
            start = record.POS - 1
            end = start + len(record.REF)
            if not contig in mutations_dict:
                mutations_dict[contig] = {}
            for allele_no, sequence in enumerate(record.gt_bases):
                mutation = Mutation(start, end, sequence)
                try:
                    mutations_dict[contig][allele_no].append(mutation)
                except KeyError:
                    mutations_dict[contig][allele_no] = [mutation]
    finally:
        vcf.close()

    return mutations_dict


def sequence_with_mutations(sequence: str, mutations: List[Mutation]) -> str:
    # Sort mutations by position
    mutations.sort()
    new_sequence = []
    start = 0

    # This implementation is wrong for overlapping mutations.
    # TODO: Make sure this is documented
    for mutation in mutations:
        new_sequence.append(sequence[start: mutation.start])
        new_sequence.append(mutation.sequence)
        start = mutation.end
    new_sequence.append(sequence[start:])
    return "".join(new_sequence)


def generate_fake_genome(sample: str,
                         reference: Path,
                         vcf_path: str,
                         ploidities: Dict[str, int]
                         ) -> Generator[SeqRecord, None, None]:
    mutations_dict = vcf_to_mutations(vcf_path, sample)
    with reference.open("rt") as reference_h:
        for seqrecord in FastaIterator(reference_h):
            ploidity = ploidities.get(seqrecord.id, 2)
            for allele_no in range(ploidity):
                new_sequence = sequence_with_mutations(str(seqrecord.seq), mutations_dict[seqrecord.id][allele_no])
                new_id = seqrecord.id + "_" + str(allele_no)
                yield SeqRecord(
                    Seq(new_sequence, seqrecord.seq.alphabet),
                    id=new_id,
                    name=new_id,
                    description=new_id)


def write_fasta(seqrecords: Iterable[SeqRecord], filepath: Path):
    filepath.parent.mkdir(parents=True)
    with filepath.open("wt") as file_h:
        for record in seqrecords:
            file_h.write(record.format("fasta"))
    return filepath


def dwgsim(*args):
    """
    Runs `dwgsim` locally with specified args.
    :param args:
    :return:
    """
    subprocess.run(["dwgsim"] + list(args))


def ploidity_file_to_dict(ploidity_file: Path) -> Dict[str, int]:
    ploidity_dict = dict()
    with ploidity_file.open("r") as file_h:
        for line in file_h:
            chromosome, ploidity = line.strip().split("\t")
            ploidity_dict[chromosome] = int(ploidity)
    return ploidity_dict


def generate_test_data(sample: str,
                       reference: Path,
                       vcf_path: str,
                       ploidity_file: Path,
                         ) -> Generator[SeqRecord, None, None]



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
    parser.add_argument("-o", "--output-dir", type=str)
    return parser



def main():
    args = argument_parser().parse_args()


if __name__ == "__main__":
    main()