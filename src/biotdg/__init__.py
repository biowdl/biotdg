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
from typing import Dict, Generator, Iterable, List, NamedTuple, Optional

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqRecord import SeqRecord

import cyvcf2


class Mutation(NamedTuple):
    start: int
    end: int
    sequence: str


def vcf_to_mutations(vcf_path: str, sample: str
                     ) -> Dict[str, Dict[int, List[Mutation]]]:
    mutations_dict: Dict[str, Dict[int, List[Mutation]]] = {}
    vcf = cyvcf2.VCFReader(vcf_path)
    try:
        vcf.set_samples([sample])
        for record in vcf:  # type: cyvcf2.Variant
            contig = record.CHROM
            start = record.start
            end = record.end
            genotypes = [record.REF] + record.ALT
            # We take the genotypes from the first sample. Because there is
            # only one sample. We take all genotypes except the last column
            # because that is a boolean that indicates if it is phased or not.
            genotype_indexes = record.genotypes[0][:-1]
            if contig not in mutations_dict.keys():
                mutations_dict[contig] = {}
            # Get all the genotypes without the boolean for phasing
            for allele_no, genotype_idx in enumerate(genotype_indexes):
                mutation = Mutation(start, end, genotypes[genotype_idx])
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
                         vcf_path: Path,
                         ploidy_dict: Dict[str, int]
                         ) -> Generator[SeqRecord, None, None]:
    mutations_dict = vcf_to_mutations(str(vcf_path), sample)
    with reference.open("rt") as reference_h:
        for seqrecord in FastaIterator(reference_h):
            ploidy = ploidy_dict.get(seqrecord.id, 2)
            for allele_no in range(ploidy):
                # Default to empty list if no mutations were listed.
                mutations = mutations_dict.get(seqrecord.id, {}
                                               ).get(allele_no, [])
                new_sequence = sequence_with_mutations(
                    sequence=str(seqrecord.seq),
                    mutations=mutations)
                new_id = seqrecord.id + "_" + str(allele_no)
                yield SeqRecord(
                    Seq(new_sequence, seqrecord.seq.alphabet),
                    id=new_id,
                    name=new_id,
                    description=new_id)


def write_fasta(seqrecords: Iterable[SeqRecord], filepath: Path):
    # Make sure output directory is present.
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("wt") as file_h:
        for record in seqrecords:
            file_h.write(record.format("fasta"))
    return filepath


def dwgsim(*args,
           in_ref_fa: str,
           out_prefix: str,
           per_base_error_rate_read1: Optional[str] = None,
           per_base_error_rate_read2: Optional[str] = None,
           length_read1: Optional[int] = None,
           length_read2: Optional[int] = None,
           mean_coverage: Optional[float] = None,
           mutation_rate: Optional[float] = None,
           probability_random_dna_read: Optional[float] = None,
           maximum_n_number: Optional[int] = None,
           random_seed: Optional[int] = None):
    argslist = list(args)
    if per_base_error_rate_read1 is not None:
        argslist.extend(["-e", str(per_base_error_rate_read1)])
    if per_base_error_rate_read2 is not None:
        argslist.extend(["-E", str(per_base_error_rate_read2)])
    if length_read1 is not None:
        argslist.extend(["-1", str(length_read1)])
    if length_read2 is not None:
        argslist.extend(["-2", str(length_read2)])
    if mean_coverage is not None:
        argslist.extend(["-C", str(mean_coverage)])
    if mutation_rate is not None:
        argslist.extend(["-r", str(mutation_rate)])
    if probability_random_dna_read is not None:
        argslist.extend(["-y", str(probability_random_dna_read)])
    if maximum_n_number is not None:
        argslist.extend(["-n", str(maximum_n_number)])
    if random_seed is not None:
        argslist.extend(["-z", str(random_seed)])
    argslist.extend([in_ref_fa, out_prefix])
    subprocess.run(["dwgsim"] + list(argslist))


def ploidy_file_to_dict(ploidy_file: Path) -> Dict[str, int]:
    ploidy_dict = dict()
    with ploidy_file.open("r") as file_h:
        for line in file_h:
            chromosome, ploidy = line.strip().split("\t")
            ploidy_dict[chromosome] = int(ploidy)
    return ploidy_dict


def generate_test_data(sample: str,
                       reference: Path,
                       vcf_path: Path,
                       ploidy_file: Path,
                       output_dir: Path,
                       random_seed: Optional[int] = None,
                       read_length: Optional[int] = None,
                       coverage: Optional[float] = None,
                       read1_error_rate: Optional[str] = None,
                       read2_error_rate: Optional[str] = None,
                       maximum_n_number: Optional[int] = None):
    # Make sure output directory is present.
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_genome = generate_fake_genome(
        sample, reference, vcf_path, ploidy_file_to_dict(ploidy_file))

    generated_genome_path = Path(output_dir, sample + ".fasta")
    write_fasta(generated_genome, generated_genome_path)
    dwgsim(in_ref_fa=str(generated_genome_path),
           out_prefix=str(Path(output_dir, sample)),
           mutation_rate=0.0,
           random_seed=random_seed,
           length_read1=read_length,
           length_read2=read_length,
           maximum_n_number=maximum_n_number,
           mean_coverage=coverage,
           per_base_error_rate_read1=read1_error_rate,
           per_base_error_rate_read2=read2_error_rate
           )


def argument_parser() -> argparse.ArgumentParser:
    """
    Creates the argument parser for biotdg
    :return: a parser object
    """
    parser = argparse.ArgumentParser(
        description="Bioinformatics Test Data Generator")
    parser.add_argument("-r", "--reference", type=Path, required=True,
                        help="Reference genome for the sample.")
    parser.add_argument("--vcf", required=True,
                        help="VCF file with mutations.")
    parser.add_argument("-p", "--ploidy-table", type=Path, required=True,
                        help="Tab-delimited file with two columns specifying "
                             "the chromosome name and its ploidy. By "
                             "default all chromosomes have a ploidy of 2")
    parser.add_argument("-s", "--sample-name", type=str, required=True,
                        help="name of the sample to generate. The sample must "
                             "be in the VCF file")
    parser.add_argument("-z", "--random-seed", type=int, default=1,
                        help="random seed for dwgsim (default: 1)")
    parser.add_argument("-l", "--read-length", type=int, default=150,
                        help="read length to be used by dwgsim")
    parser.add_argument("-C", "--coverage", type=float, default=50,
                        help="Average coverage for the generated reads. NOTE: "
                             "This is multiplied by the ploidy of the "
                             "chromosome.")
    parser.add_argument("-e", "--read1-error-rate", type=str,
                        help="Same as -e flag in dwgsim. per base/color/flow "
                             "error rate of the first read.")
    parser.add_argument("-E", "--read2-error-rate", type=str,
                        help="Same as -E flag in dwgsim. per base/color/flow "
                             "error rate of the second read.")
    parser.add_argument("-n", "--maximum-n-number", type=int,
                        help="Maximum number of Ns allowed in a given read.")
    parser.add_argument("-o", "--output-dir", type=Path)
    return parser


def main():
    args = argument_parser().parse_args()
    generate_test_data(sample=args.sample_name,
                       reference=args.reference,
                       vcf_path=args.vcf,
                       ploidy_file=args.ploidy_table,
                       output_dir=args.output_dir,
                       random_seed=args.random_seed,
                       read_length=args.read_length,
                       coverage=args.coverage,
                       read1_error_rate=args.read1_error_rate,
                       read2_error_rate=args.read2_error_rate,
                       maximum_n_number=args.maximum_n_number
                       )
