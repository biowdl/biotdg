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

import os
import shutil
import sys
import tempfile
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotdg import (Mutation, dwgsim, generate_fake_genome, generate_test_data,
                    main, ploidy_file_to_dict, sequence_with_mutations,
                    vcf_to_mutations, write_fasta)

TEST_DATA = Path(__file__).parent / Path("data")


def test_vcf_to_mutations():
    vcf = TEST_DATA / Path("test.vcf")
    vcf_mut_dict = {
        "chr1": {
            0: [
                Mutation(499, 500, "C"),
                Mutation(999, 1000, "A"),
                Mutation(1499, 1500, "A"),
                Mutation(1999, 2000, "A")
            ],
            1: [
                Mutation(499, 500, "C"),
                Mutation(999, 1000, "G"),
                Mutation(1499, 1500, "G"),
                Mutation(1999, 2000, "A")
            ]
        }
    }
    result_dict = vcf_to_mutations(str(vcf), "wgs2")
    assert vcf_mut_dict == result_dict


def test_sequence_with_mutations():
    # Function does not test if sequence is DNA, so we can use numbers and
    # symbols for easy testing.
    mutations = [Mutation(5, 6, "_"), Mutation(8, 9, "=")]
    sequence = "0123456789"
    mutated_seq = sequence_with_mutations(sequence, mutations)
    assert mutated_seq == "01234_67=9"


def test_generate_fake_genome():
    reference = TEST_DATA / Path("reference.fasta")
    vcf = TEST_DATA / Path("empty.vcf")
    ploidy_dict = {"chr1": 3, "chrX": 2, "chrY": 1}
    sequences = (generate_fake_genome("sample1",
                                      reference,
                                      vcf,
                                      ploidy_dict))
    ids_and_seqs = [(record.id, str(record.seq)) for record in sequences]
    assert ids_and_seqs == [
        ("chr1_0", "GATTACAGATTACAGATTACA"),
        ("chr1_1", "GATTACAGATTACAGATTACA"),
        ("chr1_2", "GATTACAGATTACAGATTACA"),
        ("chrX_0", "AGTCAGTCAGTC"),
        ("chrX_1", "AGTCAGTCAGTC"),
        ("chrY_0", "AGAATC"),
    ]


def test_generate_fake_genome_with_mutations():
    reference = TEST_DATA / Path("reference.fasta")
    vcf = TEST_DATA / Path("sample1.vcf")
    ploidy_dict = {"chr1": 3, "chrX": 2, "chrY": 1}
    sequences = (generate_fake_genome("sample1",
                                      reference,
                                      vcf,
                                      ploidy_dict))
    ids_and_seqs = [(record.id, str(record.seq)) for record in sequences]
    assert ids_and_seqs == [
        ("chr1_0", "GATAACAGATTACAGATTACA"),
        ("chr1_1", "GATCACTGATTACAGATTACA"),
        ("chr1_2", "GATGACAGATTACAGATTACA"),
        ("chrX_0", "AGTCAGTCAGTC"),
        ("chrX_1", "TGTCAGTCAGTC"),
        ("chrY_0", "AGACTC"),
    ]


def test_write_fasta():
    _, tmppath = tempfile.mkstemp()

    tmpfile = Path(tmppath)
    seqrecords = [SeqRecord(Seq("AAAAA"), "chr1", description=""),
                  SeqRecord(Seq("TTTTT"), "chr2", description="")]
    write_fasta(seqrecords, tmpfile)
    assert tmpfile.read_text() == ">chr1\nAAAAA\n>chr2\nTTTTT\n"
    os.remove(str(tmpfile))


def test_dwgsim():
    """
    Test if files are properly generated. Also test if mutation rate of
    0.0 does not generate mutations.
    """
    tempdir = Path(tempfile.mkdtemp())
    prefix = str(Path(tempdir, "bla"))
    reference = TEST_DATA / Path("a.fasta")
    dwgsim(in_ref_fa=str(reference), out_prefix=prefix, mutation_rate=0.0,
           per_base_error_rate_read1=0.001, per_base_error_rate_read2=0.001,
           length_read1=150, length_read2=150, random_seed=17,
           probability_random_dna_read=0.2)
    generated_files = os.listdir(str(tempdir))
    assert "bla.bwa.read1.fastq" in generated_files
    assert "bla.bwa.read2.fastq" in generated_files
    mutations_txt = tempdir / Path("bla.mutations.txt")
    assert mutations_txt.read_text() == ""
    # Manual testing was performed: A mutation rate of 0.1 does produce a
    # mutations.txt file with contents.
    # Despite having a mutation rate of 0.0, there are still read errors.
    # This is easily verified as this reference only contains A's
    # This means dwgsim is used in a valid way in this program.
    shutil.rmtree(str(tempdir))


def test_ploidy_file_to_dict():
    ploidy_file = TEST_DATA / Path("ploidy_file.tsv")
    assert ploidy_file_to_dict(ploidy_file) == {
        "chr1": 3,
        "chrX": 2,
        "chrY": 1
    }


def test_generate_test_data():
    reference = TEST_DATA / Path("reference.fasta")
    vcf = TEST_DATA / Path("sample1.vcf")
    ploidy_file = TEST_DATA / Path("ploidy_file.tsv")
    output_dir = Path(tempfile.mkdtemp())
    generate_test_data(sample="sample1",
                       reference=reference,
                       ploidy_file=ploidy_file,
                       vcf_path=vcf,
                       random_seed=17,
                       output_dir=output_dir)
    assert (output_dir / Path("sample1.fasta")).exists()
    shutil.rmtree(output_dir)


def test_main():
    output_dir = Path(tempfile.mkdtemp())
    sys.argv = ["biotdg",
                "-r", str(TEST_DATA / Path("reference.fasta")),
                "--vcf", str(TEST_DATA / Path("sample1.vcf")),
                "-p", str(TEST_DATA / Path("ploidy_file.tsv")),
                "-s", "sample1",
                "-z", "17",
                "-o", str(output_dir)]
    main()
    assert Path(output_dir, "sample1.fasta").exists()
    shutil.rmtree(output_dir)
