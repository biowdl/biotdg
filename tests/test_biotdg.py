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

from biotdg import Mutation, generate_fake_genome, sequence_with_mutations, \
    vcf_to_mutations

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
    ploidities = {"chr1": 3, "chrX": 2, "chrY": 1}
    sequences = (generate_fake_genome("sample1",
                                      reference,
                                      vcf,
                                      ploidities))
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
    ploidities = {"chr1": 3, "chrX": 2, "chrY": 1}
    sequences = (generate_fake_genome("sample1",
                                      reference,
                                      vcf,
                                      ploidities))
    ids_and_seqs = [(record.id, str(record.seq)) for record in sequences]
    assert ids_and_seqs == [
        ("chr1_0", "GATAACAGATTACAGATTACA"),
        ("chr1_1", "GATCACTGATTACAGATTACA"),
        ("chr1_2", "GATGACAGATTACAGATTACA"),
        ("chrX_0", "AGTCAGTCAGTC"),
        ("chrX_1", "TGTCAGTCAGTC"),
        ("chrY_0", "AGACTC"),
    ]