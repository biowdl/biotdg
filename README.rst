.. Badges have empty alts. So nothing shows up if they do not work.

.. image:: https://img.shields.io/pypi/v/biotdg.svg
  :target: https://pypi.org/project/biotdg/
  :alt:

.. image:: https://img.shields.io/conda/v/bioconda/biotdg.svg
  :target: http://bioconda.github.io/recipes/biotdg/README.html
  :alt:

.. image:: https://img.shields.io/pypi/pyversions/biotdg.svg
  :target: https://pypi.org/project/biotdg/
  :alt:

.. image:: https://img.shields.io/pypi/l/biotdg.svg
  :target: https://github.com/biowdl/biotdg/blob/master/LICENSE
  :alt:

.. image:: https://travis-ci.org/biowdl/biotdg.svg?branch=develop
  :target: https://travis-ci.org/biowdl/biotdg
  :alt:

.. image:: https://codecov.io/gh/biowdl/biotdg/branch/develop/graph/badge.svg
  :target: https://codecov.io/gh/biowdl/biotdg
  :alt:

biotdg: Bioinformatics Test Data Generator
==========================================

``biotdg`` can generate mutations based on vcf files for genomes where the
chromosomes have different ploidy. It was made to create test genomes for
pipelines that correctly handle the ploidy of sex chromosomes. It can also be
used to create test data for pipelines that handle triploid species such as
banana, or for pipelines that discover chromosome inbalances such as
trisomy-21 (Down syndrome) and XXY males (Klinefelter syndrome).

``biotdg`` uses a reference genome, a ploidy table and a vcf file to create a
"true genome" for a sample. For example if the ploidy table states that
``chr21`` has a ploidy of 3 then the "true genome" will have three copies
of ``chr21``. Each ``chr21`` copy will have its own mutations based on the
vcf file.

After creating the "true genome" fasta file. ``biotdg`` uses the
`dwgsim <https://github.com/nh13/dwgsim>`_ program to generate fastq reads.

Usage
-----

.. code-block:: text

    usage: biotdg [-h] -r REFERENCE --vcf VCF -p PLOIDY_TABLE -s SAMPLE_NAME
                  [-z RANDOM_SEED] [-o OUTPUT_DIR]

    Bioinformatics Test Data Generator

    optional arguments:
      -h, --help            show this help message and exit
      -r REFERENCE, --reference REFERENCE
                            Reference genome for the sample.
      --vcf VCF             VCF file with mutations.
      -p PLOIDY_TABLE, --ploidy-table PLOIDY_TABLE
                            Tab-delimited file with two columns specifying the
                            chromosome name and its ploidity. By default all
                            chromosomes have a ploidity of 2
      -s SAMPLE_NAME, --sample-name SAMPLE_NAME
                            name of the sample to generate. The sample must be in
                            the VCF file
      -z RANDOM_SEED, --random-seed RANDOM_SEED
                            random seed for dwgsim (default: 1)
      -o OUTPUT_DIR, --output-dir OUTPUT_DIR


Example
-------
Given the following ``reference.fasta`` file

.. code-block:: text

    >chr1
    GATTACA
    GATTACA
    GATTACA
    >chrX
    AGTCAGTCAGTC
    >chrY
    AGAATC

the following ploidy table.tsv

.. code-block:: text

    chr1	3
    chrX	2
    chrY	1

and the following vcf:

.. code-block:: text

    ##fileformat=VCFv4.1
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr1,length=21>
    ##contig=<ID=chrX,length=12>
    ##contig=<ID=chrY,length=6>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
    chr1	4	.	T	A,C,G	.	.	.	GT	1/2/3
    chr1	7	.	A	T	.	.	.	GT	0/1/0
    chrX	1	.	A	T	.	.	.	GT	0/1
    chrX	2	.	G	T	.	.	.	GT	0/0
    chrY	4	.	A	C	.	.	.	GT	1

        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr1,length=21>
    ##contig=<ID=chrX,length=12>
    ##contig=<ID=chrY,length=6>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
    chr1	4	.	T	A,C,G	.	.	.	GT	1/2/3
    chr1	7	.	A	T	.	.	.	GT	0/1/0
    chrX	1	.	A	T	.	.	.	GT	0/1
    chrX	2	.	G	T	.	.	.	GT	0/0
    chrY	4	.	A	C	.	.	.	GT	1

A "true genome" for sample1 looks like this:

.. code-block:: text

    >chr1_0
    GATAACAGATTACAGATTACA
    >chr1_1
    GATCACTGATTACAGATTACA
    >chr1_2
    GATGACAGATTACAGATTACA
    >chrX_0
    AGTCAGTCAGTC
    >chrX_1
    TGTCAGTCAGTC
    >chrY_0
    AGACTC

Why ``biotdg`` and not ``dwgsim``?
----------------------------------

``dwgsim`` has excellent capabilities for generating reads that are close to
real data. Therefore ``dwgsim`` is used by ``biotdg`` in this capacity.

``dwgsim`` can also generate mutations randomly and output these in VCF format.
It also has the capability to use a VCF to generate mutations. This VCF-based
method was not deemed sufficient for the following reasons:

+ Very poorly documented.
+ Only allows ploidy of 1 or 2. There is an option '3' but that does something
  different.
+ How exactly mutations are generated is unknown. Is it aware of phasing? If
  so, how does it handle it?

``biotdg`` handles the creation of the "true genome" transparently and then
uses dwgsim to generate reads.

Known limitations
-----------------
+ Overlapping mutations are not handled properly. (Probably not a concern for
  generating test data.)

