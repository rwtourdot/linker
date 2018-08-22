## Linker - Analysis tools for Long Read and Linked Read Sequencing

Table of contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Commands](#commands)
  * [Input/Output](#input/output)

Installation
------------

code requires bamtools, htslib, c++11, and zlib

  * htslib: https://github.com/samtools/htslib
  * bamtools: https://github.com/pezmaster31/bamtools

From linker directory:

Installing htslib locally
```
git clone https://github.com/samtools/htslib
cd htslib
autoheader
autoconf
./configure --prefix=/path/to/linker/packages/htslib/
make
make install
cd ..
```

Installing bamtools locally
```
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/linker/packages/bamtools/ ..
make
make install
cd ../..
```

Description
-----------

<img src="https://github.com/rwtourdot/linker/blob/master/linker_flowchart.png" width=800/>

Commands
--------

#### List Commands
```
./linker -h
```

#### Flags
  * -i: /path/to/input/bamfile
  * -v: /path/to/input/vcffile
  * -u: /path/to/input/svfile
  * -l: /path/to/haplotype_solution_file
  * -t: aneuploid (tumor) coverage file
  * -m: diploid (normal) coverage file
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -c: chromosome (chr1,chr2,...,chrX)
  * -s: (optional) second bam file /path/to/second/bamfile
  * -n: (optional) id string for output files
  * -b: (optional) binsize (default is 10kb - 10000)

#### Phase germline haplotypes from long and linked reads

This commmand takes in a vcf file and a long or linked read bam to compute phased haplotype blocks.  The vcf file should contain all germline heterozygous sites and most likely originates from a paired normal sample.  The bam file or files should be obtained with a long read technology and could originate for tumor,normal, or tumor+normal.

```
./linker phase -i ./input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
./linker phase -i ./input.bam -s ./second_input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
```
  * Output is haplotype solution file: haplotype_solution.dat

The output of this command is a file which contains the minimum energy solution to the germline haplotype.  More information on this file is described in the I/O section below.

#### Extract Heterozygous site Coverage

This command takes a vcf and bam file and extracts the read coverage of each allele. In order to count a base at a heterozygous site both the base quality and read map quality must pass a cutoff.

```
./linker coverage -i ./input.bam -v het_sites.vcf -e illumina -c chr4 -n august15
```
  * Output is heterozygous site coverage file: coverage.dat

The output of this command is a file which contains the read depth coverage of each heterozygous site for both reference and variant bases.  More information on this file is described in the I/O section below.

#### Phase Aneuploid samples based on copy number

Haplotypes can be phased further based on tumor copy number.  Copy number phasing works better in samples where aneuploidy and loss of heterozygocity is prevalent.

```
./linker cn_phase -l haplotype_solution.dat -m ./het_coverage.dat -n august15
```
  * Output is output is copy number phased haplotype file: cn_haplotype.dat

#### Phase Structural Variants (10X)

```
./linker sv_phase -i ./input.bam -v het_sites.vcf -e tenx -u ./svfile.dat -n august15
```

#### Phase Germline with HiC clone data

```
./linker hic_phase -i ./input.bam -m ./coverage.dat -c chr4 -n august15
```

#### Create Linked Read Matrix

```
./linker matrix -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```

#### Filter Het Sites by Coverage and Coverage Fraction
```
./linker filter -t tumor_coverage.dat -m normal_coverage.dat -n august15
```

#### Bin 10X sample by BX tag
```
./linker bx_bin -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```

Input/Output
--------

#### Required Files

#### Generated Files
