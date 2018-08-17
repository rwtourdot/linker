## Linker - Analysis tools for Long Read and Linked Read Sequencing

Table of contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Commands](#commands)
  * [Input](#input)
  * [Output](#output)

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

#### Phase haplotypes from long and linked reads

```
./linker phase -i ./input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
```
  * -i: /path/to/input/bamfile
  * -v: /path/to/input/vcffile
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -c: chromosome (chr1,chr2,...,chrX)
  * -s: (optional) second bam file /path/to/second/bamfile
  * -n: (optional) id string for output files

Output is haplotype solution file: haplotype_solution.dat

#### Phase Aneuploid samples based on copy number

```
./linker cn_phase -s haplotype_solution.dat -d het_coverage.dat -n august15
```
  * -s: input haplotype solution
  * -d: input het coverage file
  * -n: (optional) id string for output files
  * -b: (optional) binsize (default is 10kb - 10000)

#### Phase Structural Variants (10X)

```
./linker sv_phase -i ./input.bam -v het_sites.vcf -e tenx -s ./svfile.dat -n august15
```
  * -i: /path/to/input/bamfile
  * -v: /path/to/input/vcffile
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -s: /path/to/svfile
  * -n: (optional) id string for output files

#### Phase Germline with HiC clone data

```
./linker hic_phase -i ./input.bam -v het_sites.vcf -m ./coverage.dat -c chr4 -n august15
```
  * -i: /path/to/input/hic/bamfile
  * -v: /path/to/input/vcffile
  * -m: /path/to/coveragefile
  * -c: chromosome (chr1,chr2,...,chrX)
  * -n: (optional) id string for output files

#### Extract Heterozygous site coveragefile

```
./linker coverage -i ./input.bam -v het_sites.vcf -e tenx -c chr4 -n august15
```
  * -i: /path/to/input/bamfile
  * -v: /path/to/input/vcffile
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -c: chromosome (chr1,chr2,...,chrX)
  * -s: (optional) second bam file /path/to/second/bamfile
  * -n: (optional) id string for output files

#### Create Linked Read Matrix

```
./linker matrix -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```
  * -i: /path/to/input/bamfile
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -c: chromosome (chr1,chr2,...,chrX)
  * -b: (optional) binsize (default is 10kb - 10000)
  * -n: (optional) id string for output files

#### Filter Het Sites by Coverage and Coverage Fraction
```
./linker filter -t tumor_coverage.dat -m normal_coverage.dat -n august15
```
  * -t: aneuploid (tumor) coverage file
  * -m: diploid (normal) coverage file
  * -n: (optional) id string for output files

#### Bin 10X sample by BX tag
```
./linker bx_bin -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```
  * -i: /path/to/input/bamfile
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -c: chromosome (chr1,chr2,...,chrX)
  * -b: (optional) binsize (default is 10kb - 10000)
  * -n: (optional) id string for output files

Input
--------


Output
--------
