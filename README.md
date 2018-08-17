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

#### Phase Aneuploid samples based on copy number

```
./linker cn_phase -s haplotype_solution.dat -d het_coverage.dat -n august15
```
  * -s: /path/to/input/bamfile
  * -d: /path/to/input/vcffile
  * -n: (optional) id string for output files
  * -b: (optional) binsize (default is 10kb - 10000)

#### Phase Aneuploid samples based on copy number

```
./linker sv_phase -h
```
  * -i: /path/to/input/bamfile
  * -v: /path/to/input/vcffile
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -s: /path/to/svfile
  * -n: (optional) id string for output files

```
## phase Hi-C
./linker hic_phase -h

## Extract heterozygous site coverage
./linker coverage -h

## create linked read heatmap
./linker matrix -h

## filter het sites based on coverage and coverage fraction
./linker filter -h

## 10X specific - find droplet based genome covereage
./linker bx_bin -h
```

Input
--------


Output
--------
