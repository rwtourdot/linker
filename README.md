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

This commmand takes in a vcf file and a long or linked read bam to compute phased haplotype blocks

```
./linker phase -i ./input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
```
  * Output is haplotype solution file: haplotype_solution.dat

#### Extract Heterozygous site Coverage

```
./linker coverage -i ./input.bam -v het_sites.vcf -e tenx -c chr4 -n august15
```

#### Phase Aneuploid samples based on copy number

```
./linker cn_phase -l haplotype_solution.dat -m ./het_coverage.dat -n august15
```

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
