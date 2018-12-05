## Linker - Tools for Analyzing Long and Linked read Sequencing

Table of Contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Commands](#commands)
  * [Input/Output](#input/output)

Installation
------------

This code requires bamtools, htslib, c++11, and zlib libraries.

  * htslib: https://github.com/samtools/htslib
  * bamtools: https://github.com/pezmaster31/bamtools

From the linker directory:

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
Linker is a suite of tools useful for interpreting long or linked read sequencing of cancer genomes.  The key information long read sequencing provides is the germline haplotype of the cell line.  In cancer cell lines where Aneuplpoidy, Loss of Heterozygosity, and Structural Variation is common, haplotypes can provide a better resolution of the samples karyotype, and clarify the cancer cells genomic evolution.

Linker currently supports 10X, Pacbio, Oxford Nanopore, and HiC sequencing technologies.  Some technologies work better in the specific tools listed below, and HiC requires a distinct phasing methodology.  A diagram of the analysis workflow is shown below.

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

#### Phase Germline Haplotypes from Long and Linked Reads

This commmand takes in a vcf file and a long or linked read bam file to compute phased haplotype blocks.  The vcf file should contain all germline heterozygous sites and most likely originates from a paired normal sample.  The bam file or files should be obtained with a long read technology and could originate for tumor, normal, or tumor+normal.

```
./linker phase -i ./input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
./linker phase -i ./input.bam -s ./second_input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
```
  * Output is haplotype solution file: haplotype_solution.dat

The output of this command is a file which contains the minimum energy solution to the germline haplotype.  More information on this file is described in the I/O section below.

#### Extract Heterozygous Site Coverage

This command takes a vcf and bam file and extracts the read coverage of each allele. In order to count a base at a heterozygous site both the base quality and read map quality must pass a cutoff.

```
./linker coverage -i ./input.bam -v het_sites.vcf -e illumina -c chr4 -n august15
```
  * Output is heterozygous site coverage file: coverage.dat

The output of this command is a file which contains the read depth coverage of each heterozygous site for both reference and variant bases.  More information on this file is described in the I/O section below.

#### Phase Aneuploid Samples based on Copy Number

Haplotypes can be phased further based on tumor copy number.  Copy number phasing works better in tumor samples where aneuploidy and loss of heterozygocity is prevalent.

```
./linker cn_phase -l haplotype_solution.dat -m ./het_coverage_aneuploid.dat -n august15
```
  * Output is a copy number phased haplotype file: cn_haplotype.dat

#### Phase Structural Variants (10X/Nanopore)

Once haplotypes are found, associated Structural Variants can be phased with a 10X or Nanopore tumor sample.  Structual Variant's can be called with a SV caller and converted to the svfile input file format described in the I/O section below.

```
./linker sv_phase -i ./input.bam -v het_sites.vcf -e tenx -u ./svfile.dat -n august15
```
  * Output is a phased sv file: sv_phased.dat

#### Phase Germline with HiC (clonal sample)

HiC data can be phased in a similar fashon to SV's.  This command takes an input HiC bam file and a vcf file or coverage file and phases any chromosome contacts which overlap a het site.  HiC data is sparse and a connection of two het sites does not guarantee they are on the same allele.  HiC phasing data is therefore probibalistic and requires another long or linked read technology to phase accurately.

```
./linker hic_phase -i ./input.bam -m ./coverage.dat -c chr4 -n august15
```
  * Output is a phased hic het link file: hic_links.dat

#### Create Linked Read Matrix

A local alignment map can be extracted from any long or linked read technology. This map does not contain any allelic information but can more clearly show translocations and inversions.

```
./linker matrix -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```
  * Output is a matrix file: matrix.dat

#### Filter Het Sites by Coverage and Coverage Fraction

In order refine variant calls it can be useful to use allele fraction from a normal sample to extract a more confident subset of variants.  This command takes in a coverage file and filters het sites based on allelic fraction.  More specifically, this command checks if the fraction of variant bases is between 10 - 90 percent and passes it to a filtered output file.  This command can also filter tumor coverage data based on the allele fraction in a corresponding normal sample.

```
./linker filter -m normal_coverage.dat -n august15
./linker filter -t tumor_coverage.dat -m normal_coverage.dat -n august15
```
  * Output is a filtered coverage file: filtered_coverage.dat

#### Bin 10X Sample by BX tag

10X sequencing coverage data can be refined by looking at unique bx tag density.  This command takes a 10X bam file and bins the genome by coverage of unique bx tag.

```
./linker bx_bin -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```

Input/Output
--------

#### Required Files

##### .bam file
The bam file can originate from Illumina or another type of sequencing. The alignment method will depend on the type of sequencing, but a consistent genome reference should be used between samples.

##### .vcf file
The vcf file can be obtained with GATK (https://software.broadinstitute.org/gatk/) or another variant caller and should contain all germline and somatic SNPs.  Indels will not be considered in linkers phasing methods.  Variants should be called on the same reference as the long or linked read .bam file

##### SV file
There are several Structural Variant callers that work with paired end sequencing. The output of their calls may differ, linker takes structural variant input in the form below.

```
RAindex	chr1	pos1	str1	chr2	pos2	str2	TotalCount
20  chr1	1651206	-1	chr1	1717357	1	3
```

#### Generated Files

##### het coverage file
The heterozygous site coverage file contains information on the number of each base at every heterozygous site called with the GATK.
```
variant_id  position  variant:reference Indel Deletion  Gbase Cbase Abase Tbase Total
190204289_T_C	190204289	T:C	I|0	D|0	G|0	C|44	A|0	T|30	74
```
