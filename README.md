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
Linker is a suite of c++ tools useful for interpreting long and linked read sequencing of cancer genomes.  The most significant information long read sequencing provides is the local haplotype of a sample.  In cancer cell lines where Aneuplpoidy, Loss of Heterozygosity, and Structural Variation is common, haplotypes can provide a better resolution of the samples karyotype, and clarify the cancer cells genomic evolution.

Linker currently supports 10X, Pacbio, Oxford Nanopore, and HiC sequencing technologies.  A diagram of the germline haplotype phasing workflow is shown below.

<img src="https://github.com/rwtourdot/linker/blob/master/new_linker_flowchart.png" width=600/>

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

#### Extract Phasing Information from Long and Linked Reads

This command extracts all long and linked read phasing information from an aligned bamfile given a corresponding vcf file containing heterozygous sites.  The long read technology flag (tenx,pacbio,nanopore,hic) should correspond to the bamfile chosen.

```
./linker extract -i ./input.bam -v het_sites.vcf -e tenx -c chr21 -n sample_name
```
  * Output is a graph file: graph_variant_{}.dat

The output of this command is a graph_variant file which lists all of the unique hashes associated with each het site base.  This file can be concatenated with other graph_variant files, to combine all phasing information from multiple technologies.  

#### Solve For Germline Haplotype from Long Read Links

After extracting all phasing information into graph_variant files and combining or trimming certain hashes the samples haplotype can be solved for by:

```
./linker solve -i ./graph_variant_{}_chr21.dat -c chr21 -n sample_name
```
  * Output is the haplotype solution file: hap_solution_{}.dat

The haplotype file contains a Block Switch Energy column which can be used to define blocks.  The lower (more negative) the block energy the more likely a heterozygous site is phased correctly. More details on defining blocks can be found in the paper.

#### Generate A Whole Chromosome Haplotype Scaffold

This command takes a haplotype solution file and combines it with hic phasing information in a corresponding graph_variant_hic file to generate a full chromosome haplotype scaffold.  An energy cutoff which defines haplotype blocks is specified by the -e flag.

```
./linker scaffold -i ./hap_solution_{}_chr21.dat -g ./graph_variant_hic_chr21.dat -e -700 -c chr21 -n sample_name
```
  * Output is a haplotype scaffold file: hap_full_scaffold_{}.dat

The hap_full_scaffold file contains less heterozygous sites than the haplotype solution file but is accurate over the full length of the chromosome.

<!--
#### Phase Germline Haplotypes from Long and Linked Reads

This commmand takes in a vcf file and a long or linked read bam file to compute phased haplotype blocks.  The vcf file should contain all germline heterozygous sites and most likely originates from a paired normal sample.  The bam file or files should be obtained with a long read technology and could originate for tumor, normal, or tumor+normal.

```
./linker phase -i ./input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
./linker phase -i ./input.bam -s ./second_input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
```
  * Output is haplotype solution file: haplotype_solution.dat

The output of this command is a file which contains the minimum energy solution to the germline haplotype.  More information on this file is described in the I/O section below.
-->

<!--
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
-->

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
