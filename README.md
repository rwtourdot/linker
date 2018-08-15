## Linker - Long Read and Linked Read sequencing analysis tools

Table of contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Commands](#commands)

Installation
------------

code requires bamtools, htslib, and c++11

Description
-----------

<img src="https://github.com/rwtourdot/linker/scripts/linker_flowchart.pdf" width=800/>

Commands
--------

```
## get help
./linker -h

## phase long reads
./linker phase -h

## phase with copy number
./linker cn_phase -h

## phase structural variants
./linker sv_phase -h

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
