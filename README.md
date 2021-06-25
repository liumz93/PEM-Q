# PEM-Q

PEM-Q is a brand new pipeline built for analysis  of PEM-seq (). Comparing with superQ (https://github.com/liumz93/superQ), PEM-seq is a more powerful tool for analyzing repair outcomes of genome editing, including indels, microhomologies, large deletions, vector integrations, and translocations.

## Getting Started

PEM-Q consists of six parts:

1.Alignment

2.Random molecular barcodes extraction

3.Classify translocations

4.Classify indels

5.Events dedup

6.Statistics

### Prerequisites

Please be sure to install: 

####Python packages

1. os
2. sys
3. numpy
3. threading
4. time
5. docopt
6. pysam
7. re
8. pandas
9. Bio

####other tools

1. FLASH v1.2.11
2. bwa-0.7.12-r1034
3. samtools 1.3.1

### Installing

```
git clone https://github.com/liumz93/PEM-Q
```


## Running PEM-Q

### Basic analysis

```
PEM-Q.py genome test_sample cut-site target_chromosome primer_start primer_end primer_strand primer_sequence
```
example for test data at the c-Myc locus:

```
PEM-Q.py mm10 CC055c 61986726 chr15 61986633 61986652 + GGAAACCAGAGGGAATCCTC

```
### vector integration anlysis
```
vector_analyze.py test_sample vector.fa genmome target_chromosome primer_start sgRNA_start sgRNA_end
```
example for test:


```
vector_analyze.py test_sample Spcas9_pX330.fa mm10 chr15 + 7937 7956
```


## Output

Final results and statistic files were put into the results folder.
