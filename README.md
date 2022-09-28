# PEM-Q

PEM-Q is a brand new pipeline built for analysis  of PEM-seq (https://www.nature.com/articles/s41421-019-0088-8/). Comparing with superQ (https://github.com/liumz93/superQ), PEM-seq is a more powerful tool for analyzing repair outcomes of genome editing, including indels, microhomologies, large deletions, vector integrations, and translocations.

## Getting Started

### PEM-Q structure

1.Alignment

2.Random molecular barcodes extraction

3.Classify translocations

4.Classify indels

5.Events dedup

6.Statistics

### Prerequisites

Please be sure to install: 

#### Python packages

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

#### Other tools

1. FLASH v1.2.11
2. bwa-0.7.12-r1034
3. samtools 1.3.1

#### Build bwa index

Please add your bwa index into the environment configuration file (~/.bashrc).

```
export BWA_INDEX=YOUR/BWA_INDEX/PATH
```

For example:

```
export BWA_INDEX=/home/mengzhu/database/bwa_indexes
```

And then in $BWA_INDEX, generate a genome folder and build your bwa index. For example:

```
YOUR@SERVER:YOUR/BWA_INDEX/PATH$ mkdir mm10
YOUR@SERVER:YOUR/BWA_INDEX/PATH$ cd mm10
YOUR@SERVER:YOUR/BWA_INDEX/PATH$ bwa index -a bwtsw -p mm10 mm10.fa
YOUR@SERVER:YOUR/BWA_INDEX/PATH$ ls /home/mengzhu/database/bwa_indexes/mm10
mm10.amb  mm10.ann  mm10.bwt  mm10.pac  mm10.sa
```

### Installing PEM-Q

```
git clone https://github.com/liumz93/PEM-Q
```
Add below variables into the environment configuration file (~/.bashrc) for PEM-Q according to your  installation path.

```
export PATH="/your/path/PEM-Q:$PATH"
export PATH="/your/path/PEM-Q/main:$PATH"
export PATH="/your/path/PEM-Q/tools:$PATH"
```
## Transloc_pipeline

We have applied the TranslocPlot.R and TranslocHTMLReads.pl from transloc_pipeline (https://github.com/robinmeyers/transloc_pipeline) fast visualization, so we suggest you download the relevant scripts from their R and lib folders. Noted that if no download does not affect the running of PEM-Q and generation of the final files although errors will reporte in the last step.

For more flexible visualization, please convert result file into bdg format and view them in IGV (https://software.broadinstitute.org/software/igv/download).

## Running PEM-Q

### Basic analysis

```
PEM-Q.py genome test_sample cut-site target_chromosome primer_start primer_end primer_strand primer_sequence
```
For test:

```
PEM-Q.py mm10 CC055c 61986726 chr15 61986633 61986652 + GGAAACCAGAGGGAATCCTC

```
### Vector integration analysis
```
vector_analyze.py test_sample vector.fa genmome target_chromosome primer_start sgRNA_start sgRNA_end
```
For test:

```
vector_analyze.py CC055c Spcas9_pX330.fa mm10 chr15 + 7937 7956
```


## Output

Final results and statistic files were put into the results folder.

### statistics file

a statistic file include numbers of all editing events and provide editing efficiency. For example:

```
NoJunction	638686
Deletion	289463
Insertion	117445
Intra_Translocation.del	268
Close_inversion.del	8680
Intra_Translocation.inver	636
Inter_Translocation	11117
Translocation	12021
Editing Events	418929
Total Events	1066295
Editing Efficiency	0.39288283261198825
```
### Events file

There are separate files for deletions, insertions, inversions and intra- or inter chromosomal translocations. Each row represent a edit event. For example:

```
Qname	Bait_rname	Bait_strand	Bait_start	Bait_end	Prey_rname	Prey_strand	Prey_start	Prey_end	Rname	Strand	Junction	Sequence	B_Qstart	B_Qend	Qstart	Qend	Qlen	Insertion	Microhomolog	Prey_MQ	Barcode
ST-E00578:442:HF32JCCX2:7:2205:3254:67111	chr15	+	61986633	61986727	chr1	-	3432216	3432362	chr1	-3432362	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTGCCACTCCACTTACATAGTTGCTAAGTTGTTTGTTATACTGTACATATGTATGTGCCCATGAGTGCATGTGTATACATTTAAATTTCATATTGAAGCTTTAAATTTTGATTATTCATTCAAGATTTAGACTTAGTAGACATAAAGGAGCCACGCGTGCTCTACACGTTTATCAACGTCGT	5	99	100	246	278			60	CGTTTATCAACGTCGT	
ST-E00578:442:HF32JCCX2:7:1102:24982:42200	chr15	+	61986633	61986726	chr1	-	3432216	3432364	chr1	-3432364	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTTCCACTCCACTTACATAGTTGCTAAGTTGTTTGTTATACTGTACATATGTATGTGCCCATGAGTGCATGTGTATACATTTAAATTTCATATTGAAGCTTTAAATTTTGATTATTCATTCAAGATTTAGACTTAGTAGACATAAAGGAGCCACGCGTGCTCTACACGTTTATCAACGTCGTG	5	98	98	246	279		T	60	CGTTTATCAACGTCGTG	
ST-E00578:442:HF32JCCX2:7:2110:2727:13668	chr15	+	61986633	61986726	chr1	-	3432312	3432364	chr1	-3432364	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTTCCACTCCACTTACATAGTTGCTAAGTTGTTTGTTATACTGTACATATGTAT	5	98	98	150	150		T	60	CGTTTATCAACGTTGTG	
ST-E00578:442:HF32JCCX2:7:2218:1976:4807	chr15	+	61986633	61986726	chr1	+	4434289	4434342	chr1	4434289	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTACATATGCATGGTATATATATATATGTACATCCAGGCAAACATTCATACACA	5	98	97	150	150		CT	60	CGTATAACAGCATCGAA	
ST-E00578:442:HF32JCCX2:7:2201:26017:68148	chr15	+	61986633	61986717	chr1	+	6168405	6168467	chr1	6168405	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGAAGATTCTGGTCTGTGGTGTTCTTACTGGCCGGTCGTGAGAACGCGGCTAATAACAATTGG	5	89	88	150	150		AG	0	TTAATCTCACGATACGA	
ST-E00578:442:HF32JCCX2:7:1106:3346:40688	chr15	+	61986633	61986726	chr1	+	6168640	6168706	chr1	6168640	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTATAGCTTTACAAGGTACGCCTGGCCTTGAACTTTCTAACGAAATTCAGGACAGTCTATCAGAAGTACCACGCGTGCTCTACACACTTTCTAGGTTCGAA	5	98	98	164	197		T	0	CACTTTCTAGGTTCGAA	
ST-E00578:442:HF32JCCX2:7:1210:11221:59112	chr15	+	61986633	61986726	chr1	+	6168640	6168706	chr1	6168640	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTATAGCTTTACAAGGTACGCCTGGCCTTGAACTTTCTAACGAAATTCAGGACAGTCTATCAGAAGTACCACGCGTGCTCTACACCCTTTCTAGGTTCGAA	5	98	98	164	197		T	0	CCCTTTCTAGGTTCGAA	
ST-E00578:442:HF32JCCX2:7:1104:29883:50551	chr15	+	61986633	61986726	chr1	+	6168640	6168751	chr1	6168640	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTATAGCTTTACAAGGTACGCCTGGCCTTGAACTTTCTAACGAAATTCAGGACAGTCTATCAGAAGTAAAGTGGAAAATGGCTTTACGAGGTATGCTTGGCCTTAAACTTTCTACCACGCGTGCTCTACACTCGTGTAAGATTCCCT	5	98	98	209	243		T	0	CTCGTGTAAGATTCCCT	
ST-E00578:442:HF32JCCX2:7:2102:25327:37383	chr15	+	61986633	61986726	chr1	+	6168640	6168706	chr1	6168640	AGGAGGAAACCAGAGGGAATCCTCACATTCCTACTTGGGATCCGCGGGTATCCCTCGCGCCCCTGAATTGCTAGGAAGACTGCGGTGAGTCGTGATCTATAGCTTTACAAGGTACGCCTGGCCTTGAACTTTCTAACGAAATTCAGGACAGTCTATCAGAAGTACCACGCGTGCTCTACACTGTTTGTACACTAAGA	5	98	98	164	197		T	0	CTGTTTGTACACTAAGA
```
###Explanation of column titles

**Qname**: sequence name

**Bait_rname**: chromosome of bait

**Bait_strand**: strand of bait

**Bait_start**: chromosomal start position of bait

**Bait_end**: chromosomal end position of bait

**Prey_rname**: chromosome of prey

**Prey_strand**: strand of prey

**Prey_start**: chromosomal start position of prey

**Prey_end**: chromosomal end position of prey

**Rname**: chromosome of junction

**Strand**: strand of junction

**Junction**: chromosomal position of junction

**Sequence**: raw read sequence

**Insertion**: insertion sequence

**Microhomolog**: microhomology sequence used by deletions or transloctions

**Barcode**: random molecular barcode

### Useful tools to convert PEM-Q output into bdg format

For note, bdg format can be directly viewed in **igv**(https://igv.org/app/).

convert PEM-Q output into bdg by tab2bdg_PEMQ.py

```
tab2bdg_PEMQ.py CC055c_Translocation.tab mm10
```
convert vector output into bdg by vectorTab2bdg.py

```
vectorTab2bdg.py CC055c_all_vector.tab data/pX330_SpCas9.fa
```

