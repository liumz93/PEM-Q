#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if( length(args)<2 ){
	print("Attention !!!    usage: Rscript detectOffTarget.R <fa file> <pattern sequence> <samplename> by Mengzhu")
	q()
}

file <- args[1]
ps <- args[2]
name <- as.character(args[4])
mm <- 8
if(length(args) > 2){mm <- as.integer(args[3])}

.libPaths(c( .libPaths(), "/home/hulab/chen/R/x86_64-pc-linux-gnu-library/3.5") )
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(rtracklayer)

print("read in fa file ... ")
seqs <- readDNAStringSet(file)
# positive <- vmatchPattern("GCCTCTTTCCCACCCACCTTGGG", seqs)
# negative <- vmatchPattern(reverseComplement(DNAString("GCCTCTTTCCCACCCACCTTGGG")), seqs)
# print(positive)
# query <- DNAString("GCCTCTTTCCCACCCACCTTGGG")
query <- DNAString(ps)
max.mismatch <- mm

print("matching cut site sequence ...")
fwd <- vmatchPattern(query, seqs, max.mismatch=max.mismatch, fixed=FALSE)
print(fwd)
fwd <- as(fwd, "GRanges")
strand(fwd) <- "+"
# print(fwd)
rev <- vmatchPattern(reverseComplement(query), seqs, max.mismatch=max.mismatch, fixed=FALSE)
rev <- as(rev, "GRanges")
strand(rev) <- "-"
if (length(fwd) == 0){
	print( "Attension! No pattern found in fwd!!!")
}else{
	fwdfile <- paste(name,"_offtarget-fwd.gff",sep="")
	export(fwd, fwdfile)
	print("fwd file are generated successfully!")
}
if (length(rev) == 0){
	print("Attension! No pattern found in rev!!!")
}else{
	revfile <- paste(name,"_offtarget-rev.gff",sep="")
	export(rev, revfile)
	print("rev file are generated successfully!")
}
# # complete <- c(fwd,rev)
# print(rev)
