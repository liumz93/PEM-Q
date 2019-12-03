#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if( length(args)<4 ){
	print("usage: Rscript Cor.R Name sampleA sampleB sampleC (txtfile) Write by Mengzhu")
	q()
}

name <- args[1]
sample1 <- args[2]
sample2 <- args[3]
sample3 <- args[4]

n1 <- substring(sample1, 1, 6)
n2 <- substring(sample2, 1, 6)
n3 <- substring(sample3, 1, 6)

A <- read.table(sample1, sep="\t", header = FALSE)
B <- read.table(sample2, sep="\t", header = FALSE)
C <- read.table(sample3, sep="\t", header = FALSE)

A <- A$V4
B <- B$V4
C <- C$V4

codata <- data.frame(A,B,C)
# colnames(codata) <- c(n1,n2,n3)

TXT <- cor(codata)
print(TXT)

write.csv(cor(codata),file=paste(name,"_cor.txt",sep=""),quote=FALSE)


# write.table(txt,paste(name,"_cor.txt",sep=""),sep="\t",quote=FALSE)