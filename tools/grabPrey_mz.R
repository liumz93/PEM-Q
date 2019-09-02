#!/usr/bin/env Rscript
#Author:Liu Mengzhu

args <- commandArgs(trailingOnly = TRUE)
if( length(args)<3 ){
	print("usage: Rscript grabPrey.R <bed file> <tlx file> <output name>   (write by Mengzhu)")
	q()
}

print("Welcome use tools written by Mengzhu:)")
bed_file <- args[1]
input_file <- args[2]
output_file <- args[3]

# read in tlx file
data <- read.table(file=input_file, header=TRUE, sep="\t")
# read in bed file
bfile <- read.table(file=bed_file, header=FALSE, sep="\t")

#each chr tlx
chr <- unique(data$Rname)
datb <- data.frame()
for( i in 1:length(chr) ){
	assign(paste("data_", chr[i], sep=""), subset(data, data$Rname == chr[i]))
}

count <- 0
for(i in 1:nrow(bfile) ){
	chromosome <- bfile[i,1]
	ss <- bfile[i,2]
	es <- bfile[i,3]
	if (chromosome %in% chr == FALSE ) {next}
	#select preys
	assign(paste("select_",chromosome,sep=""),
	       subset(get(paste("data_",chromosome,sep="")), Junction >= ss & Junction < es))
	print(nrow(get(paste("select_",chromosome,sep=""))))
	datb<-rbind(datb,get(paste("select_", chromosome,sep="")))
}

# generate tlx data

# chromosome_r <- bfile[,1]
# datb <- data.frame()
# for( i in 1:length(chromosome_r) ){
# 	datb<-rbind(datb,get(paste("select_", chromosome_r[i],sep="")))
# }

datc <- subset(data,!(Qname %in% datb$Qname))

write.table(datb, file=paste(output_file,"_BED.tab",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(datc, file=paste(output_file,"_rmBED.tab",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)