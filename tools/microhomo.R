#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if( length(args)<2 ){
	print("usage: Rscript microhomo.R file name (Write by Mengzhu)")
	q()
}

miho <- function(data,type){
	miho <- data$MH_end
	miho_pat <- as.data.frame(table(miho))
	miho_pat <- miho_pat[order(-miho_pat$Freq),]
	write.table(miho_pat, paste(name,"_",type,"_homo_pat.txt",sep=""), sep="\t", quote=FALSE, row.names = FALSE)
}

#read in file
file <- args[1]
name <- args[2]
data <- read.table(file, sep="\t", header = TRUE)

# #~~microhomology analysis~~#
#total
miho(data,"All")
# #break
# datb <- data[(data$Rname == chr) & (data$Junction >= cutsite-500000) & (data$Junction <= cutsite+500000),]
# miho(datb,"Brk")
# #genomewide
# datb <- data[data$Rname != chr | (data$Rname == chr & data$Junction < cutsite-500000) | (data$Rname == chr & data$Junction > cutsite+500000),]
# miho(datb,"GenW")


