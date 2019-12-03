#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if( length(args)<2 ){
	print("\nAttention !!!    usage: Rscript LengthAnalysis.R <tab file> <sample>    by Mengzhu")
	q()
}

tab <- args[1]
sample <- args[2]
# read in tab file

# data <- read.table(file=input_file, header=TRUE, sep="\t")
data <- read.table(tab, header=TRUE, sep="\t")

# get bait and stitch length

# Bl <- nchar(as.character(data$Insertion))
Bl <- data$Bait_MH
print(head(Bl))

# length 


bl <- Bl
# bk <- c(0:300)
bk <- seq(0, 300, by=10)
a <- cut(bl,breaks=bk)
c <- table(a)
b <- data.frame(c)
x <- c(1:30)
y <- b$Freq
# y <- y/sum(y)
data<-data.frame(x,y)

print("Generating txt file for plotting!")
write.table(data,paste(sample,"_inser_len_bin10.txt",sep=""),sep="\t",quote=FALSE)

bl <- Bl
bk <- c(0:300)
a <- cut(bl,breaks=bk)
c <- table(a)
b <- data.frame(c)
x <- c(1:300)
y <- b$Freq
# y <- y/sum(y)
data<-data.frame(x,y)

print("Generating txt file for plotting!")
write.table(data,paste(sample,"_inser_len_bin1.txt",sep=""),sep="\t",quote=FALSE)