#!/usr/bin/env Rscript

library(plotly)

args <- commandArgs(trailingOnly = TRUE)
if( length(args)<2 ){
	print("\nAttention !!!    usage: Rscript LengthAnalysis.R <tab file> <sample>    by Mengzhu")
	q()
}

tlx <- args[1]
sample <- args[2]
# read in tlx file

#data <- read.table(file=input_file, header=TRUE, sep="\t")
data <- read.table(tlx, header=TRUE, sep="\t")
# data <- data[data$Type == "Suspected",]
# get bait and stitch length

Bl <- data$Bait_end-data$Bait_start+1
# Bl <- data$
# bait length 


bl <- Bl
a <- cut(bl,breaks=c(30:150))
c <- table(a)
b <- data.frame(c)
x <- c(31:150)
y <- b$Freq
y <- y/sum(y)
data<-data.frame(x,y)

print("Generating txt file for plotting!")
write.table(data,paste(sample,"_bait_length_fot_plot.txt",sep=""),sep="\t",quote=FALSE)

# print("Plotting bait length Frequecy lines!")
#
# p <- plot_ly(data, x = ~x) %>%
#         add_trace(y = ~y, name = 'bait length',type= 'scatter',mode = 'lines',line = list(color = I('blue'), width = 2)) %>%
# 		layout(title = "Bait Length",yaxis = list(title = "percentage"), xaxis = list(title = "length(bp)"))
# htmlwidgets::saveWidget(p, paste(sample,"_bait_length_graph.html",sep=""))