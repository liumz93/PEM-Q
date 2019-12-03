#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if( length(args)<1 ){
	print("\nAttention !!!    usage: Rscript plot_boxplot.R <tab file>     by Mengzhu")
	q()
}
tab <- args[1]

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

library("ggplot2")
data <- read.table(tab,header=TRUE,sep="\t")

# data <- data[which(data$Type == "Pair"),]
data <- data[which(data$Vector_inser_size >0),]
print(median(data$Vector_inser_size))

#### boxplot u6 ####

# data <- data[which(data$Type != "Discard"),]
data_summary <- function(x) {
   m <- mean(x)
   ymin <- min(x)
   ymax <- max(x)
   # ymin <- m-sd(x)
   # ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}
# data <- data[which(data$Type == "Small"),]
# data <- data[which(data$Vector_inser_size >0),]
# data <- data$Vector_inser_size
# data <- data$Qname
# write.table(data,"u6_bam_list",quote=FALSE,col.names = TRUE,row.names = FALSE)
# print(data_summary(data))
p <- ggplot(data, aes(x="", y=log10(Vector_inser_size))) +
  geom_jitter(color="black",alpha=0.5) + geom_boxplot(color="black",alpha=0,size=0.7)+ theme_bw() +
  theme(text = element_text(size = 25),
  axis.text.x= element_text(color="black"),
  axis.text.y= element_text(color="black"),
  axis.title.x= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)),#title distance to axis
  axis.title.y= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))#title distance to axis
p
ggsave("len.pdf",plot=p, width = 5, height = 5)