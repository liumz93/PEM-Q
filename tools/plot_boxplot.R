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
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
# x and y axis are transformed and formatted
data <- read.table(tab,header=TRUE,sep="\t")

# data <- data[which(data$Type == "Pair"),]
data <- data[which(data$Vector_inser_size >0),]
print(median(data$Vector_inser_size))
print(max(data$Vector_inser_size))
print(min(data$Vector_inser_size))






# #### violin plot ####
# #颜色
# # scale_fill_brewer(palette="Blues")
# #scale_fill_brewer(palette="Dark2")
# #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", ))
#
# #change order of small/medium/large
# data$Type <- factor(data$Type, levels=c("Discard","Confident","Confident_2","Half","Suspected"), ordered=TRUE)
# data <- data[which(data$Type != "Discard"),]
# data_summary <- function(x) {
#    m <- mean(x)
#    ymin <- min(x)
#    ymax <- max(x)
#    # ymin <- m-sd(x)
#    # ymax <- m+sd(x)
#    return(c(y=m,ymin=ymin,ymax=ymax))
# }
# # data <- data[which(data$Type == "Small"),]
# # data <- data[which(data$Vector_inser_size >0),]
# # data <- data$Vector_inser_size
# # print(data_summary(data))
# p <- ggplot(data, aes(x=Type, y=log10(Vector_inser_size), fill=Type)) +
#   geom_violin() + stat_summary(fun.data=data_summary)+scale_fill_brewer(palette="Blues")
# ggsave("len_lenti.pdf",plot=p)


# #### boxplot 1 ####
#
# #change order of small/medium/large
# data$Type <- factor(data$Type, levels=c("Discard","Confident","Confident_2","Half","Suspected"), ordered=TRUE)
# data <- data[which(data$Type != "Discard"),]
# data_summary <- function(x) {
#    m <- mean(x)
#    ymin <- min(x)
#    ymax <- max(x)
#    # ymin <- m-sd(x)
#    # ymax <- m+sd(x)
#    return(c(y=m,ymin=ymin,ymax=ymax))
# }
# # data <- data[which(data$Type == "Small"),]
# # data <- data[which(data$Vector_inser_size >0),]
# # data <- data$Vector_inser_size
# # print(data_summary(data))
# p <- ggplot(data, aes(x=Type, y=log10(Vector_inser_size), fill=Type)) +
#   geom_jitter(color="grey") + geom_boxplot(color="black",alpha=0,size=0.6)+ theme_bw() +
#   theme(text = element_text(size = 25),
#   axis.text.x= element_text(color="black"),
#   axis.text.y= element_text(color="black"),
#   axis.title.x= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)),#title distance to axis
#   axis.title.y= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))#title distance to axis
# p
# ggsave("len.pdf",plot=p)



#### boxplot 2 ####
# color=brewer.pal(12,"Set3")[3:5],alpha=0.5
#change order of small/medium/large
data$Type <- factor(data$Type, levels=c("Discard","Confident","Confident_2","Half","Suspected"), ordered=TRUE)
data <- data[which(data$Type != "Discard"),]
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
# print(data_summary(data))
p <- ggplot(data, aes(x=Type,y=Vector_inser_size,color=Type)) + geom_jitter(aes(alpha=0)) + geom_boxplot(color="black",alpha=0,size=1) + theme_bw() +
scale_color_manual(values= c("palevioletred2","springgreen4","grey90")) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10,9000))
p
ggsave("len.pdf",plot=p, width = 4, height = 7)









# #### boxplot u6 ####
#
# # data <- data[which(data$Type != "Discard"),]
# data_summary <- function(x) {
#    m <- mean(x)
#    ymin <- min(x)
#    ymax <- max(x)
#    # ymin <- m-sd(x)
#    # ymax <- m+sd(x)
#    return(c(y=m,ymin=ymin,ymax=ymax))
# }
# # data <- data[which(data$Type == "Small"),]
# # data <- data[which(data$Vector_inser_size >0),]
# # data <- data$Vector_inser_size
# # data <- data$Qname
# # write.table(data,"u6_bam_list",quote=FALSE,col.names = TRUE,row.names = FALSE)
# # print(data_summary(data))
# p <- ggplot(data, aes(x="", y=log10(Vector_inser_size))) +
#   geom_jitter(color="darkorange3",alpha=0.5) + geom_boxplot(color="black",alpha=0,size=0.7)+ theme_bw() +
#   theme(text = element_text(size = 25),
#   axis.text.x= element_text(color="black"),
#   axis.text.y= element_text(color="black"),
#   axis.title.x= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)),#title distance to axis
#   axis.title.y= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))#title distance to axis
# p
# ggsave("len.pdf",plot=p, width = 5, height = 5)
