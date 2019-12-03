#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if( length(args)<1 ){
	print("\nAttention !!!    usage: Rscript plot_boxplot.R <tab file>     by Mengzhu")
	q()
}
tab <- args[1]


library("ggplot2")
data <- read.table(tab,header=TRUE,sep="\t")

# data <- data[which(data$Type == "Pair"),]
data <- data[which(data$Vector_inser_size >0),]
print(median(data$Vector_inser_size))






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


#### boxplot ####

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
p <- ggplot(data, aes(x=Type, y=Vector_inser_size, fill=Type)) +
  geom_jitter(color="grey") + geom_boxplot(color="black",alpha=0,size=0.6)+ theme_bw() + 
  theme(text = element_text(size = 25),
  axis.text.x= element_text(color="black"),
  axis.text.y= element_text(color="black"),
  axis.title.x= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)),#title distance to axis
  axis.title.y= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))#title distance to axis
p 
ggsave("len.pdf",plot=p)






# #### boxplot of lenti####
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
#   geom_jitter(color="grey") + theme_bw() +
#   theme(text = element_text(size = 25),
#   axis.text.x= element_text(color="black"),
#   axis.text.y= element_text(color="black"),
#   axis.title.x= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)),#title distance to axis
#   axis.title.y= element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))#title distance to axis
# p
# ggsave("len_lenti.pdf",plot=p)
#
