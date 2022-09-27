#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "output","character", "file path to plot to"
  )
  
  OPTS <- c(
    "binfile","character","","write bin info to file",
    "binsize","integer",2000000,"bps per bin",
    "assembly","character","mm9","genome assembly",
    "strand","integer",0,"1 for positive strand, -1 for minus strand, 0 for both strands, 2 for combined strands",    
    "brkchr","character","","",
    "brksite","integer",0,"",
    "brkstrand","integer",0,"",
    "featurefile","character","","e.g. RefGene",
    "chr","character","","",
    "rstart","integer",0,"start",
    "rend","integer",0,"end",
    "rmid","integer",0,"",
    "rwindow","integer",0,"",
    "binnum","integer",0,"",
    "plottype","character","dot","",
    "plotshape","character","arrow","",
    "showM","integer",0,"",
    "showY","integer",0,"",
    "ymax","integer",0,""
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("TranslocHelper.R")
  parseArgs("TranslocPlot.R", ARGS, OPTS)
  
} else {
  source("~/TranslocPipeline/R/Rsub.R")
  source("~/TranslocPipeline/R/TranslocHelper.R")
  tlxfile <- "/Volumes//AltLab/Translocation//RawData/Alt024-20130429/NewPipelineTest/results-full/CC004_Alt024/CC004_Alt024.tlx"
  output <- "~/Working/TranslocTesting/TranslocPlot.pdf"
  binfile <- "~/Working/TranslocTesting/TranslocPlot_bins.txt"
  binsize <- 20000000
  assembly <- "mm9"
  featurefile <- "~/Desktop/mm9-TCRlocifeaturefilea.txt"
  chr <- "chr14"
  strand <- 0
  rstart <- 0
  rend <- 0
  rmid <- 54733143
  rwindow <- 200
  binnum <- 200
  showM <- 0
  showY <- 0
  plottype <- "linear"
  plotshape <- "arrow"
  ymax <- 0
  brkchr <- "chr15"
  brksite <- 61818880
  brkstrand <- 1
}

#argument checking
if (chr == "") plottype <- "dot"
if (strand == 2) plotshape <- "octogon"


suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))

# denom <- c(1,5,20,200,2000,20000,20000)
# pal <- brewer.pal(9,"Set1")
# pal <- c("black",pal[c(1,5,2,3,4,7)])

denom <- c(1,5,20,100,500,2000)
# pal <- brewer.pal(7,"Set1")
# pal <- c("black",pal[1],rgb(255,162,14,maxColorValue=255),pal[c(3,2,4)])
# pal <- c("black","light blue","pink","orange","red","purple")
pal <- brewer.pal(9,"Set1")
pal <- c(gray(3/12),topo.colors(12)[c(3)],terrain.colors(12)[c(1,6)],pal[c(5,1)])
	
# Read in cytogenetic band data
cyto <- getCytoBands(assembly)
features <- getFeatures(assembly,featurefile)

# Read in chomrosome length data
chrlen <- getChromLens(assembly)

if (chr != "") {
  if (! chr %in% names(chrlen)) stop("Error: chromosome ",chr," not found")
  chrlen <- chrlen[chr]
} else {
  if (!showM) chrlen <- chrlen[names(chrlen)!="chrM"]
  if (!showY) chrlen <- chrlen[names(chrlen)!="chrY"]
}

gr <- createGenomicRanges(chrlen,rstart=rstart,rend=rend,rmid=rmid,rwindow=rwindow,binsize=binsize,binnum=binnum)
binsize <- end(gr)[length(gr)] - start(gr)[length(gr)] + 1



columnsToRead <- c("Rname","Junction","Strand")
tlx <- readTLX(tlxfile,columnsToRead)



tlxtot <- nrow(tlx)

if (tlxtot < 1) quit()

tlx <- tlx[tlx$Rname %in% names(chrlen),]
if (strand == 1 || strand == -1) {
  tlx <- tlx[tlx$Strand == strand,]
} else if (strand == 2) {
  tlx$Strand <- 1
}

tlxgr <- tlxToGR(tlx,chrlen)

gr$hits <- countOverlaps(gr,tlxgr)

gr$mids <- end(gr) - binsize/2

tlxdisp <- sum(gr$hits)

tmphits <- gr$hits
gr$hitvec <- list(c())
for (i in length(denom):1) {
  gr$hitvec <- mapply(function(x,y) { if (y > 0) c(x,rep(i,y)) else x },gr$hitvec,tmphits%/%denom[i])
  tmphits <- tmphits - tmphits%/%denom[i]*denom[i]
}
gr$hitveclen <- unlist(lapply(gr$hitvec,length))


pdf(output,width=11,height=8.5)

marginsize <- 0.5
marginunits <- "inches"

pageVP <- viewport(name="page",width=unit(1,"npc")-unit(2*marginsize,marginunits),height=unit(1,"npc")-unit(2*marginsize,marginunits))
pushViewport(pageVP)

headerheight <- 3
headerunits <- "lines"

headerVP <- viewport(y=unit(1,"npc"),name="header",width=unit(1,"npc"),height=unit(headerheight,headerunits),just="top")
pushViewport(headerVP)

printHeader(tlxfile,tlxdisp,tlxtot,assembly,chr,denom,pal,plottype,plotshape)

popViewport()
plotVP <- viewport(name="plot",y=0,height=unit(1,"npc")-unit(headerheight,headerunits),just="bottom")
pushViewport(plotVP)


if (length(chrlen) > 1) {

  rotateVP <- 1
  
  chrwidth <- 1
  chrwidthunit <- "mm"
  
  chrpos <- rep(0,length(chrlen))
  names(chrpos) <- names(chrlen)
  gr$ypos <- 0
  for (i in 1:length(chrlen)) {
    gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="-"]$ypos <- rev(1:length(gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="-"]))
    gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="+"]$ypos <- rev(1:length(gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="+"]))
    
    
    if (i == 1) {
      chrpos[i] <- max(gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="-"]$hitveclen)
      
    } else {
      
      negbins <- rev(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "-"])
      
      posbins <- rev(gr[seqnames(gr) == names(chrlen)[i-1] & strand(gr) == "+"])
      
      bothbins <- min(length(posbins),length(negbins))
      
#       negbins <- negbins[negbins$ypos <=  min(max(negbins$ypos),max(posbins$ypos))]
#       posbins <- posbins[posbins$ypos <=  min(max(negbins$ypos),max(posbins$ypos))]
      
      chrpos[i] <- chrpos[i-1] + max(negbins$hitveclen[1:bothbins] + posbins$hitveclen[1:bothbins],
                                     negbins$hitveclen[1:min(bothbins,length(posbins)-1)]+posbins$hitveclen[2:min(bothbins+1,length(posbins))],
                                     negbins$hitveclen[2:min(bothbins+1,length(negbins))]+posbins$hitveclen[1:min(bothbins,length(negbins)-1)]) + 1
    }
  }
  chrposlim <- chrpos[length(chrpos)] + max(gr[seqnames(gr) == names(chrlen)[length(chrlen)] & strand(gr) == "+"]$hitveclen)
  
  totalwidthmm <- convertX(unit(1,"npc"),"mm",valueOnly=T)
  chrwidthnative <- chrwidth*chrposlim/(totalwidthmm-length(chrlen)*chrwidth)
  chrposlim <- chrposlim+chrwidthnative*length(chrlen)
  
  genomeVP <- viewport(name="genome",y=unit(2,"lines"),height=unit(1,"npc")-unit(2,"lines"),just="bottom",xscale=c(0,chrposlim),yscale=c(1,max(chrlen)))
  pushViewport(genomeVP)
  
  chrpos <- 0:(length(chrlen)-1)*chrwidthnative + chrpos
  
  grid.rect(x=unit(chrpos,"native"),y=unit(0,"npc"),width=unit(chrwidth,chrwidthunit),height=unit(chrlen,"native"),just=c("left","bottom"))
  cyto$Xpos <-chrpos[match(cyto$Chr,names(chrpos))]
  grid.rect(x=unit(cyto$Xpos,"native"),y=unit(chrlen[cyto$Chr] - cyto$Start,"native"),width=unit(chrwidth,chrwidthunit),height=unit(cyto$End-cyto$Start,"native"),just=c("left","top"),gp=gpar(fill=cyto$Color,lty=0))
  
  
  if (brkchr %in% names(chrpos) && brksite > 0) {
    arrowlen <- 5
    arrowlenunit <- "mm"
    arrowlennative <- convertHeight(unit(arrowlen,arrowlenunit),"native",valueOnly=T)
    
    
    grid.rect(x=unit(chrpos[brkchr],"native"),y=unit(chrlen[brkchr]-brksite,"native"),width=unit(chrwidth,chrwidthunit),just="left",height=unit(1,"mm"),gp=gpar(fill="yellow",linejoin="mitre"))
    if (brkstrand == 1 || brkstrand == -1) {
      ypoints <- unit(chrlen[brkchr]-c(brksite-brkstrand*arrowlennative,brksite-brkstrand*arrowlennative/2,brksite-brkstrand*arrowlennative/2,brksite,brksite-brkstrand*arrowlennative/2,brksite-brkstrand*arrowlennative/2,brksite-brkstrand*arrowlennative),"native")
      xpoints <- unit(c(chrpos[brkchr]+chrwidthnative*3/4,chrpos[brkchr]+chrwidthnative*3/4,chrpos[brkchr]+chrwidthnative,chrpos[brkchr]+chrwidthnative/2,chrpos[brkchr],chrpos[brkchr]+chrwidthnative*1/4,chrpos[brkchr]+chrwidthnative*1/4),"native") 
      grid.polygon(x=xpoints,y=ypoints,gp=gpar(fill="yellow",linejoin="mitre"))
    }
  }
  
  chrVPs <- list()
  
  for (i in 1:length(chrlen)) {
    negVP <- viewport(x=unit(chrpos[i],"native"), y=unit(0,"npc"), 
                      width=unit(max(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "-"]$hitveclen),"native"),
                      height=unit(chrlen[i],"native"), just=c("right","bottom"),
                      xscale=c(1,0),yscale=c(chrlen[i],1),clip="off")
    posVP <- viewport(x=unit(chrpos[i]+chrwidthnative,"native"), y=unit(0,"npc"), 
                      width=unit(max(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "+"]$hitveclen),"native"),
                      height=unit(chrlen[i],"native"), just=c("left","bottom"),
                      xscale=c(0,1),yscale=c(chrlen[i],1),clip="off")
    
    
    chrVPs[[i]] <- list(posVP,negVP)
  }
  pushViewport(viewport(name="chrlabel",y=unit(0,"npc"),height=unit(2,"lines"),just="top",xscale=c(0,chrposlim)))
  grid.text(sub("chr","",names(chrlen)),x=unit(chrpos+0.75/2,"native"),gp=gpar(cex=1.25))
  popViewport()
  
} else {
  
  rotateVP <- 0
  
  rstart <- min(start(gr))
  rend <- max(end(gr))
  
  if (plottype == "dot") {
    yscalewidth = unit(0,"lines")
  } else {
    yscalewidth = unit(2,"lines")
  }
  
  xscaleVP <- viewport(name="xscale",x=yscalewidth,width=unit(1,"npc")-yscalewidth,y=unit(1,"npc"),height=unit(2,"lines"),just=c("left","top"),xscale=c(rstart,rend),clip="off")
  pushViewport(xscaleVP)
  plotXScale(chr,rstart,rend)
  popViewport()
  
 
  
  genomeVP <- viewport(name="genome",x=yscalewidth,width=unit(1,"npc")-yscalewidth,height=unit(1,"npc")-unit(4,"lines"),just=c("left"),xscale=c(rstart,rend),yscale=c(-1,1))
  pushViewport(genomeVP)
  
  chrwidth <- 3
  chrwidthunit <- "mm"
  
  grid.rect(y=unit(0,"native"),height=unit(chrwidth,chrwidthunit))  
    
  if (rend - rstart < 2000000 || featurefile != "") {
    features <- subset(features,Chr == chr & End >= rstart & Start <= rend)
    features <- features[!duplicated(features$Name),]
    features <- features[with(features,order(Start)),]
    features$Start <- ifelse(features$Start < rstart, rstart, features$Start)
    features$End <- ifelse(features$End > rend, rend, features$End)
    if (nrow(features) > 0) {
      grid.rect(x=unit(features$Start,"native"),y=unit(0,"native"),width=unit(features$End-features$Start,"native"),height=unit(chrwidth,chrwidthunit),just="left",gp=gpar(fill=getCytoColor()["gpos25"]))
      featureVP <- viewport(name="feature",y=unit(0,"npc"),height=unit(2,"lines"),just="top",xscale=c(rstart,rend),clip="off")
      pushViewport(featureVP)
      plotFeatures(features,chr,rstart,rend)
		  #,col="gray")
      popViewport()
    }
  } else {
    cyto <- subset(cyto, Chr == chr & End >= rstart & Start <= rend)
    cyto$Start <- ifelse(cyto$Start < rstart, rstart, cyto$Start)
    cyto$End <- ifelse(cyto$End > rend, rend, cyto$End)
    grid.rect(x=unit(cyto$Start,"native"),y=unit(0,"native"),width=unit(cyto$End-cyto$Start,"native"),height=unit(chrwidth,chrwidthunit),just="left",gp=gpar(fill=cyto$Color,lty=0))
  }
  
  
  
  if (names(chrlen)[1] == brkchr && brksite >= rstart && brksite <= rend) {
    grid.lines(x=unit(brksite,"native"),y=c(unit(0.05,"npc"),unit(0.95,"npc")),gp=gpar(lty=3))
    grid.rect(x=unit(brksite,"native"),y=unit(0,"native"),width=unit(1,"mm"),height=unit(chrwidth,chrwidthunit),gp=gpar(fill="yellow",linejoin="mitre"))
    if (brkstrand == 1 || brkstrand == -1) {
      arrowlen <- 5
      arrowlenunit <- "mm"
      arrowlennative <- convertWidth(unit(arrowlen,arrowlenunit),"native",valueOnly=T)
      

      xpoints <- unit(c(brksite-brkstrand*arrowlennative,brksite-brkstrand*arrowlennative/2,brksite-brkstrand*arrowlennative/2,brksite,brksite-brkstrand*arrowlennative/2,brksite-brkstrand*arrowlennative/2,brksite-brkstrand*arrowlennative),"native")
      ypoints <- unit(0,"native")+unit(c(chrwidth/4,chrwidth/4,chrwidth/2,0,-chrwidth/2,-chrwidth/4,-chrwidth/4),chrwidthunit)
      
      grid.polygon(x=xpoints,y=ypoints,gp=gpar(fill="yellow",linejoin="mitre"))
    }
  }
  
  ymin <- 0
  if (plottype == "dot") {
    ymax <- 1
  } else {
    if (ymax == 0) ymax <- max(gr$hits)
    if (plottype == "log") {
      ymax <- log10(ymax)
      ymin <- -1
    }
  }
    
  posVP <- viewport(y=unit(0,"native")+unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just="bottom",xscale=c(rstart,rend),yscale=c(ymin,ymax),clip="on")
  negVP <- viewport(y=unit(0,"native")-unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just="top",xscale=c(rstart,rend),yscale=c(ymax,ymin),clip="on")
  chrVPs <- list(list(posVP,negVP))
  
  if (plottype != "dot") {
    pushViewport(viewport(x=unit(-2,"mm"),y=unit(0,"native")+unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just=c("left","bottom"),yscale=c(ymin,ymax),clip="off"))
    grid.yaxis(gp=gpar(cex=0.75))
    popViewport()
    
    pushViewport(viewport(x=unit(-2,"mm"),y=unit(0,"native")-unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just=c("left","top"),yscale=c(ymax,ymin),clip="off"))
    grid.yaxis(gp=gpar(cex=0.75))
    popViewport()
  }
  
}

for (i in 1:length(chrlen)) {
  posVP <- chrVPs[[i]][[1]]
  negVP <- chrVPs[[i]][[2]]
  
  
  
  pushViewport(posVP)
  plotJunctions(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "+",],binsize,strand=1,plottype=plottype,plotshape=plotshape,pal=pal,rotateVP=rotateVP)
  popViewport()
  
  pushViewport(negVP)
  plotJunctions(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "-",],binsize,strand=-1,plottype=plottype,plotshape=plotshape,pal=pal,rotateVP=rotateVP)
  popViewport()
  
}
  


dev.off()

if (binfile != "") {
  if (strand == 2) {
    bin_output <- data.frame(Rname=as.vector(seqnames(gr)),Rstart=start(gr),Rend=end(gr),Hits=gr$hits)
  } else {
    bin_output <- data.frame(Rname=as.vector(seqnames(gr)),Rstart=start(gr),Rend=end(gr),Strand=as.vector(strand(gr)),Hits=gr$hits)
  }
  
  write.table(bin_output,binfile,quote=F,sep="\t",na="",row.names=F,col.names=T)
}

