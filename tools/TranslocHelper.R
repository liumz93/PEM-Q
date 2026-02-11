
split.by.chr <- function(gr) {
  split(gr,seqnames(gr))
}

myCountOverlaps <- function(bins, tlx.cumsum) {
  return(as.numeric(tlx.cumsum[end(bins)] - tlx.cumsum[pmax(start(bins)-1, 1)]))
}

calculate_local_significance <- function(bins, tlx.cumsum, bin.width, bg.width) {
    if (length(bins) < 1) return(numeric())
    bins.hits <- myCountOverlaps(bins, tlx.cumsum)
    bg.bins <- suppressWarnings(resize(bins, width=bg.width, fix="center") %>%
                                    trim())
    bg.bins.hits <- myCountOverlaps(bg.bins,tlx.cumsum)

    p.local <- maply(cbind(k=bins.hits,
                           N=bg.bins.hits,
                           w=pmax(width(bins), bin.width) / width(bg.bins)),
                     function(k, N, w) {
                         2 * sum(dbinom(k:N, N, w)) + (k/w - N - 1) * dbinom(k, N, w)
                     }, .expand=F)

    p.local <- ifelse(p.local < 0 | p.local > 1, 1, p.local)
    return(p.local)
}

calculate_chr_significance <- function(bins, tlx.cumsum, bin.width) {
    if (length(bins) < 1) return(numeric())

    bins.hits <- myCountOverlaps(bins,tlx.cumsum)
    chr.hits <- as.numeric(tlx.cumsum[length(tlx.cumsum)])

    w.mode <- as.numeric(names(which.max(table(width(bins))))) /
        seqlengths(bins)[as.character(seqnames(bins)[1])]

    binom.calc.table <- dbinom(1:chr.hits, chr.hits, w.mode)
    binom.calc.cumsum <- rev(cumsum(rev(binom.calc.table)))
    scan.stat.table <- 2*binom.calc.cumsum +
        (1:chr.hits/w.mode - chr.hits - 1) * binom.calc.table

    p.chr <- maply(cbind(k=bins.hits,
                         N=chr.hits,
                         w=pmax(width(bins), bin.width) /
                             seqlengths(bins)[as.character(seqnames(bins))]),
                   function(k, N, w) {
                       if (w == w.mode) {
                           return( scan.stat.table[k] )
                       } else {
                           return(2 * sum(dbinom(k:N, N, w)) +
                                      (k/w - N - 1) * dbinom(k, N, w))
                       }
                   }, .expand=F)

    p.chr <- ifelse(p.chr < 0 | p.chr > 1, 1, p.chr)

    return(p.chr)
}


calculate.local.significance <- function(bins,tlx.cumsum,bin.width,bg.width,cores=4) {

  if (length(bins) < 1) return(bins)

#   print(as.character(seqnames(bins)[1]))

#   tlx.tree <- GIntervalTree(tlx.gr)

  bins$hits <- myCountOverlaps(bins,tlx.cumsum)
#   bins$hits <- countOverlaps(bins,tlx.tree)

  bg.bins <- suppressWarnings(resize(bins, width=bg.width, fix="center") %>%
                                  trim())
  bg.bins$hits <- myCountOverlaps(bg.bins,tlx.cumsum)
  bins$bg.local <- bg.bins$hits



  bins$p.local <- mapply(function(k,N,w) {

      return( 2* sum(dbinom(k:N,N,w)) + (k/w - N - 1) * dbinom(k,N,w))

  },bins$hits,
  bg.bins$hits,
  pmax(width(bins), bin.width) / width(bg.bins),
#   mc.cores=cores,
  SIMPLIFY=T)

  bins$p.local[bins$p.local < 0 | bins$p.local > 1] <- 1


  #   bg.bins$prob <- width(bg.bins)/(width(bg.bins)+bg.bins$hits)
  #   bins$p <- pnbinom(bins$hits-1,size=width(bins),prob=bg.bins$prob,lower.tail=F)
  return(bins)
}


calculate.chr.significance <- function(bins,tlx.cumsum,bin.width,cores=4) {

  if (length(bins) < 1) return(bins)

  bins$hits <- myCountOverlaps(bins,tlx.cumsum)


  chr.hits <- as.numeric(tlx.cumsum[length(tlx.cumsum)])
  bins$bg.chr <- chr.hits


  w.mode <- as.numeric(names(which.max(table(width(bins)))))/seqlengths(bins)[as.character(seqnames(bins)[1])]

  binom.calc.table <- dbinom(1:chr.hits,chr.hits,w.mode)
  binom.calc.cumsum <- rev(cumsum(rev(binom.calc.table)))
  scan.stat.table <- 2*binom.calc.cumsum + (1:chr.hits/w.mode - chr.hits - 1) * binom.calc.table

#   w.mode <- as.numeric(names(which.max(table(width(bins1)))))/seqlengths(bins1)[as.character(seqnames(bins1)[1])]
#   binom.calc.table <- rep(NA,nrow = max(bins$hits),ncol = max(chr.hit.table))



  bins$p.chr <- mapply(function(k,N,w) {
    if (w == w.mode) {
      return( scan.stat.table[k] )
    } else {
      return( 2* sum(dbinom(k:N,N,w)) + (k/w - N - 1) * dbinom(k,N,w))
    }
  },bins$hits,
  chr.hits,
  pmax(width(bins),bin.width)/seqlengths(bins)[as.character(seqnames(bins))],
#   mc.cores=cores,
  SIMPLIFY=T)

  bins$p.chr[bins$p.chr < 0 | bins$p.chr > 1] <- 1

  #   bg.bins$prob <- width(bg.bins)/(width(bg.bins)+bg.bins$hits)
  #   bins$p <- pnbinom(bins$hits-1,size=width(bins),prob=bg.bins$prob,lower.tail=F)
  return(bins)
}


trim_peaks <- function(peaks, tlx.gr) {


    if (length(peaks) < 1) return(peaks)

    start(peaks) <- start(tlx.gr[findOverlaps(peaks, tlx.gr, select="first")])
    end(peaks) <- end(tlx.gr[findOverlaps(peaks, tlx.gr, select="last")])

    return(peaks)
}



trimPeaks <- function(peaks, tlx.gr) {


  if (length(peaks) < 1) return(peaks)

  tlx.tree <- GIntervalTree(tlx.gr)



  start(peaks) <- start(tlx.tree[findOverlaps(peaks,tlx.tree,select="first")])
  end(peaks) <- end(tlx.tree[findOverlaps(peaks,tlx.tree,select="last")])



  return(peaks)
}


removeHotSpots <- function(tlxgr,hotspotfile) {
  hotspots <- import(hotspotfile)
  tlxgr <- tlxgr[tlxgr %outside% hotspots]
  return(tlxgr)
}

tlxToBW <- function(tlxgr,mapfile,chrlen,binsize,strand=0) {

  mapbw <- BigWigFile(normalizePath(mapfile))
  gr <- createGenomicRanges(chrlen,binsize=binsize,strand=strand)



  gr$hits <- countOverlaps(gr,tlxgr)

  gr$mapability <- unlist(lapply(split(gr,seqnames(gr)),function(gr,bwf) {
    print(paste("calculating mapability on",seqnames(gr)[1]))
    unlist(mclapply(split(gr,1:length(gr)),calculateMapability,bwf,mc.cores=cores))
  }, mapbw))

  gr <- gr[gr$mapability > 0]


  gr$dens <- gr$hits/(gr$mapability*(end(gr)-start(gr)+1))
  gr$adjHits <- gr$dens*binsize
  gr$score <- gr$adjHits

  return(gr)
}



readTLX <- function(tlxfile,columnsToRead = c(),columnsToIgnore = c()) {
  if (length(columnsToRead) > 0) {
    header <- readHeader(tlxfile)
    colClasses <- proColumnClasses(header,columnsToRead)
    tlx <- read.delim(tlxfile,header=T,colClasses=colClasses,as.is=T)
  } else if (length(columnsToIgnore) > 0 ) {
    header <- readHeader(tlxfile)
    colClasses <- antiColumnClasses(header,columnsToIgnore)
    tlx <- read.delim(tlxfile,header=T,colClasses=colClasses,as.is=T)
  } else {
    tlx <- read.delim(tlxfile,header=T,as.is=T)
  }

  return(tlx)
}

tlxToGR <- function(tlx,chrlen,strand=0) {
  tlx <- tlx[tlx$Rname %in% names(chrlen),]

  if (strand == 1 || strand == -1) {
    tlx <- tlx[tlx$Strand == strand,]
  }
  gr <- GRanges(seqnames=tlx$Rname,ranges=IRanges(start=tlx$Junction,end=tlx$Junction),strand=ifelse(tlx$Strand==1,"+","-"),seqlengths=chrlen)
  return(gr)
}

calculateMapability <- function(range,bwf) {
  mapbw <- import(bwf,selection=range)
  if (length(mapbw) < 1) return(0)
  mapability <- sum((end(mapbw) - start(mapbw) + 1) * mapbw$score) / (tail(end(mapbw),n=1) - start(mapbw[1]))
  return(mapability)
}

plotJunctions <- function (gr,binsize,strand=1,plottype="dot",plotshape="arrow",pal=NULL,rotateVP=0) {

  if (sum(gr$hits) < 1) return()

  if (plottype == "dot") {
    dots <- data.frame(value=unlist(gr$hitvec),mids=rep(gr$mids,gr$hitveclen),stackpos=unlist(lapply(gr$hitveclen,function(x) { if (x>0) seq(1,x) })))
    if (rotateVP) {
      dotwidth <- convertHeight(-1*unit(binsize,"native"),"mm",valueOnly=T)
      vpheight <- convertWidth(unit(1,"npc"),"mm",valueOnly=T)
      dotheight <- convertWidth(strand*unit(min(dotwidth,vpheight/max(gr$hitveclen)),"mm"),"native",valueOnly=T)

    } else {
      dotwidth <- convertWidth(unit(binsize,"native"),"mm",valueOnly=T)
      vpheight <- convertHeight(unit(1,"npc"),"mm",valueOnly=T)
      dotheight <- convertHeight(strand*unit(min(dotwidth,vpheight/max(gr$hitveclen)),"mm"),"native",valueOnly=T)

    }

    if (plotshape == "arrow") {
      xpoints <- unit(c(dots$mids-strand*binsize/2,dots$mids-strand*binsize/5,dots$mids-strand*binsize/5,dots$mids+strand*binsize/2,dots$mids-strand*binsize/5,dots$mids-strand*binsize/5,dots$mids-strand*binsize/2),"native")
      ypoints <- unit(c(dots$stackpos-1/3,dots$stackpos-1/3,dots$stackpos,dots$stackpos-1/2,dots$stackpos-1,dots$stackpos-2/3,dots$stackpos-2/3)*dotheight,"native")
      vertices <- 7
    } else if (plotshape == "triangle") {
      xpoints <- unit(c(dots$mids-strand*binsize/2,dots$mids-strand*binsize/2,dots$mids+strand*binsize/2),"native")
      ypoints <- unit(c(dots$stackpos-1,dots$stackpos,dots$stackpos-0.5)*dotheight,"native")
      vertices <- 3
    } else if (plotshape == "octogon") {
      xpoints <- unit(c(dots$mids+strand*binsize/5,dots$mids+strand*binsize/2,dots$mids+strand*binsize/2,dots$mids+strand*binsize/5,
                        dots$mids-strand*binsize/5,dots$mids-strand*binsize/2,dots$mids-strand*binsize/2,dots$mids-strand*binsize/5),"native")
      ypoints <- unit(c(dots$stackpos,dots$stackpos-3/10,dots$stackpos-7/10,dots$stackpos-1,
                        dots$stackpos-1,dots$stackpos-7/10,dots$stackpos-3/10,dots$stackpos)*dotheight,"native")
      vertices <- 8
    }


    if (rotateVP) {
      tmp <- xpoints
      xpoints <- ypoints
      ypoints <- tmp
    }

    grid.polygon(x=xpoints,y=ypoints,id=rep(1:nrow(dots),vertices),gp=gpar(fill=pal[dots$value],lty=0))




  } else {
    xpoints <- gr$mids
    if (plottype == "linear") {
      ypoints <- gr$hits
    } else {
      ypoints <- log10(gr$hits)
      ypoints <- ifelse(ypoints == -Inf,-1,ypoints)
    }
    if (rotateVP) {
      tmp <- xpoints
      xpoints <- ypoints
      ypoints <- tmp
    }
    if (strand == 1) {
      linecolor <- pal[5]
    } else {
      linecolor <- pal[2]
    }
    grid.lines(x=unit(xpoints,"native"),y=unit(ypoints,"native"),gp=gpar(col=linecolor,lwd=2))
  }
}

plotXScale <- function(chr,rstart,rend) {

  sizeArray <- c(50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000)
  grid.xaxis(at=c(rstart,rend),label=c(" "," "))
  grid.text(label=formatBP(rend-rstart))
  grid.text(label=prettyNum(rstart,big.mark=","),x=unit(0,"npc"),just="left")
  grid.text(label=prettyNum(rend,big.mark=","),x=unit(1,"npc"),just="right")
}


plotFeatures <- function(features,chr,rstart,rend) {
  grid.text(features$Name,x=unit((features$Start+features$End)/2,"native"),y=unit(1-(((0:(nrow(features)-1))%%3))*0.33,"npc"),just="top",gp=gpar(cex=0.75))
#   grid.text(features$Name,x=unit((features$Start+features$End)/2,"native"),y=unit(0,"npc"),just="bottom",gp=gpar(cex=0.75))

}

printHeader <- function(tlxfile,tlxdisp,tlxtot,assembly,chr,denom,pal,plottype,plotshape) {

  if (plottype == "dot") {
    pushViewport(viewport(x=unit(1,"npc"),width=unit(1.5,"inches"),xscale=c(0,2),just="right"))
    x_legend <- unit(floor((1:length(denom)-1)/3) + 0.5,"native")
    y_legend <- unit((1:length(denom)-1)%%3 + 0.5,"lines")
    grid.text(label=denom,x=unit(x_legend,"native"),y=unit(y_legend,"lines"),just="left")

    if (plotshape == "arrow") {
      vertices <- 7
      xpoints <- rep(x_legend-unit(1,"mm"),each=vertices) - unit(c(4,8/3,8/3,0,8/3,8/3,4),"mm")
      ypoints <- rep(y_legend,each=vertices) + unit(c(4/5,4/5,2,0,-2,-4/5,-4/5),"mm")
    } else if (plotshape == "triangle") {
      vertices <- 3
      xpoints <- rep(x_legend-unit(1,"mm"),each=vertices) - unit(c(4,0,4),"mm")
      ypoints <- rep(y_legend,each=vertices) + unit(c(2,0,-2),"mm")
    } else if (plotshape == "octogon") {
      vertices <- 8
      xpoints <- rep(x_legend-unit(1,"mm"),each=vertices) - unit(c(3/10,0,0,3/10,7/10,1,1,7/10)*4,"mm")
      ypoints <- rep(y_legend,each=vertices) + unit(c(1/2,1/5,-1/5,-1/2,-1/2,-1/5,1/5,1/2)*4,"mm")
    }

    grid.polygon(x=xpoints,y=ypoints,id=rep(1:length(denom),each=vertices),gp=gpar(fill=pal,lty=0))
    popViewport()

    textwidth <- unit(1,"npc")-unit(1.5,"inches")

  } else {
    textwidth <- unit(1,"npc")
  }

  spacer <- unit(3,"mm")

  pushViewport(viewport(x=unit(0.5,"npc"),width=unit(0.5,"npc"),just="right"))
  titletext <- sub(".tlx","",basename(tlxfile))
  grid.text(titletext,x=unit(1,"npc")-spacer,just="right",gp=gpar(cex=min(2.5,1/convertHeight(stringHeight(titletext),"npc",valueOnly=T),(1-convertWidth(spacer,"npc",valueOnly=T))/convertWidth(stringWidth(titletext),"npc",valueOnly=T))))
  popViewport()

  pushViewport(viewport(x=unit(0.5,"npc"),width=textwidth,y=unit(1,"npc"),height=unit(0.5,"npc"),just=c("left","top")))
  hittext <- paste("Displaying",prettyNum(tlxdisp,big.mark=","),"of",prettyNum(tlxtot,big.mark=","),"hits")
  grid.text(hittext,x=spacer,just="left",gp=gpar(cex=min(1.25,1/convertHeight(stringHeight(hittext),"npc",valueOnly=T),(1-convertWidth(spacer,"npc",valueOnly=T))/convertWidth(stringWidth(hittext),"npc",valueOnly=T))))
  popViewport()

  pushViewport(viewport(x=unit(0.5,"npc"),width=unit(0.5,"npc"),y=unit(0,"npc"),height=unit(0.5,"npc"),just=c("left","bottom")))

  if (chr != "") {
    displaytext <- paste(assembly,"-",chr,"-",formatBP(binsize,1),"bins")
  } else {
    displaytext <- paste(assembly,"-",formatBP(binsize,1),"bins")
  }
  grid.text(displaytext,x=spacer,just="left",gp=gpar(cex=min(1.25,1/convertHeight(stringHeight(displaytext),"npc",valueOnly=T),(1-convertWidth(spacer,"npc",valueOnly=T))/convertWidth(stringWidth(displaytext),"npc",valueOnly=T))))
  popViewport()




#
#
#
#   bintext <- paste(formatBP(binsize,1),"bins")
#
#
#
#   titletext <- paste(sub(".tlx","",basename(tlxfile))," - Displaying ",prettyNum(tlxdisp,big.mark=",")," of ",prettyNum(tlxtot,big.mark=",")," hits - ",formatBP(binsize,1)," bins",sep="")
#   grid.text(titletext,gp=gpar(cex=min(2,1/convertHeight(stringHeight(titletext),"npc",valueOnly=T),1/convertWidth(stringWidth(titletext),"npc",valueOnly=T))))
}

createGenomicRanges <- function (chrlen,rstart=0,rend=0,rmid=0,rwindow=0,binsize=0,binnum=0,strand=0) {
  if (length(chrlen) > 1) {
    rends <- unlist(lapply(chrlen,function(x){rev(seq(from=x,to=1,by=-binsize,))}))
    rstarts <- unlist(lapply(rends,function(rend){ max(1,rend-binsize+1) } ))
    chrs <- Rle(names(chrlen),c(ceiling((chrlen)/binsize)))
  } else if (rmid != 0 && rwindow != 0) {
    if (binnum != 0) binsize <- ceiling(2*rwindow/binnum)
    rstarts <- c(rev(seq(from=rmid-binsize,to=rmid-rwindow,by=-binsize)),seq(from=rmid,to=rmid+rwindow-binsize,by=binsize))
    rends <- rstarts + binsize - 1
    chrs <- Rle(names(chrlen),length(rstarts))
  } else if (rstart != 0 && rend != 0 ) {
    if (binnum != 0) binsize <- ceiling((rend-rstart+1)/binnum)
    rstarts <- seq(from=rstart,to=rend-1,by=binsize)
    rends <- rstarts + binsize - 1
    chrs <- Rle(names(chrlen),length(rstarts))
  } else {
    rstart <- 1
    rend <- chrlen
    if (binnum != 0) binsize <- ceiling((rend-rstart+1)/binnum)
    rends <- rev(seq(from=rend,to=rstart,by=-binsize))
    rstarts <- unlist(lapply(rends,function(rend){max(1,rend-binsize+1)}))
    chrs <- Rle(names(chrlen),length(rstarts))
  }

  if (strand == 0) {
    strands <- Rle(c("+","-"),rep(length(chrs),2))
    gr <- GRanges(seqnames=rep(chrs,2),ranges=IRanges(start=rep(rstarts,2),end=rep(rends,2)),strand=strands,seqlengths=chrlen)
  } else {
    strands <- Rle("*",length(chrs))
    gr <- GRanges(seqnames=chrs,ranges=IRanges(start=rstarts,end=rends),strand=strands,seqlengths=chrlen)
  }
  seqlevels(gr) <- names(chrlen)
  return(gr)
}

getChromLens <- function (assembly) {
  chrlen <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/ChromInfo.txt',sep="/"),header=F,as.is=T,col.names=c('Name','Length'))
  # Convert data frame to vector with element names set as chromosome names
  chrlen <- structure(chrlen$Length,names=as.character(chrlen$Name))

  # Organize chromosome names into the correct order
  chrnum <- as.numeric(sub("chr","",names(chrlen[grep("chr[0-9]+",names(chrlen),perl=T)])))
  chrlet <- sub("chr","",names(chrlen[grep("chr[A-Z]+",names(chrlen),perl=T)]))
  chrnum <- paste("chr",chrnum[order(chrnum)],sep="")
  chrlet <- paste("chr",chrlet[match(c("X","Y","M"),chrlet)],sep="")
  chrlevels <- names(chrlen)[c(match(chrnum,names(chrlen)),match(chrlet,names(chrlen)))]
  # Order the chrlen object by this order
  chrlen <- chrlen[chrlevels]
  return(chrlen)
}

getCytoBands <- function (assembly) {
  cyto <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/cytoBand.txt',sep="/"),header=F,as.is=T,col.names=c('Chr','Start','End','ID','Type'))
  cytocolor <- getCytoColor()
  # Create new color column in cytoband data
  cyto$Color <- cytocolor[match(cyto$Type,names(cytocolor))]
  return(cyto)
}

getFeatures <- function (assembly,featurefile) {
  if (featurefile == "") featurefile <- paste(Sys.getenv('GENOME_DB'),assembly,'annotation/refGene.bed',sep="/")
  features <- read.delim(featurefile,header=F,as.is=T,col.names=c('Chr','Start','End','Name'))
  return(features)
}

getCytoColor <- function() {
  cytocolor <- c()

  cytocolor["gpos100"]  <- rgb(0,0,0,max=255)
  cytocolor["gpos"]     <- rgb(0,0,0,max=255)
  cytocolor["gpos75"]   <- rgb(130,130,130,max=255)
  cytocolor["gpos66"]   <- rgb(160,160,160,max=255)
  cytocolor["gpos50"]   <- rgb(200,200,200,max=255)
  cytocolor["gpos33"]   <- rgb(210,210,210,max=255)
  cytocolor["gpos25"]   <- rgb(200,200,200,max=255)
  cytocolor["gvar"]     <- rgb(220,220,220,max=255)
  cytocolor["gneg"]     <- rgb(255,255,255,max=255)
  cytocolor["acen"]     <- rgb(217,47,39,max=255)
  cytocolor["stalk"]    <- rgb(100,127,164,max=255)

  return(cytocolor)
}

formatBP <- function(x,pts=1) {
  unlist(lapply(x, function (bp) {
    if (bp > 1000000000) paste(round(bp/1000000000,pts),"Gb") else
    if (bp > 1000000) paste(round(bp/1000000,pts),"Mb") else
    if (bp > 1000) paste(round(bp/1000,pts),"Kb") else
    paste(bp,"bp")
  }))
}

