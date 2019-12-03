#Useful R subroutines
library(tools)

# first argument is a start and end coordinate pair defining a span
# second argument is dataframe with Tstart and Tend columns
# returns the number of basepairs in the span covered by the dataframe 
bedIntersect <- function(x,df) {
  df <- df[df$Tend >= x[1] & df$Tstart <= x[2],]

  overlaps <- apply(df[,c("Tstart","Tend")],1,function(d) {
      start <- max(x[1],d[1])
      end <- min(x[2],d[2])
      overlap <- max(0,end-start+1)
      return(overlap)
    })
  return(sum(overlaps))
}


# first argument is a start and end coordinate pair defining a span
# second argument is dataframe with Tstart column
# returns the number of elements in dataframe with Tstart in the span 
bedCount <- function(x,df) {
  return(nrow(df[df$Tstart >= x[1] & df$Tstart <= x[2],]))
}


parseArgs <- function(scriptname, ARGS, OPTS=NULL) {
  if (length(ARGS) %% 3 != 0)
    stop("Error: ARGS must have length divisable by 3")
  if (!is.null(OPTS) && length(OPTS) %% 4 != 0)
    stop("Error: OPTS must have length divisable by 4")
  ARGS <- matrix(ARGS, nrow=length(ARGS)/3, ncol=3, byrow=T, dimnames=list(c(),c("name","type","description")))
  OPTS <- if (!is.null(OPTS)) matrix(OPTS, nrow=length(OPTS)/4, ncol=4, byrow=T, dimnames=list(c(),c("name","type","default","description")))
  
  ARGV <- commandArgs(trailingOnly=T)
  
  catargs <- paste(ARGS[,"name"],ARGS[,"description"],sep=" - ",collapse="\n")
  catopts <- if (!is.null(OPTS)) paste(OPTS[,"name"],"(",OPTS[,"default"],") - ", OPTS[,"description"], sep=" ", collapse="\n")
  
  usage <- paste("\nUsage: Rscript ",scriptname," ",paste(ARGS[,"name"],collapse=" ")," [OPTS] ","\n\nArguments (defaults):\n\n", catargs,"\n",catopts,"\n\n",sep="")
  
  if (length(ARGV) < nrow(ARGS))
    stop("Not enough arguments.\n",usage)
  
  types <- c("character","numeric","integer","logical")
  
  for (i in 1:nrow(ARGS)) {
    if (! ARGS[i,"type"] %in% types)
      stop("Unrecognized argument type ",ARGS[i,"type"]," for argument ",ARGS[i,"name"],"\n",usage)
    cmd <- paste(ARGS[i,"name"],"<<-as.",ARGS[i,"type"],"(\"",ARGV[i],"\")",sep="")
    eval(parse(text=cmd))
  }
  
  if (!is.null(OPTS)) for (k in 1:nrow(OPTS)) {
    if (!OPTS[k,"type"] %in% types)
      stop("Unrecognized argument type ",OPTS[k,"type"]," for optional argument ",OPTS[k,"name"],"\n",usage)
    cmd <- paste(OPTS[k,"name"],"<<-as.",OPTS[k,"type"],"(\"",OPTS[k,"default"],"\")",sep="")
    eval(parse(text=cmd))
  }
  
  if (length(ARGV) > i) for (j in (i+1):length(ARGV)) {
    opt <- unlist(strsplit(ARGV[j],"="))
    if (length(opt) != 2)
      stop("opt=value expression expected instead of ",ARGV[j],"\n",usage)
    if (is.null(OPTS) || ! opt[1] %in% OPTS[,"name"])
      stop("optional argument ",opt[1]," not recognized\n",usage)
    idx <- match(opt[1],OPTS[,"name"])
    cmd <- paste(opt[1],"<<-as.",OPTS[idx,"type"],"(\"",opt[2],"\")",sep="")
    eval(parse(text=cmd))
  }
}


read.delim.save <- function(file, do.save=T, force.read=F, ext=".RData", verbose=F, ...) {
  binfile <- paste(file_path_sans_ext(file),ext,sep="")
  if ( !force.read && file.access(binfile,mode=4)==0 ) {
    catv(verbose, "\nLoading binary version: ",binfile,"...",sep="")
    var <- load(binfile)
    catv(verbose, "Done\n")
    return(eval(parse(text=var)))
  } else {
    catv(verbose, "\nReading ",file," into memory...",sep="")
    dat <- read.delim(file,...)
    catv(verbose, "Done\n")
  }
  if ( do.save ) {
    catv(verbose, "\nSaving binary version to ",binfile,"...",sep="")
    save(dat,file=binfile)
    catv(verbose, "Done\n")
  }
  return(dat)
}

catv <- function(verbose, ...) {
  if (verbose == T) {
    cat(...);
  }
}

readHeader <- function (filename) {
  con  <- file(filename, open = "r")
  header <- unlist(strsplit(readLines(con, n = 1),"\t"))
  close(con)
  return(header)
}

proColumnClasses <- function (header, columnsToRead) {
  colClasses <- rep("NULL",length(header))
  colClasses[match(columnsToRead,header)] <- NA
  return(colClasses)
}

antiColumnClasses <- function (header, columnsToIgnore) {
  colClasses <- rep("NULL",length(header))
  colClasses[-match(columnsToRead,header)] <- NA
  return(colClasses)
}