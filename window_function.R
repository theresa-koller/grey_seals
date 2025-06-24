# Genomic Window functions

# This script is highly inspired by Olkkonen & Löytynoja (2023).

# packages and wd

library(tidyverse)
library(data.table)
library(IRanges)
library(scales)

# note that the "work" directory needs to be specified

# mask and index

m <- read_table('work/maskdir/joint_mask.bed',col_names=c("ctg","start","end"), show_col_types=F)
m$start <- m$start+1
m$end <- m$end+1
f <- read.table("work/reference/ribbon.fa.fai")
colnames(f)[1] <- "ctg"
m_order <-  merge(f, m, by = "ctg", sort = FALSE)[,-c(2,3,4,5)]

# function to calculate windows

windows <- function(dat){
  colnames(dat) <- c("ctg", "pos", "pi")
  pos <- dat[,c(1,2)]
  colnames(pos) <- c("ctg","pos")
  ww <- 250000                                                     # window size
  dc <- c()
  d <- dat
  for(ctg in unique(pos$ctg)){
    #print(ctg)
    d1 <- d[as.character(d$ctg)==ctg,]                             # select one contig at a time
    d2 <- cbind(d1,floor(d1$pos/ww)+1)                             # add a new column indicating the window (each 250kbp in size) "true size"
    colnames(d2)[dim(d2)[2]] <- "window"                           # label it "window"
    # keep 2 df: with and without nan sites
    d2.nan <- d2                                                   
    d2.clean <- na.omit(d2)                                        
    d2.nan <- d2.nan[d2.nan$window %in% d2.clean$window,]          
    if(dim(d2.clean)[1]>2){                                        # if more than two positions
      # pi sum for each window
      pis <- d2.clean$pi
      d3 <- aggregate(pis,by=list(d2.clean$window),FUN=sum)        # sum per window
      colnames(d3) <- c("window", "pi")
      # effective length of window - joint mask
      m1 <- m_order[m_order$ctg==ctg,]                             # select joint masked regions for this contig
      m2 <- IRanges(start=m1$start,end=m1$end)                     # convert into IRanges
      w1 <- seq(1,f$V2[f$ctg==ctg],ww)                             # make 250kbp windows for this contig
      w2 <- IRanges(start=w1,width=ww)                             # convert into IRanges
      s <- rep(0,dim(d3)[1])
      t <- rep(0,dim(d3)[1])
      for(i in end(w2)/ww){                                        # iterate over each window ("i" being the index)
        s[i] <- sum(width(restrict(m2,start(w2)[i],end(w2)[i])))   # sum of joint masked regions
        t[i] <- i                                                  # window index
      }
      s <- cbind(t,s)
      colnames(s) <- c("window","nsites")
      s <- as.data.frame(s)
      s <- s[s$window %in% d2.clean$window,]                       
      # subtract nan sites from effective site of window
      nan <- d2.nan %>% group_by(window) %>% summarise(n_na = sum(is.na(pi)))
      nan <- as.data.frame(nan)
      s$nsites <- s$nsites-nan$n_na
      # put summed pis and effective length of window together
      d4 <- merge(d3,s,by.x=1,by.y=1)                              # merge the pi/dxy counts and joint mask counts
      dc <- rbind.data.frame(dc, cbind(ctg, d4))                   # add them to the final table, written to the disk in the end
    }
  }
  dc_mean <- dc
  dc_mean$pi <- dc_mean$pi/dc_mean$nsites
  return(dc_mean)
}


# function to calculate bins for line plot

bins.plot <- function(dat) {
  dat.bins <- hist(dat[dat$nsites>100000,pi],breaks=c(-0.001,seq(0.0001,0.0052,0.0001),0.030),plot=F)$counts
  dat.bins <- dat.bins/sum(dat.bins)*100    # % windows
  return(dat.bins)
}
