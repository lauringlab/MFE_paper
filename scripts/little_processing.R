require(plyr)

## This script take in a difference file and adjusts the positions to represent the distance from the translation start site.
wsn33.bed<-read.csv("../data/diff_wsn33_HA.csv",stringsAsFactors = F,comment.char = "#")
data<-read.csv("../data/MFE_HA.csv")
coding.adjust<-function(x){
  chr<-unique(x$chr)
  start<-wsn33.bed$off.5[match(x$chr,wsn33.bed$chr)]

  mutate(x,coding.pos=pos-start)
}

data<-ddply(data,~chr,coding.adjust)

write.csv(subset(data,pos<1699),"../data/MFE_HA.codingpos.csv")
