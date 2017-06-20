#!/usr/bin/env Rscript

#----------------------------------------------#
# Script to convert diploid structure output   #
#                                              #
# to haploid BAPS input			       #
#                                              #
#----------------------------------------------#

#f <- file("stdin")
outfile<-"./Baps_"
f<-"./batch_1.structure.tsv"
structure<-read.table(f,skip=2,header=F,row.names=NULL)
structure<-structure[seq(1,dim(structure)[1],2),]
pops<-structure[,1:2]
structure<-structure[,-c(1:2)]
structure[structure[]==0]<--9
structure[structure[]==1]<-0
structure[structure[]==2]<-1
structure[structure[]==3]<-2
structure[structure[]==4]<-3

write.table(cbind(structure,1:dim(structure)[1]),paste(outfile,"single.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(cbind(structure,as.numeric(pops[,2])),paste(outfile,"pops.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
write.table(cbind(pops,as.numeric(pops[,2])),paste(outfile,"populations.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
