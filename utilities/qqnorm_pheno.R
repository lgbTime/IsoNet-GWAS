#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
data<-read.table(args[1], row.names=1)
name=args[2]
pdfname = paste(args[2], ".pdf", sep="")
print(pdfname)
pdf(pdfname)
a<-qqnorm(data[,1])
b<-as.data.frame(a$x)
row.names(b)<-row.names(data)
write.table(b,name,sep="\t",quote=F,col.names=F)
dev.off()
