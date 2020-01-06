setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/cldn6")

memo="lncRNA"

dgefile<-"lncRNA.dge.tcga.pancancer.smd.meta.pvalue.csv"
osfile<-"lncRNA.os.hr.tcga.pancancer.meta.pvalue.csv"

dge<-read.csv(dgefile)
os<-read.csv(osfile)

out<-merge(dge,os,by="X")

pick<-subset(out,pval<10^-6 & pval.fixed<0.05 & pval.random<0.05 & beta*(TE.fixed-1)>0)

dim(pick)

write.csv(out,file=paste(memo,".tcga.pancancer.chol.overall.rnaseq.dmg.smd",".os.hr.csv",sep=""),quote=F)
write.table(out,file=paste(memo,".tcga.pancancer.chol.overall.rnaseq.dmg.smd",".os.hr.txt",sep=""),quote=F,col.names=NA,row.names=T,sep="\t")

write.csv(pick,file=paste(memo,".tcga.pancancer.chol.pick.rnaseq.dmg.smd","os.hr.csv",sep=""),quote=F)
write.table(pick,file=paste(memo,".tcga.pancancer.chol.pick.rnaseq.dmg.smd","os.hr.txt",sep=""),quote=F,col.names=NA,row.names=T,sep="\t")

dir.create("pick")

# be careful to run it or else it will takes quite long time

for(i in 1:nrow(pick)){
x<-paste("cp *-",pick[i,1],"* ./pick",sep="")
system(x)
}

# install.packages("CMplot")

library("CMplot")
ensg<-read.table("~/hpc/db/hg19/ENSG.hg19.bed")
ensgrlt<-ensg[match(unlist(lapply(strsplit(as.character(out$X),"[.]"),function(x) x[1])),ensg[,5]),]
output<-data.frame(ensgrlt,out)

cminput<-na.omit(data.frame(SNP=output$V5,Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval,stringsAsFactors=F))
CMplot(cminput,plot.type="b",memo=paste(memo,".fix",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.table(cminput,file=paste(memo,".pval.manhattan.qqplot.meta.dge.txt",sep=""),sep="\t",quote=F,row.name=T,col.names=NA)

cminput<-data.frame(SNP=output$V5,Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval.fix)
CMplot(cminput,plot.type="b",memo=paste(memo,".fix",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.table(cminput,file=paste(memo,".pval.fix.manhattan.qqplot.meta.OS.HR.txt",sep=""),sep="\t",quote=F,row.name=T,col.names=NA)

cminput<-data.frame(SNP=output$V5,Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval.random)
CMplot(cminput,plot.type="b",memo=paste(memo,".random",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.table(cminput,file=paste(memo,".pval.random.manhattan.qqplot.meta.OS.HR.txt",sep=""),sep="\t",quote=F,row.name=T,col.names=NA)






