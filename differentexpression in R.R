setwd("C:/Users/Administrator/desktop/result")
input_data=read.table("C:/Users/Administrator/Desktop/result/diff.csv",sep=",",header = TRUE,row.names = 1)
input_data<-input_data[,-3]
input_data<-input_data[,1:4]
input_data

input_data<-as.matrix(input_data)
condition=factor(c(rep('ctl',2),rep('l',2)))
coldata<-data.frame(row.names = colnames(input_data),condition)
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = input_data,colData = coldata,design = ~condition)

dds<-DESeq(dds)
res<-results(dds,alpha = 0.05)
summary(res)
res=res[order(res$padj),]
resdata<-merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)), by='row.names',sort=FALSE)
names(resdata)[1]<-'Gene'
write.table(resdata,file='differentresults_E12.txt',sep='\t',quote=F,row.names = F)
plotMA(res)

maplot<-function(res,thresh=0.05,labelsig=TRUE,...){
  with(res,plot(baseMean,log2FoldChange,pch=20,cex=.5,log="x",...))
  with(subset(res,padj<thresh),points(baseMean,log2FoldChange,col="red",pch=20))
}
png("diffexpr-maplot.png",1500,1000,pointsize=20)
maplot(resdata,main="MAPlot")
dev.off()


#install.packages('ggrepel')
library(ggplot2)
library(ggrepel)


resdata$significant<-as.factor(resdata$padj<0.05 & abs(resdata$log2FoldChange)>1)
ggplot(data=resdata,aes(x=log2FoldChange,y=-log10(padj),color=significant))+
  geom_point()+
  ylim (0,8)+
  scale_color_manual(values =c('black','red'))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = 'Volcano plot',x='log2(fold change)',y='-log10(padj)')+
  theme(plot.title=element_text(hjust=0.5))+
  geom_text_repel(data=subset(resdata,-log10(padj)>6),aes(label=Gene),col='black',alpha=0.8)



resdata=read.table('differentresults_L1.txt',header=T)
head(resdata)
up<-subset(resdata, resdata$log2FoldChange>1 & resdata$padj<0.05)
write.table(up,'L1_up.txt')
down<-subset(resdata, resdata$log2FoldChange<-1 & resdata$padj<0.05)
write.table(down,'L1_down.txt')


