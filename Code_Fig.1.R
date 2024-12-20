#####Fragmentomics data from ONT (NGDC accession: OMIX007720) will be released on the date of publication.
CR1_frag=read.table("CR1.stats")
WB1_frag=read.table("WB1.stats")
CR2_frag=read.table("CR2.stats")
CR3_frag=read.table("CR3.stats")
CR4_frag=read.table("CR4.stats")
HR2_frag=read.table("HR2.stats")
HR1_frag=read.table("HR1.stats")
CR5_frag=read.table("CR5.stats")
CR6_frag=read.table("CR6.stats")
library(ggplot2)
library(ggsci)
library(scales)
library(gghalves)
CR1_frag=CR1_frag$V15
WB1_frag=WB1_frag$V15
CR3_frag=CR3_frag$V15
CR2_frag=CR2_frag$V15
CR4_frag=CR4_frag$V15
HR1_frag=HR1_frag$V15
HR2_frag=HR2_frag$V15
CR5_frag=CR5_frag$V15
CR6_frag=CR6_frag$V15

NAME=c("WB1","CR1","CR2","CR3","CR4","CR5","CR6","HR1","HR2")
res=matrix(NA,9,3)
for (t in 1:9){
  LS=get(stringr::str_c(NAME[t],"_frag"))
  res[t,1]=length(LS[LS<=1000])/length(LS)
  res[t,2]=length(LS[LS>1000 & LS<8000])/length(LS)
  res[t,3]=length(LS[LS>=8000])/length(LS)
}
rownames(res)=NAME
group=c("WBC",rep("Cancer",6),"Normal","Normal")
rownames(res)=c("WB1","CR1","CR2","CR3","CR4","CR5","CR6","HR1","HR2")
colnames(res)=c("<1kb","1~8kb",">8kb")
res=as.data.frame(res)
res$group=group
LS=reshape2::melt(res)
##############Fig.1C
ggplot(LS,aes(fill=group,y=value,x=variable))+geom_half_boxplot(outlier.alpha = 0)+ggthemes::theme_few()+theme(panel.grid.major = element_blank(),axis.text = element_text(face= "italic",color="black",size=10),axis.title = element_text(size=12,face= "bold"),title = element_text(size=15,face= "bold"),axis.text.x = element_text(angle=90),panel.border = element_blank(),axis.line.x=element_blank(),axis.line.y = element_blank(),panel.background = element_blank(),axis.ticks.y  = element_line(linetype = 8,linewidth = 1,lineend = "square"))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL)) +labs(x="",y="Fragment proportion")+geom_half_point(aes(color=group),side='r',size=0.5,alpha=0.2)+ggsci::scale_fill_lancet()+ggsci::scale_color_lancet()

##############Fig.1D
library(ggpubr)
P1<-ggplot(NULL)+geom_histogram(aes(x=CR1_frag[CR1_frag<1000]),binwidth = 10,color=pal_lancet("lanonc")(9)[2],fill="white")+geom_vline(xintercept = c(170,370,570,770,970),linetype=2,size=0.8,color=pal_lancet("lanonc")(9)[2])+theme_classic()+theme(axis.text = element_text(color="black"))+labs(x="CR1")
P2<-ggplot(NULL)+geom_histogram(aes(x=CR4_frag[CR4_frag<1000]),binwidth = 10,color=pal_lancet("lanonc")(9)[5],fill="white")+geom_vline(xintercept = c(170,370,570,770,970),linetype=2,size=0.8,color=pal_lancet("lanonc")(9)[5])+theme_classic()+theme(axis.text = element_text(color="black"))+labs(x="CR4")
P3<-ggplot(NULL)+geom_histogram(aes(x=CR6_frag[CR6_frag<1000]),binwidth = 10,color=pal_lancet("lanonc")(9)[7],fill="white")+geom_vline(xintercept = c(170,370,570,770,970),linetype=2,size=0.8,color=pal_lancet("lanonc")(9)[7])+theme_classic()+theme(axis.text = element_text(color="black"))+labs(x="CR6")
ggarrange(P1,P2,P3,ncol = 3,nrow = 1)

#########Fig.1E
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
CR1_frag=read.table("CR1.stats")
WB1_frag=read.table("WB1.stats")
CR2_frag=read.table("CR2.stats")
CR3_frag=read.table("CR3.stats")
CR4_frag=read.table("CR4.stats")
HR2_frag=read.table("HR2.stats")
HR1_frag=read.table("HR1.stats")
CR5_frag=read.table("CR5.stats")
CR6_frag=read.table("CR6.stats")
NAME=c("WB1","CR1","CR2","CR3","CR4","CR5","CR6","HR1","HR2")
peaklist=list()
for (t in 1:9){
  LS=get(stringr::str_c(NAME[t],"_frag"))
  colnames(LS)[1:3]=c("chr","start","end")
  LS=LS[,-7:-14]
  LS$NB=1
  peaklist[[t]]=GRanges(LS)
  
}
names(peaklist)=NAME
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaklist, getTagMatrix, windows=promoter)
names(tagMatrixList)=c("WB1","CR1","CR2","CR3","CR4","CR5","CR6","HR1","HR2")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))+theme(panel.grid.major = element_blank(),axis.text = element_text(face= "italic",color="black",size=10),axis.title = element_text(size=12,face= "bold"),title = element_text(size=10,face= "bold"),legend.text = element_text(size=10),axis.text.x = element_text(angle=90),axis.line.x=element_blank(),axis.line.y = element_blank())+ guides(colour=guide_legend(title=NULL))

#########Fig.1F
peaklist=list()
for (t in 1:9){
  LS=get(stringr::str_c(NAME[t],"_frag"))
  colnames(LS)[1:3]=c("chr","start","end")
  LS=LS[,-7:-14]
  LS$NB=1
  LS=LS[LS$V15<=1000,]
  peaklist[[t]]=GRanges(LS)
  
}
names(peaklist)=NAME
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaklist, getTagMatrix, windows=promoter)
names(tagMatrixList)=c("WB1","CR1","CR2","CR3","CR4","CR5","CR6","HR1","HR2")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))+theme(panel.grid.major = element_blank(),axis.text = element_text(face= "italic",color="black",size=10),axis.title = element_text(size=12,face= "bold"),title = element_text(size=10,face= "bold"),legend.text = element_text(size=10),axis.text.x = element_text(angle=90),axis.line.x=element_blank(),axis.line.y = element_blank())+ guides(colour=guide_legend(title=NULL))

############Fig.1G
MOTIF_1G=read.delim("clipboard") ##Copy Fig. 1G data from the source Data.XLSX file, and import it into R using the clipboard
library(ggplot2)
ggplot(MOTIF_1G,aes(x=factor(name,levels = c("CR2","CR6","CR3","CR5","CR4","CR1","WB1","HR1","HR2")),y=value,fill=Var1))+geom_bar(stat = "identity",position = "stack")+theme(panel.grid.major = element_blank(),axis.text = element_text(face= "italic",color="black",size=10),axis.title = element_text(size=12,face= "bold"),title = element_text(size=15,face= "bold"),axis.text.x = element_text(angle=90),panel.border = element_blank(),axis.line.x=element_blank(),axis.line.y = element_blank(),panel.background = element_blank(),legend.position ="top",axis.ticks.y  = element_line(linetype = 8,linewidth = 1,lineend = "square"))+labs(x="",y="Contribution of profiles")+ggsci::scale_fill_lancet()+scale_y_continuous(expand = c(0,0))+guides(fill = guide_legend(title = 'Profiles'))

############Fig.1H
coef_matrix=as.matrix(read.delim("clipboard",row.names = NULL)) ##Copy Fig. 1H data from the source Data.XLSX file, and import it into R using the clipboard
LS=reshape2::melt(coef_matrix)
LS$group=NA
LS$group=stringr::str_sub(LS$Var2,1,1)
LSdata=LS[LS$Var1==1,]
LSdata$F1=stringr::str_sub(LSdata$Var2,1,1)
LSdata$F2=stringr::str_sub(LSdata$Var2,2,2)
LSdata$F3=stringr::str_sub(LSdata$Var2,3,3)
LSdata$F4=stringr::str_sub(LSdata$Var2,4,4)
MOFdata=matrix(NA,4,4)
rownames(MOFdata)=c("A","C","G","T")
colnames(MOFdata)=1:4
for (m in 1:4){
  for (n in 1:4) {
    LS1=LSdata[LSdata[,stringr::str_c("F",as.character(n))]==(c("A","C","G","T")[m]),]
    MOFdata[m,n]=sum(LS1$value)
  }
}
MOFdata=MOFdata*10
MOFdata=round(MOFdata)
library("ggseqlogo")
P1=ggplot(LS[LS$Var1==1,],aes(x=Var2,y=value,fill=group))+geom_bar(stat="identity")+theme(panel.grid.major = element_blank(),axis.text = element_text(face= "italic",color="black",size=10),axis.title = element_text(size=12,face= "bold"),title = element_text(size=15,face= "bold"),axis.text.x = element_text(angle=90),panel.border = element_blank(),axis.line.x=element_blank(),axis.line.y = element_blank(),panel.background = element_blank(),axis.ticks.y  = element_line(linetype = 8,linewidth = 1,lineend = "square"))+scale_x_discrete(breaks=c("AAAA","ACAA","AGAA","ATAA","CAAA","CCAA","CGAA","CTAA","GAAA","GCAA","GGAA","GTAA","TAAA","TCAA","TGAA","TTAA"))+labs(x="4-mer motif",y="Frequency by NMF",title = "End-motif profile I")+ggsci::scale_fill_lancet()+scale_y_continuous(expand = c(0,0))+coord_cartesian(ylim = c(0,0.2))+guides(fill=F)
P3=ggplot()+geom_logo(MOFdata)+theme_logo()+theme(panel.grid.major = element_blank(),axis.text = element_text(color="black",size=10),axis.title = element_text(size=10,face= "bold"),title = element_text(size=10,face= "bold"),axis.text.x = element_text(),panel.border = element_blank(),axis.line.x=element_blank(),axis.line.y = element_blank(),panel.background = element_blank(),axis.ticks.y  = element_line(linetype = 8,linewidth = 1,lineend = "square"))+scale_y_continuous(expand = c(0,0))+coord_cartesian(ylim = c(0,0.2))
LS2=aggregate(value~group,LSdata,sum)
LS2$group=stringr::str_c(LS2$group,"-end")
P2=ggplot(LS2,aes(y=group,x=value,fill=group))+geom_bar(stat="identity")+theme(panel.grid.major = element_blank(),axis.text = element_text(color="black",size=10),axis.text.y = element_text(face= "bold",color="black",size=10),axis.title = element_text(size=10),title = element_text(size=15),axis.text.x = element_text(),panel.border = element_blank(),axis.line.x=element_blank(),axis.line.y = element_blank(),panel.background = element_blank(),axis.ticks.y  = element_line(linetype = 8,linewidth = 1,lineend = "square"))+labs(y="",x="Frequency by NMF")+ggsci::scale_fill_lancet()+guides(fill=F)
library(patchwork)
design <- "
  111112
  111113
"
P1 + P2 + P3 + plot_layout(design =design)





