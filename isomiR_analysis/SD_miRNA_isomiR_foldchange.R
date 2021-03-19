library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)

dataset <- "catholic" # catholic TCGA_rev1 Tsinghua

miRNA.fp <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD//Total/Frequency/AlleleFrequencyDistAll_Mature_RPM_FP.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep='\t')
isomir.fp <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD_test/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_RPM_FP.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep="\t")
#isomir.fp.p <- read.table("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1_RPM_FP.txt",header=T,row.names=1,check.names=F,quote="",sep="\t")[1:15]
fold.fp <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_FP.txt",sep=""),header=T,row.names = 1,check.names=F,quote="",sep="\t")
#fold.fp.p <- read.table("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1_Frequency_FP.txt",header=T,row.names = 1,check.names=F,quote="",sep="\t")[1:15]
miRNA.tp <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD//Total/Frequency/AlleleFrequencyDistAll_Mature_RPM_TP.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep='\t')
isomir.tp <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD_test/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_RPM_TP.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep="\t")
#isomir.fp.p <- read.table("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1_RPM_FP.txt",header=T,row.names=1,check.names=F,quote="",sep="\t")[1:15]
fold.tp <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_TP.txt",sep=""),header=T,row.names = 1,check.names=F,quote="",sep="\t")

nrow(miRNA.fp)
nrow(isomir.fp)
#nrow(isomir.fp)
nrow(fold.fp)
#nrow(fold.fp.p)

#setdiff(rownames(isomir.fp.p),rownames(isomir.fp.m))

cv <- function(exp){
  return(sd(exp))# / mean(exp))
}

calculate.cv <- function(data){
  temp.frame <- data.frame()
  for (name in rownames(data)){
    temp.exp <- as.numeric(as.character(data[rownames(data) == name,]))
    temp.frame <- rbind(temp.frame, data.frame(matrix(c(name,mean(temp.exp),sd(temp.exp)/length(temp.exp), cv(temp.exp)),nrow=1,ncol=4),check.names=F))
  }
  colnames(temp.frame) <- c("Name","Mean","SD","CV")
  temp.frame$CV <- as.numeric(as.character(temp.frame$CV))
  temp.frame$Mean <- as.numeric(as.character(temp.frame$Mean))
  temp.frame$SD <- as.numeric(as.character(temp.frame$SD))
  return(temp.frame) 
}

create.dotplot <- function(mature, iso, ratio){
  set <- intersect(rownames(iso),rownames(ratio))
  iso.filtered <- iso[rownames(iso) %in% set,]
  ratio.filtered <- ratio[rownames(ratio) %in% set,]
  
  total.frame <- data.frame()
  
  mature.cv <- calculate.cv(mature)
  iso.cv <- calculate.cv(iso.filtered)
  ratio.cv <- calculate.cv(ratio.filtered)
  
  mature.cv$Data <- "Mature"
  iso.cv$Data <- "IsomiR"
  ratio.cv$Data <- "Ratio"
  print(iso.cv)
  total.frame <- rbind(mature.cv, iso.cv)
  total.frame <- rbind(total.frame, ratio.cv)
  
  total.frame$Data <- factor(total.frame$Data, c("Mature","IsomiR","Ratio"))
  #"Mature",
  #"#8f3a39",
  temp.dotplot <- ggplot(total.frame,aes(x=Data,y=CV+1,fill=Data)) + 
    #geom_dotplot(binaxis="y",stackdir="center",dotsize=0.01,binwidth = 1) +
    #geom_violin(trim=F) +
    stat_boxplot(geom='errorbar',lwd=.3,position = position_dodge(width=.3),width=.3) + 
    geom_boxplot(position=position_dodge(width=.6),width=.6, fill="white") +
    scale_fill_manual(values=c("black","#61bc3f","#5d7db0"),guide=F) +
    #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #              labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits=c(10^-1,10^5)) +
    theme_bw() + 
    theme(panel.grid=element_blank(),axis.title = element_text(size=12,color="black"), axis.text = element_text(size=12,color="black"),axis.ticks = element_line(color="black"),axis.line=element_line(color="black")) +
    xlab("") + ylab("Standard deviation") + ggtitle(paste("n = ",length(set),sep=""))# + ylim(0,5)
  temp.scatterplot <- ggplot(mature.cv,aes(x=Mean,y=SD)) + geom_point() #+ ylim(0,100) + xlim(0,1000)
  return(temp.dotplot)
}
#mature <- miRNA.fp
#iso <- isomir.fp
#ratio <- fold.fp
a <- create.dotplot(miRNA.fp, isomir.fp, 2^fold.fp)
b <- create.dotplot(miRNA.tp, isomir.tp, 2^fold.tp)

c <- ggarrange(a,b,labels=c("FP","TP"),nrow=1,ncol=2)
ggsave(c,file=paste("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD//Total/Frequency/Figure_Correlation/CV_",dataset,"_notlog.pdf",sep=""),dpi=300,width=10,height=5)
