library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)
library(data.table)

dataset <- "TCGA_rev1" # catholic TCGA_rev1 Tsinghua

miRNA <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA//Total/Frequency/AlleleFrequencyDistAll_Mature_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep='\t')
isomir <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1Minus1_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep="\t")
fold <- read.table(paste("Project/LiverCancer/",dataset,"/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_All.txt",sep=""),header=T,row.names = 1,check.names=F,quote="",sep="\t")

miRNA.fp <- miRNA[rownames(miRNA) %like% "FP",]
miRNA.tp <- miRNA[rownames(miRNA) %like% "TP",]
isomir.fp <- isomir[rownames(isomir) %like% "FP",]
isomir.tp <- isomir[rownames(isomir) %like% "TP",]
fold.fp <- fold[rownames(fold) %like% "FP",]
fold.tp <- fold[rownames(fold) %like% "TP",]

#setdiff(rownames(isomir.fp.p),rownames(isomir.fp.m))

cv <- function(exp){
  return(sd(exp) / mean(exp))
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
                  labels = trans_format("log10", math_format(10^.x)), limits=c(10^-2,10^2)) +
    theme_bw() + 
    theme(panel.grid=element_blank(),axis.title = element_text(size=12,color="black"), axis.text = element_text(size=12,color="black"),axis.ticks = element_line(color="black"),axis.line=element_line(color="black")) +
    xlab("") + ylab("Standard deviation") + ggtitle(paste("n = ",length(set),sep=""))# + ylim(0,5)
  temp.scatterplot <- ggplot(mature.cv,aes(x=Mean,y=SD)) + geom_point() #+ ylim(0,100) + xlim(0,1000)
  return(temp.dotplot)
}
#mature <- miRNA.fp
#iso <- isomir.fp
#ratio <- fold.fp
a <- create.dotplot(log10(miRNA.fp+0.1), log10(isomir.fp+0.1), fold.fp)
b <- create.dotplot(log10(miRNA.tp+0.1), log10(isomir.tp+0.1), fold.tp)
#b <- create.dotplot(log10(miRNA.tp+0.1), log10(isomir.normal+0.1), fold.normal)
#c <- create.dotplot(log10(miRNA.cancer+0.1), log10(isomir.cancer+0.1), fold.cancer)
c <- ggarrange(a,b,labels=c("FP","TP"),nrow=1,ncol=2)
ggsave(c,file=paste("New/LiverCancer/Catholic/5.IsomiR/CV_",dataset,"_alllog.pdf",sep=""),dpi=300,width=10,height=5)
