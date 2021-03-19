library(ggplot2)
library(ggrepel)
library(ggpubr)

outputD <- "New/LiverCancer/Catholic/1.miRNA/Figure_Correlation/"
dir.create(file.path(outputD), showWarnings = FALSE)

data.cat <- read.table(file="New/LiverCancer/Catholic/1.miRNA/Correlation_Mature_MajorIsomiR.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
data.tcga <- read.table(file="New/LiverCancer/TCGA/1.miRNA/Correlation_Mature_MajorIsomiR.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
data.tsinghua <- read.table(file="New/LiverCancer/Tsinghua/1.miRNA/Correlation_Mature_MajorIsomiR.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
arm <- "FP"

for (arm in c("FP","TP")){
  data.cat1 <- data.cat[data.cat$Arm == arm,]
  data.cat1$outlier <- "Others"
  data.tcga1 <- data.tcga[data.tcga$Arm == arm,]
  data.tcga1$outlier <- "Others"
  data.tsinghua1 <- data.tsinghua[data.tsinghua$Arm == arm,]
  data.tsinghua1$outlier <- "Others"
  if (arm == "FP"){
    data.cat1[data.cat1$Precursor %in% c("hsa-mir-192","hsa-mir-122") & data.cat1$Arm == "FP",]$outlier <- "Outliers"
    data.tcga1[data.tcga1$Precursor %in% c("hsa-mir-192","hsa-mir-122") & data.tcga1$Arm == "FP",]$outlier <- "Outliers"
    data.tsinghua1[data.tsinghua1$Precursor %in% c("hsa-mir-192","hsa-mir-122") & data.tsinghua1$Arm == "FP",]$outlier <- "Outliers"
  }
  #data1$outlier <- factor(data1$outlier, levels=c("Others","Outliers"))
  
  temp.data.cat1 <- data.cat1[data.cat1$Precursor %in% data.tcga1$Precursor,]
  temp.data.tcga1 <- data.tcga1[data.tcga1$Precursor %in% data.cat1$Precursor,]
  temp.data.tcga1 <- temp.data.tcga1[match(temp.data.cat1$Precursor,temp.data.tcga1$Precursor),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data1 <- cbind(temp.data.cat1[,colnames(temp.data.cat1) %in% c("Precursor","Arm","Median_Ratio_Major_1","SD_Ratio_Major_1")], temp.data.tcga1[,colnames(temp.data.tcga1) %in% c("Median_Ratio_Major_1","SD_Ratio_Major_1","outlier")])
  #head(temp.data1)
  colnames(temp.data1) <- c("Precursor","Arm","Cat_Frequency","Cat_Frequency_SD","TCGA_Frequency","TCGA_Frequency_SD","outlier")
  temp.data1 <- temp.data1[!is.na(temp.data1$Cat_Frequency) & !is.na(temp.data1$TCGA_Frequency),]
  p1 <- ggplot(temp.data1,aes(x=Cat_Frequency,y=TCGA_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=TCGA_Frequency-TCGA_Frequency_SD,ymax=TCGA_Frequency+TCGA_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","red")) + 
    ggtitle(label=bquote("n ="~.(nrow(temp.data1))~";"~italic(rho)~"="~.(round(cor(temp.data1$Cat_Frequency,temp.data1$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data1$Cat_Frequency,temp.data1$TCGA_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Catholic)") + 
    ylab("5' heterogeneity frequency (TCGA)") + 
    coord_fixed()
  ggsave(plot=p1,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Cat_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)

  temp.data.cat2 <- data.cat1[data.cat1$Precursor %in% data.tsinghua1$Precursor,]
  temp.data.tsinghua2 <- data.tsinghua1[data.tsinghua1$Precursor %in% data.cat1$Precursor,]
  temp.data.tsinghua2 <- temp.data.tsinghua2[match(temp.data.cat2$Precursor,temp.data.tsinghua2$Precursor),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data2 <- cbind(temp.data.cat2[,colnames(temp.data.cat2) %in% c("Precursor","Arm","Median_Ratio_Major_1","SD_Ratio_Major_1")], temp.data.tsinghua2[,colnames(temp.data.tsinghua2) %in% c("Median_Ratio_Major_1","SD_Ratio_Major_1","outlier")])
  #head(temp.data1)
  colnames(temp.data2) <- c("Precursor","Arm","Cat_Frequency","Cat_Frequency_SD","Tsinghwa_Frequency","Tsinghua_Frequency_SD","outlier")
  
  temp.data2 <- temp.data2[!is.na(temp.data2$Cat_Frequency) & !is.na(temp.data2$Tsinghwa_Frequency),]
  p2 <- ggplot(temp.data2,aes(x=Cat_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","red")) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data2))~";"~italic(rho)~"="~.(round(cor(temp.data2$Cat_Frequency,temp.data2$Tsinghwa_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data2$Cat_Frequency,temp.data2$Tsinghwa_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Catholic)") + 
    ylab("5' heterogeneity frequency (Tsinghua)") + 
    coord_fixed()
  ggsave(plot=p2,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Cat_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)

  temp.data.tcga3 <- data.tcga1[data.tcga1$Precursor %in% data.tsinghua1$Precursor,]
  temp.data.tsinghua3 <- data.tsinghua1[data.tsinghua1$Precursor %in% data.tcga1$Precursor,]
  temp.data.tsinghua3 <- temp.data.tsinghua3[match(temp.data.tcga3$Precursor,temp.data.tsinghua3$Precursor),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data3 <- cbind(temp.data.tcga3[,colnames(temp.data.tcga3) %in% c("Precursor","Arm","Median_Ratio_Major_1","SD_Ratio_Major_1")], temp.data.tsinghua3[,colnames(temp.data.tsinghua3) %in% c("Median_Ratio_Major_1","SD_Ratio_Major_1","outlier")])
  #head(temp.data1)
  colnames(temp.data3) <- c("Precursor","Arm","TCGA_Frequency","TCGA_Frequency_SD","Tsinghwa_Frequency","Tsinghua_Frequency_SD","outlier")
  
  temp.data3 <- temp.data3[!is.na(temp.data3$TCGA_Frequency) & !is.na(temp.data3$Tsinghwa_Frequency),]

  p3 <- ggplot(temp.data3,aes(x=TCGA_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=TCGA_Frequency-TCGA_Frequency_SD,xmax=TCGA_Frequency+TCGA_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","red")) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data3))~";"~italic(rho)~"="~.(round(cor(temp.data3$Tsinghwa_Frequency,temp.data3$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data3$Tsinghwa_Frequency,temp.data3$TCGA_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (TCGA)") + 
    ylab("5' heterogeneity frequency (Tsinghua)") + 
    coord_fixed()
  ggsave(plot=p3,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_TCGA_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)

  temp.data.cat4 <- data.cat1[data.cat1$Precursor %in% data.tcga1$Precursor,]
  temp.data.tcga4 <- data.tcga1[data.tcga1$Precursor %in% data.cat1$Precursor,]
  temp.data.tcga4 <- temp.data.tcga4[match(temp.data.cat4$Precursor,temp.data.tcga4$Precursor),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data4 <- cbind(temp.data.cat4[,colnames(temp.data.cat4) %in% c("Precursor","Arm","Median_Ratio_Major_2","SD_Ratio_Major_2")], temp.data.tcga4[,colnames(temp.data.tcga4) %in% c("Median_Ratio_Major_2","SD_Ratio_Major_2","outlier")])
  #head(temp.data1)
  colnames(temp.data4) <- c("Precursor","Arm","Cat_Frequency","Cat_Frequency_SD","TCGA_Frequency","TCGA_Frequency_SD","outlier")
  temp.data4 <- temp.data4[!is.na(temp.data4$Cat_Frequency) & !is.na(temp.data4$TCGA_Frequency),]
  p4 <- ggplot(temp.data4,aes(x=Cat_Frequency,y=TCGA_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=TCGA_Frequency-TCGA_Frequency_SD,ymax=TCGA_Frequency+TCGA_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","red")) + 
    ggtitle(label=bquote("n ="~.(nrow(temp.data4))~";"~italic(rho)~"="~.(round(cor(temp.data4$Cat_Frequency,temp.data4$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data4$Cat_Frequency,temp.data4$TCGA_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Catholic)") + 
    ylab("5' heterogeneity frequency (TCGA)") + 
    coord_fixed()
  ggsave(plot=p4,file=paste(outputD,"/CorrelationAll_",arm,"_SecondMajor_frequency_Cat_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)

  temp.data.cat5 <- data.cat1[data.cat1$Precursor %in% data.tsinghua1$Precursor,]
  temp.data.tsinghua5 <- data.tsinghua1[data.tsinghua1$Precursor %in% data.cat1$Precursor,]
  temp.data.tsinghua5 <- temp.data.tsinghua5[match(temp.data.cat5$Precursor,temp.data.tsinghua5$Precursor),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data5 <- cbind(temp.data.cat5[,colnames(temp.data.cat5) %in% c("Precursor","Arm","Median_Ratio_Major_2","SD_Ratio_Major_2")], temp.data.tsinghua5[,colnames(temp.data.tsinghua5) %in% c("Median_Ratio_Major_2","SD_Ratio_Major_2","outlier")])
  #head(temp.data1)
  colnames(temp.data5) <- c("Precursor","Arm","Cat_Frequency","Cat_Frequency_SD","Tsinghwa_Frequency","Tsinghua_Frequency_SD","outlier")
  
  temp.data5 <- temp.data5[!is.na(temp.data5$Cat_Frequency) & !is.na(temp.data5$Tsinghwa_Frequency),]
  p5 <- ggplot(temp.data5,aes(x=Cat_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","red")) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data5))~";"~italic(rho)~"="~.(round(cor(temp.data5$Cat_Frequency,temp.data5$Tsinghwa_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data5$Cat_Frequency,temp.data5$Tsinghwa_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Catholic)") + 
    ylab("5' heterogeneity frequency (Tsinghua)") + 
    coord_fixed()
  ggsave(plot=p5,file=paste(outputD,"/CorrelationAll_",arm,"_SecondMajor_frequency_Cat_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)

  temp.data.tcga6 <- data.tcga1[data.tcga1$Precursor %in% data.tsinghua1$Precursor,]
  temp.data.tsinghua6 <- data.tsinghua1[data.tsinghua1$Precursor %in% data.tcga1$Precursor,]
  temp.data.tsinghua6 <- temp.data.tsinghua6[match(temp.data.tcga6$Precursor,temp.data.tsinghua6$Precursor),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data6 <- cbind(temp.data.tcga6[,colnames(temp.data.tcga6) %in% c("Precursor","Arm","Median_Ratio_Major_2","SD_Ratio_Major_2")], temp.data.tsinghua6[,colnames(temp.data.tsinghua6) %in% c("Median_Ratio_Major_2","SD_Ratio_Major_2","outlier")])
  #head(temp.data1)
  colnames(temp.data6) <- c("Precursor","Arm","TCGA_Frequency","TCGA_Frequency_SD","Tsinghwa_Frequency","Tsinghua_Frequency_SD","outlier")
  
  temp.data6 <- temp.data6[!is.na(temp.data6$TCGA_Frequency) & !is.na(temp.data6$Tsinghwa_Frequency),]
  
  p6 <- ggplot(temp.data6,aes(x=TCGA_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=TCGA_Frequency-TCGA_Frequency_SD,xmax=TCGA_Frequency+TCGA_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","red")) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data6))~";"~italic(rho)~"="~.(round(cor(temp.data6$Tsinghwa_Frequency,temp.data6$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data6$Tsinghwa_Frequency,temp.data6$TCGA_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (TCGA)") + 
    ylab("5' heterogeneity frequency (Tsinghua)") + 
    coord_fixed()
  ggsave(plot=p6,file=paste(outputD,"/CorrelationAll_",arm,"_SecondMajor_frequency_TCGA_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)

  ################################################################################################################################
  #colnames(data.cat1)
  
  anova.sigs2 <- intersect(intersect(data.cat1[data.cat1$Anova_Major_1 <= 0.05 & ! is.na(data.cat1$Anova_Major_1),]$Precursor, data.tcga1[data.tcga1$Anova_Major_1 <= 0.05 & ! is.na(data.tcga1$Anova_Major_1),]$Precursor), data.tsinghua1[data.tsinghua1$Anova_Major_1 <= 0.05 & ! is.na(data.tsinghua1$Anova_Major_1),]$Precursor)
  anova.sigs <- intersect(data.cat1[data.cat1$Anova_Major_1 <= 0.05 & ! is.na(data.cat1$Anova_Major_1),]$Precursor, data.tsinghua1[data.tsinghua1$Anova_Major_1 <= 0.05 & ! is.na(data.tsinghua1$Anova_Major_1),]$Precursor)
  
  temp.data7 <- data.cat1[!is.na(data.cat1$Median_Ratio_Major_1),]
  temp.data7$Anova <- "ns"
  temp.data7[temp.data7$Precursor %in% anova.sigs,]$Anova <- "Sig"
  temp.data7[temp.data7$Precursor %in% anova.sigs2,]$Anova <- "Sig2"

  if (arm == "FP"){
    temp.data7[temp.data7$Precursor == "hsa-mir-21" & temp.data7$Arm == "FP",]$Anova <- "Sig3"
  }
  temp.data7$Anova <- factor(temp.data7$Anova, levels=c("ns","Sig","Sig2","Sig3"))
  p7 <- ggplot(temp.data7,aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_1-SD_Ratio_Cancer_Major_1,ymax=Median_Ratio_Cancer_Major_1+SD_Ratio_Cancer_Major_1),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_1-SD_Ratio_Normal_Major_1,xmax=Median_Ratio_Normal_Major_1+SD_Ratio_Normal_Major_1),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data7[temp.data7$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data7))~";"~italic(rho)~"="~.(round(cor(temp.data7$Median_Ratio_Normal_Major_1,temp.data7$Median_Ratio_Cancer_Major_1),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data7$Median_Ratio_Normal_Major_1,temp.data7$Median_Ratio_Cancer_Major_1)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p7,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p7.inset <- ggplot(temp.data7[temp.data7$Anova != "ns",],aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_1-SD_Ratio_Cancer_Major_1,ymax=Median_Ratio_Cancer_Major_1+SD_Ratio_Cancer_Major_1),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_1-SD_Ratio_Normal_Major_1,xmax=Median_Ratio_Normal_Major_1+SD_Ratio_Normal_Major_1),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data7, aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data7[temp.data7$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()

  ggsave(plot=p7.inset,file=paste(outputD,"/Correlation_",arm,"_FirstMajor_frequency_inset_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)

  temp.data8 <- data.tcga1[!is.na(data.tcga1$Median_Ratio_Major_1),]
  temp.data8$Anova <- "ns"
  temp.data8[temp.data8$Precursor %in% anova.sigs,]$Anova <- "Sig"
  temp.data8[temp.data8$Precursor %in% anova.sigs2,]$Anova <- "Sig2"
  
  if (arm == "FP"){
    temp.data8[temp.data8$Precursor == "hsa-mir-21" & temp.data8$Arm == "FP",]$Anova <- "Sig3"
  }
  temp.data8$Anova <- factor(temp.data8$Anova, levels=c("ns","Sig","Sig2","Sig3"))
  p8 <- ggplot(temp.data8,aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_1-SD_Ratio_Cancer_Major_1,ymax=Median_Ratio_Cancer_Major_1+SD_Ratio_Cancer_Major_1),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_1-SD_Ratio_Normal_Major_1,xmax=Median_Ratio_Normal_Major_1+SD_Ratio_Normal_Major_1),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data8[temp.data8$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data8))~";"~italic(rho)~"="~.(round(cor(temp.data8$Median_Ratio_Normal_Major_1,temp.data8$Median_Ratio_Cancer_Major_1),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data8$Median_Ratio_Normal_Major_1,temp.data8$Median_Ratio_Cancer_Major_1)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p8,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p8.inset <- ggplot(temp.data8[temp.data8$Anova != "ns",],aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_1-SD_Ratio_Cancer_Major_1,ymax=Median_Ratio_Cancer_Major_1+SD_Ratio_Cancer_Major_1),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_1-SD_Ratio_Normal_Major_1,xmax=Median_Ratio_Normal_Major_1+SD_Ratio_Normal_Major_1),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data8, aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data8[temp.data8$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","black","black","black"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p8.inset,file=paste(outputD,"/Correlation_",arm,"_FirstMajor_frequency_inset_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data9 <- data.tsinghua1[!is.na(data.tsinghua1$Median_Ratio_Major_1),]
  temp.data9$Anova <- "ns"
  temp.data9[temp.data9$Precursor %in% anova.sigs,]$Anova <- "Sig"
  temp.data9[temp.data9$Precursor %in% anova.sigs2,]$Anova <- "Sig2"
  
  if (arm == "FP"){
    temp.data9[temp.data9$Precursor == "hsa-mir-21" & temp.data9$Arm == "FP",]$Anova <- "Sig3"
  }
  temp.data9$Anova <- factor(temp.data9$Anova, levels=c("ns","Sig","Sig2","Sig3"))
  p9 <- ggplot(temp.data9,aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_1-SD_Ratio_Cancer_Major_1,ymax=Median_Ratio_Cancer_Major_1+SD_Ratio_Cancer_Major_1),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_1-SD_Ratio_Normal_Major_1,xmax=Median_Ratio_Normal_Major_1+SD_Ratio_Normal_Major_1),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data9[temp.data9$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data9))~";"~italic(rho)~"="~.(round(cor(temp.data9$Median_Ratio_Normal_Major_1,temp.data9$Median_Ratio_Cancer_Major_1),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data9$Median_Ratio_Normal_Major_1,temp.data9$Median_Ratio_Cancer_Major_1)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p9,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p9.inset <- ggplot(temp.data9[temp.data9$Anova != "ns",],aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_1-SD_Ratio_Cancer_Major_1,ymax=Median_Ratio_Cancer_Major_1+SD_Ratio_Cancer_Major_1),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_1-SD_Ratio_Normal_Major_1,xmax=Median_Ratio_Normal_Major_1+SD_Ratio_Normal_Major_1),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data9, aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data9[temp.data9$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_1,y=Median_Ratio_Cancer_Major_1,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","black","black","black"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p9.inset,file=paste(outputD,"/Correlation_",arm,"_FirstMajor_frequency_inset_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  ################################################################################################################################
  #colnames(data.cat1)
  
  anova.sigs2 <- intersect(intersect(data.cat1[data.cat1$Anova_Major_2 <= 0.05 & ! is.na(data.cat1$Anova_Major_2),]$Precursor, data.tcga1[data.tcga1$Anova_Major_2 <= 0.05 & ! is.na(data.tcga1$Anova_Major_2),]$Precursor), data.tsinghua1[data.tsinghua1$Anova_Major_2 <= 0.05 & ! is.na(data.tsinghua1$Anova_Major_2),]$Precursor)
  anova.sigs <- intersect(data.cat1[data.cat1$Anova_Major_2 <= 0.05 & ! is.na(data.cat1$Anova_Major_2),]$Precursor, data.tsinghua1[data.tsinghua1$Anova_Major_2 <= 0.05 & ! is.na(data.tsinghua1$Anova_Major_2),]$Precursor)
  
  temp.data7 <- data.cat1[!is.na(data.cat1$Median_Ratio_Major_2),]
  temp.data7$Anova <- "ns"
  temp.data7[temp.data7$Precursor %in% anova.sigs,]$Anova <- "Sig"
  temp.data7[temp.data7$Precursor %in% anova.sigs2,]$Anova <- "Sig2"
  
  if (arm == "FP"){
    temp.data7[temp.data7$Precursor == "hsa-mir-21" & temp.data7$Arm == "FP",]$Anova <- "Sig3"
  }
  temp.data7$Anova <- factor(temp.data7$Anova, levels=c("ns","Sig","Sig2","Sig3"))
  p7 <- ggplot(temp.data7,aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_2-SD_Ratio_Cancer_Major_2,ymax=Median_Ratio_Cancer_Major_2+SD_Ratio_Cancer_Major_2),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_2-SD_Ratio_Normal_Major_2,xmax=Median_Ratio_Normal_Major_2+SD_Ratio_Normal_Major_2),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data7[temp.data7$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data7))~";"~italic(rho)~"="~.(round(cor(temp.data7$Median_Ratio_Normal_Major_2,temp.data7$Median_Ratio_Cancer_Major_2),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data7$Median_Ratio_Normal_Major_2,temp.data7$Median_Ratio_Cancer_Major_2)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p7,file=paste(outputD,"/CorrelationAll_",arm,"_SecondMajor_frequency_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p7.inset <- ggplot(temp.data7[temp.data7$Anova != "ns",],aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_2-SD_Ratio_Cancer_Major_2,ymax=Median_Ratio_Cancer_Major_2+SD_Ratio_Cancer_Major_2),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_2-SD_Ratio_Normal_Major_2,xmax=Median_Ratio_Normal_Major_2+SD_Ratio_Normal_Major_2),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data7, aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data7[temp.data7$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p7.inset,file=paste(outputD,"/Correlation_",arm,"_SecondMajor_frequency_inset_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data8 <- data.tcga1[!is.na(data.tcga1$Median_Ratio_Major_2),]
  temp.data8$Anova <- "ns"
  temp.data8[temp.data8$Precursor %in% anova.sigs,]$Anova <- "Sig"
  temp.data8[temp.data8$Precursor %in% anova.sigs2,]$Anova <- "Sig2"
  
  if (arm == "FP"){
    temp.data8[temp.data8$Precursor == "hsa-mir-21" & temp.data8$Arm == "FP",]$Anova <- "Sig3"
  }
  temp.data8$Anova <- factor(temp.data8$Anova, levels=c("ns","Sig","Sig2","Sig3"))
  p8 <- ggplot(temp.data8,aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_2-SD_Ratio_Cancer_Major_2,ymax=Median_Ratio_Cancer_Major_2+SD_Ratio_Cancer_Major_2),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_2-SD_Ratio_Normal_Major_2,xmax=Median_Ratio_Normal_Major_2+SD_Ratio_Normal_Major_2),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data8[temp.data8$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data8))~";"~italic(rho)~"="~.(round(cor(temp.data8$Median_Ratio_Normal_Major_2,temp.data8$Median_Ratio_Cancer_Major_2),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data8$Median_Ratio_Normal_Major_2,temp.data8$Median_Ratio_Cancer_Major_2)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p8,file=paste(outputD,"/CorrelationAll_",arm,"_SecondMajor_frequency_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p8.inset <- ggplot(temp.data8[temp.data8$Anova != "ns",],aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_2-SD_Ratio_Cancer_Major_2,ymax=Median_Ratio_Cancer_Major_2+SD_Ratio_Cancer_Major_2),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_2-SD_Ratio_Normal_Major_2,xmax=Median_Ratio_Normal_Major_2+SD_Ratio_Normal_Major_2),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data8, aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data8[temp.data8$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","black","black","black"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p8.inset,file=paste(outputD,"/Correlation_",arm,"_SecondMajor_frequency_inset_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data9 <- data.tsinghua1[!is.na(data.tsinghua1$Median_Ratio_Major_2),]
  temp.data9$Anova <- "ns"
  temp.data9[temp.data9$Precursor %in% anova.sigs,]$Anova <- "Sig"
  temp.data9[temp.data9$Precursor %in% anova.sigs2,]$Anova <- "Sig2"
  
  if (arm == "FP"){
    temp.data9[temp.data9$Precursor == "hsa-mir-21" & temp.data9$Arm == "FP",]$Anova <- "Sig3"
  }
  temp.data9$Anova <- factor(temp.data9$Anova, levels=c("ns","Sig","Sig2","Sig3"))
  p9 <- ggplot(temp.data9,aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_2-SD_Ratio_Cancer_Major_2,ymax=Median_Ratio_Cancer_Major_2+SD_Ratio_Cancer_Major_2),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_2-SD_Ratio_Normal_Major_2,xmax=Median_Ratio_Normal_Major_2+SD_Ratio_Normal_Major_2),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data9[temp.data9$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data9))~";"~italic(rho)~"="~.(round(cor(temp.data9$Median_Ratio_Normal_Major_2,temp.data9$Median_Ratio_Cancer_Major_2),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data9$Median_Ratio_Normal_Major_2,temp.data9$Median_Ratio_Cancer_Major_2)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p9,file=paste(outputD,"/CorrelationAll_",arm,"_SecondMajor_frequency_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p9.inset <- ggplot(temp.data9[temp.data9$Anova != "ns",],aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer_Major_2-SD_Ratio_Cancer_Major_2,ymax=Median_Ratio_Cancer_Major_2+SD_Ratio_Cancer_Major_2),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal_Major_2-SD_Ratio_Normal_Major_2,xmax=Median_Ratio_Normal_Major_2+SD_Ratio_Normal_Major_2),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data9, aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data9[temp.data9$Anova != "ns",], aes(x=Median_Ratio_Normal_Major_2,y=Median_Ratio_Cancer_Major_2,label=Precursor),size=3.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
    scale_color_manual(values=c("black","black","black","black"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p9.inset,file=paste(outputD,"/Correlation_",arm,"_SecondMajor_frequency_inset_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
}

