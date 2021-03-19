library(ggplot2)
library(ggrepel)
library(ggpubr)

outputD <- "New/LiverCancer/Catholic/1.miRNA/Figure_Correlation2/"
dir.create(file.path(outputD), showWarnings = FALSE)

data.cat <- read.table(file="New/LiverCancer/Catholic/1.miRNA/Correlation_Mature_AllIsomiRs.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
data.tcga <- read.table(file="New/LiverCancer/TCGA/1.miRNA/Correlation_Mature_AllIsomiRs.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
data.tsinghua <- read.table(file="New/LiverCancer/Tsinghua/1.miRNA/Correlation_Mature_AllIsomiRs.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
arm <- "TP"

for (arm in c("FP","TP")){
  data.cat1 <- data.cat[data.cat$Arm == arm,]
  data.cat1$outlier <- "Others"
  data.tcga1 <- data.tcga[data.tcga$Arm == arm,]
  data.tcga1$outlier <- "Others"
  data.tsinghua1 <- data.tsinghua[data.tsinghua$Arm == arm,]
  data.tsinghua1$outlier <- "Others"
  if (arm == "FP"){
    data.cat1[data.cat1$isomiR %in% c("hsa-miR-192-5p_GACCUAU_1","hsa-miR-122-5p_GAGUGUG_1") & data.cat1$Arm == "FP",]$outlier <- "Outliers"
    data.tcga1[data.tcga1$isomiR %in% c("hsa-miR-192-5p_GACCUAU_1","hsa-miR-122-5p_GAGUGUG_1") & data.tcga1$Arm == "FP",]$outlier <- "Outliers"
    data.tsinghua1[data.tsinghua1$isomiR %in% c("hsa-miR-192-5p_GACCUAU_1","hsa-miR-122-5p_GAGUGUG_1") & data.tsinghua1$Arm == "FP",]$outlier <- "Outliers"
  }
  #data1$outlier <- factor(data1$outlier, levels=c("Others","Outliers"))
  
  temp.data.cat1 <- data.cat1[data.cat1$isomiR %in% data.tcga1$isomiR,]
  temp.data.tcga1 <- data.tcga1[data.tcga1$isomiR %in% data.cat1$isomiR,]
  temp.data.tcga1 <- temp.data.tcga1[match(temp.data.cat1$isomiR,temp.data.tcga1$isomiR),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data1 <- cbind(temp.data.cat1[,colnames(temp.data.cat1) %in% c("isomiR","Precursor","Arm","Median_Ratio","SD_Ratio")], temp.data.tcga1[,colnames(temp.data.tcga1) %in% c("Median_Ratio","SD_Ratio","outlier")])
  #head(temp.data1)
  colnames(temp.data1) <- c("isomiR","Precursor","Arm","Cat_Frequency","Cat_Frequency_SD","TCGA_Frequency","TCGA_Frequency_SD","outlier")
  temp.data1 <- temp.data1[!is.na(temp.data1$Cat_Frequency) & !is.na(temp.data1$TCGA_Frequency),]
  p1 <- ggplot(temp.data1,aes(x=Cat_Frequency,y=TCGA_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=TCGA_Frequency-TCGA_Frequency_SD,ymax=TCGA_Frequency+TCGA_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","red")) + 
    ggtitle(label=bquote("n ="~.(nrow(temp.data1))~";"~italic(rho)~"="~.(round(cor(temp.data1$Cat_Frequency,temp.data1$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data1$Cat_Frequency,temp.data1$TCGA_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Catholic)") + 
    ylab("5' heterogeneity frequency (TCGA)") + 
    coord_fixed()
  ggsave(plot=p1,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Cat_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data.cat2 <- data.cat1[data.cat1$isomiR %in% data.tsinghua1$isomiR,]
  temp.data.tsinghua2 <- data.tsinghua1[data.tsinghua1$isomiR %in% data.cat1$isomiR,]
  temp.data.tsinghua2 <- temp.data.tsinghua2[match(temp.data.cat2$isomiR,temp.data.tsinghua2$isomiR),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data2 <- cbind(temp.data.cat2[,colnames(temp.data.cat2) %in% c("isomiR","Precursor","Arm","Median_Ratio","SD_Ratio")], temp.data.tsinghua2[,colnames(temp.data.tsinghua2) %in% c("Median_Ratio","SD_Ratio","outlier")])
  #head(temp.data1)
  colnames(temp.data2) <- c("isomiR","Precursor","Arm","Cat_Frequency","Cat_Frequency_SD","Tsinghwa_Frequency","Tsinghua_Frequency_SD","outlier")
  
  temp.data2 <- temp.data2[!is.na(temp.data2$Cat_Frequency) & !is.na(temp.data2$Tsinghwa_Frequency),]
  p2 <- ggplot(temp.data2,aes(x=Cat_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","red")) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data2))~";"~italic(rho)~"="~.(round(cor(temp.data2$Cat_Frequency,temp.data2$Tsinghwa_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data2$Cat_Frequency,temp.data2$Tsinghwa_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Catholic)") + 
    ylab("5' heterogeneity frequency (Tsinghua)") + 
    coord_fixed()
  ggsave(plot=p2,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Cat_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data.tcga3 <- data.tcga1[data.tcga1$isomiR %in% data.tsinghua1$isomiR,]
  temp.data.tsinghua3 <- data.tsinghua1[data.tsinghua1$isomiR %in% data.tcga1$isomiR,]
  temp.data.tsinghua3 <- temp.data.tsinghua3[match(temp.data.tcga3$isomiR,temp.data.tsinghua3$isomiR),]
  #rownames(temp.data.cat1) == rownames(temp.data.tcga1)
  temp.data3 <- cbind(temp.data.tcga3[,colnames(temp.data.tcga3) %in% c("isomiR","Precursor","Arm","Median_Ratio","SD_Ratio")], temp.data.tsinghua3[,colnames(temp.data.tsinghua3) %in% c("Median_Ratio","SD_Ratio","outlier")])
  #head(temp.data1)
  colnames(temp.data3) <- c("isomiR","Precursor","Arm","TCGA_Frequency","TCGA_Frequency_SD","Tsinghwa_Frequency","Tsinghua_Frequency_SD","outlier")
  
  temp.data3 <- temp.data3[!is.na(temp.data3$TCGA_Frequency) & !is.na(temp.data3$Tsinghwa_Frequency),]
  
  p3 <- ggplot(temp.data3,aes(x=TCGA_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
    geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=TCGA_Frequency-TCGA_Frequency_SD,xmax=TCGA_Frequency+TCGA_Frequency_SD),height=0.2,color="gray70") +
    geom_point(size=1.5) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","red")) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data3))~";"~italic(rho)~"="~.(round(cor(temp.data3$Tsinghwa_Frequency,temp.data3$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data3$Tsinghwa_Frequency,temp.data3$TCGA_Frequency)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (TCGA)") + 
    ylab("5' heterogeneity frequency (Tsinghua)") + 
    coord_fixed()
  ggsave(plot=p3,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_TCGA_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  ################################################################################################################################
  
  anova.sigs2 <- intersect(intersect(data.cat1[data.cat1$Anova <= 0.05 & ! is.na(data.cat1$Anova),]$isomiR, data.tcga1[data.tcga1$Anova <= 0.05 & ! is.na(data.tcga1$Anova),]$isomiR), data.tsinghua1[data.tsinghua1$Anova <= 0.05 & ! is.na(data.tsinghua1$Anova),]$isomiR)
  anova.sigs <- intersect(data.cat1[data.cat1$Anova <= 0.05 & ! is.na(data.cat1$Anova),]$isomiR, data.tsinghua1[data.tsinghua1$Anova <= 0.05 & ! is.na(data.tsinghua1$Anova),]$isomiR)
  
  temp.data7 <- data.cat1[!is.na(data.cat1$Median_Ratio),]
  temp.data7$Anova <- "ns"
  temp.data7[temp.data7$isomiR %in% anova.sigs,]$Anova <- "Sig"
  temp.data7[temp.data7$isomiR %in% anova.sigs2,]$Anova <- "Sig2"
  
  #if (arm == "FP"){
  #  temp.data7[temp.data7$Precursor == "hsa-mir-21" & temp.data7$Arm == "FP",]$Anova <- "Sig3"
  #}
  temp.data7$Anova <- factor(temp.data7$Anova, levels=c("ns","Sig","Sig2"))
  p7 <- ggplot(temp.data7,aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer-SD_Ratio_Cancer,ymax=Median_Ratio_Cancer+SD_Ratio_Cancer),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal-SD_Ratio_Normal,xmax=Median_Ratio_Normal+SD_Ratio_Normal),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data7[temp.data7$Anova != "ns",], aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,label=isomiR),size=1) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","Green","red"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data7))~";"~italic(rho)~"="~.(round(cor(temp.data7$Median_Ratio_Normal,temp.data7$Median_Ratio_Cancer),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data7$Median_Ratio_Normal,temp.data7$Median_Ratio_Cancer)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p7,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p7.inset <- ggplot(temp.data7[temp.data7$Anova != "ns",],aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer-SD_Ratio_Cancer,ymax=Median_Ratio_Cancer+SD_Ratio_Cancer),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal-SD_Ratio_Normal,xmax=Median_Ratio_Normal+SD_Ratio_Normal),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data7, aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data7[temp.data7$Anova != "ns",], aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,label=isomiR),size=1) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p7.inset,file=paste(outputD,"/Correlation_",arm,"_FirstMajor_frequency_inset_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data8 <- data.tcga1[!is.na(data.tcga1$Median_Ratio),]
  temp.data8$Anova <- "ns"
  temp.data8[temp.data8$isomiR %in% anova.sigs,]$Anova <- "Sig"
  temp.data8[temp.data8$isomiR %in% anova.sigs2,]$Anova <- "Sig2"
  
  temp.data8$Anova <- factor(temp.data8$Anova, levels=c("ns","Sig","Sig2"))
  p8 <- ggplot(temp.data8,aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer-SD_Ratio_Cancer,ymax=Median_Ratio_Cancer+SD_Ratio_Cancer),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal-SD_Ratio_Normal,xmax=Median_Ratio_Normal+SD_Ratio_Normal),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data8[temp.data8$Anova != "ns",], aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,label=isomiR),size=1) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","Green","red","blue"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data8))~";"~italic(rho)~"="~.(round(cor(temp.data8$Median_Ratio_Normal,temp.data8$Median_Ratio_Cancer),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data8$Median_Ratio_Normal,temp.data8$Median_Ratio_Cancer)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p8,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p8.inset <- ggplot(temp.data8[temp.data8$Anova != "ns",],aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer-SD_Ratio_Cancer,ymax=Median_Ratio_Cancer+SD_Ratio_Cancer),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal-SD_Ratio_Normal,xmax=Median_Ratio_Normal+SD_Ratio_Normal),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data8, aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data8[temp.data8$Anova != "ns",], aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,label=isomiR),size=1) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","black","black","black"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p8.inset,file=paste(outputD,"/Correlation_",arm,"_FirstMajor_frequency_inset_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  temp.data9 <- data.tsinghua1[!is.na(data.tsinghua1$Median_Ratio),]
  temp.data9$Anova <- "ns"
  temp.data9[temp.data9$isomiR %in% anova.sigs,]$Anova <- "Sig"
  temp.data9[temp.data9$isomiR %in% anova.sigs2,]$Anova <- "Sig2"
  
  temp.data9$Anova <- factor(temp.data9$Anova, levels=c("ns","Sig","Sig2"))
  p9 <- ggplot(temp.data9,aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer-SD_Ratio_Cancer,ymax=Median_Ratio_Cancer+SD_Ratio_Cancer),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal-SD_Ratio_Normal,xmax=Median_Ratio_Normal+SD_Ratio_Normal),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data9[temp.data9$Anova != "ns",], aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,label=isomiR),size=1) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","Green","red"),guide=F) +
    ggtitle(label=bquote("n ="~.(nrow(temp.data9))~";"~italic(rho)~"="~.(round(cor(temp.data9$Median_Ratio_Normal,temp.data9$Median_Ratio_Cancer),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data9$Median_Ratio_Normal,temp.data9$Median_Ratio_Cancer)$p.value,digits=3)))) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  ggsave(plot=p9,file=paste(outputD,"/CorrelationAll_",arm,"_FirstMajor_frequency_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
  
  p9.inset <- ggplot(temp.data9[temp.data9$Anova != "ns",],aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,color=Anova)) +
    geom_errorbar(aes(ymin=Median_Ratio_Cancer-SD_Ratio_Cancer,ymax=Median_Ratio_Cancer+SD_Ratio_Cancer),width=0.2,color="gray70") +
    geom_errorbarh(aes(xmin=Median_Ratio_Normal-SD_Ratio_Normal,xmax=Median_Ratio_Normal+SD_Ratio_Normal),height=0.2,color="gray70") +
    geom_point(size=1.5) +
    geom_smooth(inherit.aes = F, data = temp.data9, aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer), method=lm,color="lightslateblue",fill="lightslateblue") +
    geom_abline(slope = 1, linetype ="dashed") + 
    geom_text_repel(inherit.aes = F, data = temp.data9[temp.data9$Anova != "ns",], aes(x=Median_Ratio_Normal,y=Median_Ratio_Cancer,label=isomiR),size=1) + 
    scale_x_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6.5,2.5)) + 
    scale_color_manual(values=c("black","black","black","black"),guide=F) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
    xlab("5' heterogeneity frequency (Normal)") + 
    ylab("5' heterogeneity frequency (HCC)") + 
    coord_fixed()
  
  ggsave(plot=p9.inset,file=paste(outputD,"/Correlation_",arm,"_FirstMajor_frequency_inset_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
  
}
