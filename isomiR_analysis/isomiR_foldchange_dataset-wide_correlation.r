library(ggplot2)
library(ggrepel)
library(ggpubr)

#types <- c("Normal","Cancer","Total")
types <- c("Total")
#arm <- "FP"
#type <- "Total"
for (type in types){
  inputD <- paste("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/",type,"/Frequency/",sep="")
  data <- read.table(file=paste(inputD,"AlleleFrequencyDist_OffsetPlusMiNum1_Frequency_IntersectSet_FDR.txt",sep=""), header=T,row.names=NULL,check.names=F,sep='\t',quote="")
  for (arm in levels(data$Arm)){
    data1 <- data[data$Arm == arm,]
    data1 <- data1[data1$Precursor != "hsa-mir-215",]
    data1$outlier <- "Others"
    if (arm == "FP"){
      data1[data1$Precursor %in% c("hsa-mir-192","hsa-mir-122") & data1$Arm == "FP",]$outlier <- "Outliers"
    }
    data1$outlier <- factor(data1$outlier, levels=c("Others","Outliers"))
    temp.data1 <- data1[!is.na(data1$Cat_Frequency) & !is.na(data1$TCGA_Frequency),]
    p1 <- ggplot(temp.data1,aes(x=Cat_Frequency,y=TCGA_Frequency, color=outlier)) +
      geom_errorbar(aes(ymin=TCGA_Frequency-TCGA_Frequency_SD,ymax=TCGA_Frequency+TCGA_Frequency_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      scale_x_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red")) + 
      ggtitle(label=bquote("n ="~.(nrow(temp.data1))~";"~italic(rho)~"="~.(round(cor(temp.data1$Cat_Frequency,temp.data1$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data1$Cat_Frequency,temp.data1$TCGA_Frequency)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Catholic)") + 
      ylab("5' heterogeneity frequency (TCGA)") + 
      coord_fixed()
    ggsave(plot=p1,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Cat_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p1,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Cat_TCGA.png",sep=""),dpi=300,width=5.5,height=4)
    temp.data2 <- data1[!is.na(data1$Cat_Frequency) & !is.na(data1$Tsinghwa_Frequency),]
    p2 <- ggplot(temp.data2,aes(x=Cat_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
      geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=Cat_Frequency-Cat_Frequency_SD,xmax=Cat_Frequency+Cat_Frequency_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      scale_x_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red")) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data2))~";"~italic(rho)~"="~.(round(cor(temp.data2$Cat_Frequency,temp.data2$Tsinghwa_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data2$Cat_Frequency,temp.data2$Tsinghwa_Frequency)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Catholic)") + 
      ylab("5' heterogeneity frequency (Tsinghua)") + 
      coord_fixed()
    ggsave(plot=p2,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Cat_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p2,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Cat_Tsinghwa.png",sep=""),dpi=300,width=5.5,height=4)
    temp.data3 <- data1[!is.na(data1$TCGA_Frequency) & !is.na(data1$Tsinghwa_Frequency),]
    p3 <- ggplot(temp.data3,aes(x=TCGA_Frequency,y=Tsinghwa_Frequency, color=outlier)) +
      geom_errorbar(aes(ymin=Tsinghwa_Frequency-Tsinghua_Frequency_SD,ymax=Tsinghwa_Frequency+Tsinghua_Frequency_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=TCGA_Frequency-TCGA_Frequency_SD,xmax=TCGA_Frequency+TCGA_Frequency_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      scale_x_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red")) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data3))~";"~italic(rho)~"="~.(round(cor(temp.data3$Tsinghwa_Frequency,temp.data3$TCGA_Frequency),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data3$Tsinghwa_Frequency,temp.data3$TCGA_Frequency)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (TCGA)") + 
      ylab("5' heterogeneity frequency (Tsinghua)") + 
      coord_fixed()
    ggsave(plot=p3,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_TCGA_Tsinghwa.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p3,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_TCGA_Tsinghwa.png",sep=""),dpi=300,width=5.5,height=4)
    
    temp.data4 <- data1[!is.na(data1$Cat_Frequency),]
    temp.data4$Anova <- "ns"
    temp.data4[temp.data4$OneWayAnova_Catholic <= 0.05 & temp.data4$OneWayAnova_TCGA <= 0.05 & temp.data4$OneWayAnova_Tsinghua <= 0.05 & !is.na(temp.data4$OneWayAnova_Catholic) & !is.na(temp.data4$OneWayAnova_TCGA) & !is.na(temp.data4$OneWayAnova_Tsinghua) & !is.na(temp.data4$Cat_Frequency) & temp.data4$Catholic_Correlation == TRUE & !is.na(temp.data4$TCGA_Frequency) & temp.data4$TCGA_Correlation == TRUE & !is.na(temp.data4$Tsinghwa_Frequency) & temp.data4$Tsinghwa_Correlation == TRUE,]$Anova <- "Sig"
    if (arm == "FP"){
      temp.data4[temp.data4$Precursor == "hsa-mir-21" & temp.data4$Arm == "FP",]$Anova <- "Sig2"
    }
    temp.data4$Anova <- factor(temp.data4$Anova, levels=c("ns","Sig","Sig2"))
    p4 <- ggplot(temp.data4,aes(x=Cat_Frequency_N,y=Cat_Frequency_C,color=Anova)) +
      geom_errorbar(aes(ymin=Cat_Frequency_C-Cat_Frequency_C_SD,ymax=Cat_Frequency_C+Cat_Frequency_C_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=Cat_Frequency_N-Cat_Frequency_N_SD,xmax=Cat_Frequency_N+Cat_Frequency_N_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) +
      geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1, linetype ="dashed") + 
      geom_text_repel(inherit.aes = F, data = temp.data4[temp.data4$Anova != "ns",], aes(x=Cat_Frequency_N,y=Cat_Frequency_C,label=Precursor),size=3.5) + 
      scale_x_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data4))~";"~italic(rho)~"="~.(round(cor(temp.data4$Cat_Frequency_N,temp.data4$Cat_Frequency_C),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data4$Cat_Frequency_N,temp.data4$Cat_Frequency_C)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Normal)") + 
      ylab("5' heterogeneity frequency (HCC)") + 
      coord_fixed()
    ggsave(plot=p4,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p4,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Cat.png",sep=""),dpi=300,width=5.5,height=4)
    p4.inset <- ggplot(temp.data4[temp.data4$Anova != "ns",],aes(x=Cat_Frequency_N,y=Cat_Frequency_C,color=Anova)) +
      geom_errorbar(aes(ymin=Cat_Frequency_C-Cat_Frequency_C_SD,ymax=Cat_Frequency_C+Cat_Frequency_C_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=Cat_Frequency_N-Cat_Frequency_N_SD,xmax=Cat_Frequency_N+Cat_Frequency_N_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) +
      geom_smooth(inherit.aes = F, data=temp.data4,aes(x=Cat_Frequency_N,y=Cat_Frequency_C),method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1, linetype ="dashed") + 
      geom_text_repel(inherit.aes = F, data = temp.data4[temp.data4$Anova != "ns",], aes(x=Cat_Frequency_N,y=Cat_Frequency_C,label=Precursor),size=3.5) + 
      scale_x_continuous(breaks=c(-5,-4,-3,-2,-1,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data4))~";"~italic(rho)~"="~.(round(cor(temp.data4$Cat_Frequency_N,temp.data4$Cat_Frequency_C),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data4$Cat_Frequency_N,temp.data4$Cat_Frequency_C)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Normal)") + 
      ylab("5' heterogeneity frequency (HCC)") + 
      coord_fixed()
    ggsave(plot=p4.inset,file=paste(inputD,"/Figure_Correlation/Correlation_",arm,"_offset_frequency_inset_Cat.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p4.inset,file=paste(inputD,"/Figure_Correlation/Correlation_",arm,"_offset_frequency_inset_Cat.png",sep=""),dpi=300,width=5.5,height=4)
    #temp.data5 <- data1[!is.na(data1$TCGA_Frequency) & data1$TCGA_Correlation == TRUE,]
    temp.data5 <- data1[!is.na(data1$TCGA_Frequency),]
    temp.data5$Anova <- "ns"
    temp.data5[temp.data5$OneWayAnova_Catholic <= 0.05 & temp.data5$OneWayAnova_TCGA <= 0.05 & temp.data5$OneWayAnova_Tsinghua <= 0.05 & !is.na(temp.data5$OneWayAnova_Catholic) & !is.na(temp.data5$OneWayAnova_TCGA) & !is.na(temp.data5$OneWayAnova_Tsinghua) & !is.na(temp.data5$Cat_Frequency) & temp.data5$Catholic_Correlation == TRUE & !is.na(temp.data5$TCGA_Frequency) & temp.data5$TCGA_Correlation == TRUE & !is.na(temp.data5$Tsinghwa_Frequency) & temp.data5$Tsinghwa_Correlation == TRUE,]$Anova <- "Sig"
    if (arm == "FP"){
      temp.data5[temp.data5$Precursor == "hsa-mir-21" & temp.data5$Arm == "FP",]$Anova <- "Sig2"
    }
    temp.data5$Anova <- factor(temp.data5$Anova, levels=c("ns","Sig","Sig2"))
    p5 <- ggplot(temp.data5,aes(x=TCGA_Frequency_N,y=TCGA_Frequency_C,color=Anova)) +
      geom_errorbar(aes(ymin=TCGA_Frequency_C-TCGA_Frequency_C_SD,ymax=TCGA_Frequency_C+TCGA_Frequency_C_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=TCGA_Frequency_N-TCGA_Frequency_N_SD,xmax=TCGA_Frequency_N+TCGA_Frequency_N_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1, linetype ="dashed") + 
      geom_text_repel(inherit.aes = F, data = temp.data5[temp.data5$Anova != "ns",], aes(x=TCGA_Frequency_N,y=TCGA_Frequency_C,label=Precursor),size=3.5) + 
      scale_x_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data5))~";"~italic(rho)~"="~.(round(cor(temp.data5$TCGA_Frequency_N,temp.data5$TCGA_Frequency_C),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data5$TCGA_Frequency_N,temp.data5$TCGA_Frequency_C)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Normal)") + 
      ylab("5' heterogeneity frequency (HCC)") + 
      coord_fixed()
    p5.inset <- ggplot(temp.data5[temp.data5$Anova != "ns",],aes(x=TCGA_Frequency_N,y=TCGA_Frequency_C,color=Anova)) +
      geom_errorbar(aes(ymin=TCGA_Frequency_C-TCGA_Frequency_C_SD,ymax=TCGA_Frequency_C+TCGA_Frequency_C_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=TCGA_Frequency_N-TCGA_Frequency_N_SD,xmax=TCGA_Frequency_N+TCGA_Frequency_N_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      geom_smooth(inherit.aes = F, data = temp.data5, aes(x=TCGA_Frequency_N,y=TCGA_Frequency_C), method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1, linetype ="dashed") + 
      geom_text_repel(inherit.aes = F, data = temp.data5[temp.data5$Anova != "ns",], aes(x=TCGA_Frequency_N,y=TCGA_Frequency_C,label=Precursor),size=3.5) + 
      scale_x_continuous(breaks=c(-5,-4,-3,-2,-1,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data5))~";"~italic(rho)~"="~.(round(cor(temp.data5$TCGA_Frequency_N,temp.data5$TCGA_Frequency_C),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data5$TCGA_Frequency_N,temp.data5$TCGA_Frequency_C)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Normal)") + 
      ylab("5' heterogeneity frequency (HCC)") + 
      coord_fixed()
    ggsave(plot=p5,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p5,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_TCGA.png",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p5.inset,file=paste(inputD,"/Figure_Correlation/Correlation_",arm,"_offset_frequency_inset_TCGA.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p5.inset,file=paste(inputD,"/Figure_Correlation/Correlation_",arm,"_offset_frequency_inset_TCGA.png",sep=""),dpi=300,width=5.5,height=4)
    temp.data6 <- data1[!is.na(data1$Tsinghwa_Frequency),]
    temp.data6$Anova <- "ns"
    
    temp.data6[temp.data6$OneWayAnova_Catholic <= 0.05 & temp.data6$OneWayAnova_TCGA <= 0.05 & temp.data6$OneWayAnova_Tsinghua <= 0.05 & !is.na(temp.data6$OneWayAnova_Catholic) & !is.na(temp.data6$OneWayAnova_TCGA) & !is.na(temp.data6$OneWayAnova_Tsinghua) & !is.na(temp.data6$Cat_Frequency) & temp.data6$Catholic_Correlation == TRUE & !is.na(temp.data6$TCGA_Frequency) & temp.data6$TCGA_Correlation == TRUE & !is.na(temp.data6$Tsinghwa_Frequency) & temp.data6$Tsinghwa_Correlation == TRUE,]$Anova <- "Sig"
    if (arm == "FP"){
      temp.data6[temp.data6$Precursor == "hsa-mir-21" & temp.data6$Arm == "FP",]$Anova <- "Sig2"
    }
    temp.data6$Anova <- factor(temp.data6$Anova, levels=c("ns","Sig","Sig2"))
    p6 <- ggplot(temp.data6,aes(x=Tsinghua_Frequency_N,y=Tsinghua_Frequency_C,color=Anova)) +
      geom_errorbar(aes(ymin=Tsinghua_Frequency_C-Tsinghua_Frequency_C_SD,ymax=Tsinghua_Frequency_C+Tsinghua_Frequency_C_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=Tsinghua_Frequency_N-Tsinghua_Frequency_N_SD,xmax=Tsinghua_Frequency_N+Tsinghua_Frequency_N_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1, linetype ="dashed") + 
      geom_text_repel(inherit.aes = F, data = temp.data6[temp.data6$Anova != "ns",], aes(x=Tsinghua_Frequency_N,y=Tsinghua_Frequency_C,label=Precursor),size=3.5) + 
      scale_x_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-4,-2,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data6))~";"~italic(rho)~"="~.(round(cor(temp.data6$Tsinghua_Frequency_N,temp.data6$Tsinghua_Frequency_C),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data6$Tsinghua_Frequency_N,temp.data6$Tsinghua_Frequency_C)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Normal)") + 
      ylab("5' heterogeneity frequency (HCC)") + 
      coord_fixed()
    p6.inset <- ggplot(temp.data6[temp.data6$Anova != "ns",],aes(x=Tsinghua_Frequency_N,y=Tsinghua_Frequency_C,color=Anova)) +
      geom_errorbar(aes(ymin=Tsinghua_Frequency_C-Tsinghua_Frequency_C_SD,ymax=Tsinghua_Frequency_C+Tsinghua_Frequency_C_SD),width=0.5,color="gray70") +
      geom_errorbarh(aes(xmin=Tsinghua_Frequency_N-Tsinghua_Frequency_N_SD,xmax=Tsinghua_Frequency_N+Tsinghua_Frequency_N_SD),height=0.5,color="gray70") +
      geom_point(size=1.5) + 
      geom_smooth(inherit.aes = F, data = temp.data6, aes(x=Tsinghua_Frequency_N,y=Tsinghua_Frequency_C), method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1, linetype ="dashed") + 
      geom_text_repel(inherit.aes = F, data = temp.data6[temp.data6$Anova != "ns",], aes(x=Tsinghua_Frequency_N,y=Tsinghua_Frequency_C,label=Precursor),size=3.5) + 
      scale_x_continuous(breaks=c(-5,-4,-3,-2,-1,0),limits=c(-5.5,.5)) + 
      scale_y_continuous(breaks=c(-5,-4,-3,-2,-1,0),limits=c(-5.5,.5)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data6))~";"~italic(rho)~"="~.(round(cor(temp.data6$Tsinghua_Frequency_N,temp.data6$Tsinghua_Frequency_C),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data6$Tsinghua_Frequency_N,temp.data6$Tsinghua_Frequency_C)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (Normal)") + 
      ylab("5' heterogeneity frequency (HCC)") + 
      coord_fixed()
    ggsave(plot=p6,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p6,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_offset_frequency_Tsinghua.png",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p6.inset,file=paste(inputD,"/Figure_Correlation/Correlation_",arm,"_offset_frequency_inset_Tsinghua.pdf",sep=""),dpi=300,width=5.5,height=4)
    ggsave(plot=p6.inset,file=paste(inputD,"/Figure_Correlation/Correlation_",arm,"_offset_frequency_inset_Tsinghua.png",sep=""),dpi=300,width=5.5,height=4)
    
    temp.data7 <- data1[!is.na(data1$Cat_Frequency_Down) & !is.na(data1$Cat_Frequency_Up),]
    head(temp.data7)
    
    p7 <- ggplot(temp.data7,aes(x=log10(Cat_Frequency_Down),y=log10(Cat_Frequency_Up))) +
      geom_point(size=1.5) + 
      geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1,linetype ="dashed") + 
      geom_abline(slope = 1, intercept = c(1,-1), linetype="dashed", color="red")+
      scale_x_continuous(breaks=c(-6,-4,-2,0),limits=c(-6,0)) + 
      scale_y_continuous(breaks=c(-6,-4,-2,0),limits=c(-6,0)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data7))~";"~italic(rho)~"="~.(round(cor(temp.data7$Cat_Frequency_Down,temp.data7$Cat_Frequency_Up),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data7$Cat_Frequency_Down,temp.data7$Cat_Frequency_Up)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (offset: +1)") + 
      ylab("5' heterogeneity frequency (offset: -1)") + 
      coord_fixed()
    
    temp.data8 <- data1[!is.na(data1$TCGA_Frequency_Down) & !is.na(data1$TCGA_Frequency_Up),]
    
    p8 <- ggplot(temp.data8,aes(x=log10(TCGA_Frequency_Down),y=log10(TCGA_Frequency_Up))) +
      geom_point(size=1.5) + 
      geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1,linetype ="dashed") + 
      geom_abline(slope = 1, intercept = c(1,-1), linetype="dashed", color="red")+
      scale_x_continuous(breaks=c(-6,-4,-2,0),limits=c(-6,0)) + 
      scale_y_continuous(breaks=c(-6,-4,-2,0),limits=c(-6,0)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data8))~";"~italic(rho)~"="~.(round(cor(temp.data8$TCGA_Frequency_Down,temp.data8$TCGA_Frequency_Up),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data8$TCGA_Frequency_Down,temp.data8$TCGA_Frequency_Up)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (offset: +1)") + 
      ylab("5' heterogeneity frequency (offset: -1)") + 
      coord_fixed()
    
    temp.data9 <- data1[!is.na(data1$Tsinghwa_Frequency_Down) & !is.na(data1$Tsinghwa_Frequency_Up),]
    
    p9 <- ggplot(temp.data9,aes(x=log10(Tsinghwa_Frequency_Down),y=log10(Tsinghwa_Frequency_Up))) +
      geom_point(size=1.5) + 
      geom_smooth(method=lm,color="lightslateblue",fill="lightslateblue") +
      geom_abline(slope = 1,linetype ="dashed") + 
      geom_abline(slope = 1, intercept = c(1,-1), linetype="dashed", color="red")+
      scale_x_continuous(breaks=c(-6,-4,-2,0),limits=c(-6,0)) + 
      scale_y_continuous(breaks=c(-6,-4,-2,0),limits=c(-6,0)) + 
      scale_color_manual(values=c("black","red","blue"),guide=F) +
      ggtitle(label=bquote("n ="~.(nrow(temp.data9))~";"~italic(rho)~"="~.(round(cor(temp.data9$Tsinghwa_Frequency_Down,temp.data9$Tsinghwa_Frequency_Up),2))~";"~italic(P)~"="~.(signif(cor.test(temp.data9$Tsinghwa_Frequency_Down,temp.data9$Tsinghwa_Frequency_Up)$p.value,digits=3)))) + 
      theme_bw() + 
      theme(axis.text.x=element_text(size=12), axis.title=element_text(size=15), axis.text.y=element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
      xlab("5' heterogeneity frequency (offset: +1)") + 
      ylab("5' heterogeneity frequency (offset: -1)") + 
      coord_fixed()
    p.offset <- ggarrange(p7,p8,p9,nrow=1,ncol=3,labels=c("Catholic","TCGA","Tsinghua"))
    ggsave(plot=p.offset,file=paste(inputD,"/Figure_Correlation/CorrelationAll_",arm,"_BetweenOffset_frequency.pdf",sep=""),dpi=300,width=10,height=4)
  }
}
