library(ggplot2)
library(ggrepel)

## isomiRs
data.cat.isomiR <- read.table("New/LiverCancer/Catholic/1.miRNA/IsomiR_position_information_v2.txt",header=T, row.names=NULL,check.names=F,quote="",sep="\t")
data.tcga.isomiR <- read.table("New/LiverCancer/TCGA/1.miRNA/IsomiR_position_information_v2.txt",header=T, row.names=NULL,check.names=F,quote="",sep="\t")
data.tsinghua.isomiR <- read.table("New/LiverCancer/Tsinghua/1.miRNA/IsomiR_position_information_v2.txt",header=T, row.names=NULL,check.names=F,quote="",sep="\t")

## Anova FDR <= 0.05 in all cohorts
data.cat <- read.table("New/LiverCancer/Catholic/1.miRNA/DE_isomiR_ratio.txt",header=T, row.names=NULL,check.names=F,quote="",sep="\t")
data.tcga <- read.table("New/LiverCancer/TCGA/1.miRNA/DE_isomiR_ratio.txt",header=T, row.names=NULL,check.names=F,quote="",sep="\t")
data.tsinghua <- read.table("New/LiverCancer/Tsinghua/1.miRNA/DE_isomiR_ratio.txt",header=T, row.names=NULL,check.names=F,quote="",sep="\t")

#head(data.cat$isomiR)
data.cat.filtered <- data.cat[data.cat$Anova <= 0.05,]
data.tcga.filtered <- data.tcga[data.tcga$Anova <= 0.05,]
data.tsinghua.filtered <- data.tsinghua[data.tsinghua$Anova <= 0.05,]

total <- intersect(data.cat.filtered$isomiR, data.tsinghua.filtered$isomiR)

### DE-mi abs(log2(fold change)) >= log2(1.5) and FDR <= 0.05 in Catholic & Tsinghua

data.cat <- read.table("New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_TotalDEG.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
data.tcga <- read.table("New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_TotalDEG.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
data.tsinghua <- read.table("New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered_TotalDEG.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
#head(data.cat)
data.cat$Source <- "Catholic"
data.tcga$Source <- "TCGA"
data.tsinghua$Source <- "Tsinghua"

data.cat$isomiR <- "N"
data.tcga$isomiR <- "N"
data.tsinghua$isomiR <- "N"


data.cat[data.cat$Gene %in% data.cat.isomiR[data.cat.isomiR$Position != 0,]$miRNA,]$isomiR <- "Y"
data.tcga[data.tcga$Gene %in% data.tcga.isomiR[data.tcga.isomiR$Position != 0,]$miRNA,]$isomiR <- "Y"
data.tsinghua[data.tsinghua$Gene %in% data.tsinghua.isomiR[data.tsinghua.isomiR$Position != 0,]$miRNA,]$isomiR <- "Y"

data.cat$Sig <- "ns"
data.cat[data.cat$Foldchange >= log2(1.5) & data.cat$FDR <= 0.05,]$Sig <- "sig"
data.cat[data.cat$Foldchange <= -log2(1.5) & data.cat$FDR <= 0.05,]$Sig <- "sig"

data.tcga$Sig <- "ns"
data.tcga[data.tcga$Foldchange >= log2(1.5) & data.tcga$FDR <= 0.05,]$Sig <- "sig"
data.tcga[data.tcga$Foldchange <= -log2(1.5) & data.tcga$FDR <= 0.05,]$Sig <- "sig"

data.tsinghua$Sig <- "ns"
data.tsinghua[data.tsinghua$Foldchange >= log2(1.5) & data.tsinghua$FDR <= 0.05,]$Sig <- "sig"
data.tsinghua[data.tsinghua$Foldchange <= -log2(1.5) & data.tsinghua$FDR <= 0.05,]$Sig <- "sig"

sigs2 <- intersect(data.cat[data.cat$Sig == "sig",]$Gene, intersect(data.tcga[data.tcga$Sig == "sig",]$Gene, data.tsinghua[data.tsinghua$Sig == "sig",]$Gene))
sigs1 <- intersect(data.cat[data.cat$Sig == "sig",]$Gene, data.tsinghua[data.tsinghua$Sig == "sig",]$Gene)

intersect(total, sigs2)

data <- rbind(data.cat, rbind(data.tcga, data.tsinghua))

data$Sig <- "ns"
data[data$Foldchange >= log2(1.5) & data$FDR <= 0.05,]$Sig <- "sig"
data[data$Foldchange <= -log2(1.5) & data$FDR <= 0.05,]$Sig <- "sig"
data[data$Gene %in% intersect(total, sigs2),]$Sig <- "sig2"
#data[data$Gene == "hsa-miR-21-5p_GCUUAUC_1",]$Sig <- "sig2"
#data[data$Gene == "hsa-miR-21-5p_UAGCUUA_-1",]$Sig <- "sig2"

data$Label <- ""
data[data$Gene %in% intersect(total, sigs2),]$Label <- as.character(data[data$Gene %in% intersect(total, sigs2),]$Gene)
#data[data$Gene == "hsa-miR-21-5p_GCUUAUC_1",]$Label <- "hsa-miR-21-5p_GCUUAUC_1"
#data[data$Gene == "hsa-miR-21-5p_UAGCUUA_-1",]$Label <- "hsa-miR-21-5p_UAGCUUA_-1"
#levels(factor(data$Label))
head(data)
data[data$Sig == "sig" & data$isomiR == "Y" ,]
data$Sig <- factor(data$Sig, levels=c("ns","sig","sig2"))
data$isomiR <- factor(data$isomiR, levels=c("N","Y"))
#nrow(data[data$Source == "Catholic",])
a <- ggplot(data, aes(x=Foldchange,y=-log10(FDR), color=Sig, shape=isomiR)) + 
  geom_point(size=1.5) +
  geom_text_repel(aes(label=Label),size=1,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines")) +
  scale_x_continuous(breaks=seq(-10,10,by=5),limits=c(-11,11)) + 
  #scale_y_continuous(breaks=seq(-6,2,by=2),limits=c(-6,2.5)) + 
  geom_hline(yintercept = -log10(0.05), color="blue", linetype = "dashed") +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color="blue", linetype="dashed") +
  scale_color_manual(values=c("gray","red","blue"),guide=F) +
  scale_shape_manual(values=c(32,16),guide=F) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=7), axis.title=element_text(size=7), axis.text.y=element_text(size=7), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank(), aspect.ratio = 1) + 
  xlab("Fold change (log2)") + 
  ylab("FDR (-log10)") + 
  facet_wrap(~Source, scales = "free_y")

ggsave(a,file="New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_TotalDEG.pdf",dpi=300,width=12,height=5)
