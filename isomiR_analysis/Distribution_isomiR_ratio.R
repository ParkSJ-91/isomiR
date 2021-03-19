library(ggplot2)
library(ggrepel)
library(ggpubr)

outputD <- "New/LiverCancer/Catholic/1.miRNA/Figure_Correlation/"
dir.create(file.path(outputD), showWarnings = FALSE)

data.cat <- read.table(file="New/LiverCancer/Catholic/1.miRNA/Correlation_Mature_MajorIsomiR.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
data.tcga <- read.table(file="New/LiverCancer/TCGA/1.miRNA/Correlation_Mature_MajorIsomiR.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")
data.tsinghua <- read.table(file="New/LiverCancer/Tsinghua/1.miRNA/Correlation_Mature_MajorIsomiR.txt", header=T,row.names=NULL,check.names=F,sep='\t',quote="")

colnames(data.cat)

data.cat <- data.cat[! is.na(data.cat$Median_Ratio_Major_1),]
data.tcga <- data.tcga[! is.na(data.tcga$Median_Ratio_Major_1),]
data.tsinghua <- data.tsinghua[! is.na(data.tsinghua$Median_Ratio_Major_1),]

data.cat$Source <- "Catholic"
data.tcga$Source <- "TCGA"
data.tsinghua$Source <- "Tsinghua"

data <- rbind(data.cat, rbind(data.tcga, data.tsinghua))
data$Source <- factor(data$Source, levels=c("Catholic","TCGA","Tsinghua"))

a <- ggplot(data, aes(x=Median_Ratio_Major_1,y=..density..)) + 
  geom_histogram(position="stack",na.rm=T, color="black",fill="darkgreen") + 
  geom_density(alpha=0.2, na.rm=T, fill="#FF6666") +
  theme_bw() + theme(axis.text=element_text(size=7,color="black"), axis.title=element_text(size=7,color="black"), panel.grid=element_blank(), legend.text=element_text(size=7,color="black"), legend.title=element_text(size=7,color="black"),aspect.ratio = 1) + 
  facet_wrap(~Source)

ggsave(a, file=paste(outputD,"Distribution_IsomiR_Ratio.pdf",sep=""),dpi=300,width=12,height=5)
