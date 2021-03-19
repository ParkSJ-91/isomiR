library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)
library(data.table)

data.cat <- read.table("New/LiverCancer/Catholic/1.miRNA/Correlation_Mature_AllIsomiRs_AllSamples.txt",header=T,row.names=NULL,check.names=F,quote="",sep='\t')
data.tcga <- read.table("New/LiverCancer/TCGA/1.miRNA/Correlation_Mature_AllIsomiRs_AllSamples.txt",header=T,row.names=NULL,check.names=F,quote="",sep='\t')
data.tsinghua <- read.table("New/LiverCancer/Tsinghua/1.miRNA/Correlation_Mature_AllIsomiRs_AllSamples.txt",header=T,row.names=NULL,check.names=F,quote="",sep='\t')

#data.cat$Ratio <- 10**data.cat$Ratio
#data.tcga$Ratio <- 10**data.tcga$Ratio
#data.tsinghua$Ratio <- 10**data.tsinghua$Ratio

set <- intersect(data.cat$IsomiR, data.tsinghua$IsomiR)
#head(total.data.cat.filtered)
total.data.cat.filtered <- data.cat[data.cat$IsomiR %in% set,]
total.data.tcga.filtered <- data.tcga[data.tcga$IsomiR %in% set,]
total.data.tsinghua.filtered <- data.tsinghua[data.tsinghua$IsomiR %in% set,]

#total.data <- rbind(rbind(total.data.cat.filtered, total.data.tcga.filtered), total.data.tsinghua.filtered)
total.data <- rbind(total.data.cat.filtered, total.data.tsinghua.filtered)
total.data$Source <- factor(total.data$Source)
total.data$Arm <- factor(total.data$Arm)

total.data$Sample <- as.numeric(total.data$Sample)
total.data$Precursor <- as.numeric(total.data$Precursor)
total.data$Source <- as.numeric(total.data$Source)
total.data$Arm <- as.numeric(total.data$Arm)
head(total.data)

total.data$IsomiRExpNorm <- (total.data$IsomiRExp - min(total.data$IsomiRExp)) / (max(total.data$IsomiRExp) - min(total.data$IsomiRExp))
total.data$RatioNorm <- (total.data$Ratio - min(total.data$Ratio)) / (max(total.data$Ratio) - min(total.data$Ratio))
total.data$PreExpNorm <- (total.data$PreExp - min(total.data$PreExp)) / (max(total.data$PreExp) - min(total.data$PreExp))

temp.logit <- lm(Ratio ~ Sample + Precursor + Arm + Offset + IsomiRExp + PreExp + Source, total.data)
#temp.logit <- lm(Ratio ~ Sample + Precursor + log10(IsomiRExp+0.1) + log10(CanoExp+0.1) + Source, total.data)
summary(temp.logit)$coefficients
#levels(factor(total.data$Source))
#temp.logit <- lm(IsomiRExp ~ Sample + Precursor + Arm + Offset + Ratio + CanoExp + Source, total.data)
#summary(temp.logit)

temp.logit <- lm(IsomiRExp ~ Ratio + PreExp + Source, total.data)
summary(temp.logit)$coefficients

data.cat$Sample <- as.numeric(data.cat$Sample)
data.cat$Precursor <- as.numeric(data.cat$Precursor)
data.cat$IsomiRExpNorm <- (data.cat$IsomiRExp - min(data.cat$IsomiRExp)) / (max(data.cat$IsomiRExp) - min(data.cat$IsomiRExp))
data.cat$RatioNorm <- (data.cat$Ratio - min(data.cat$Ratio)) / (max(data.cat$Ratio) - min(data.cat$Ratio))
data.cat$PreExpNorm <- (data.cat$PreExp - min(data.cat$PreExp)) / (max(data.cat$PreExp) - min(data.cat$PreExp))

temp.logit.cat <- lm(Ratio ~ Sample + Precursor + Arm + Offset + IsomiRExp + PreExp, data.cat)
#temp.logit.cat <- lm(IsomiRExp ~ Sample + Precursor + Ratio + CanoExp, data.cat)
summary(temp.logit.cat)$coefficients

temp.logit.cat <- lm(IsomiRExpNorm ~ RatioNorm + PreExpNorm, data.cat)
summary(temp.logit.cat)$coefficients

data.tcga$Sample <- as.numeric(data.tcga$Sample)
data.tcga$Precursor <- as.numeric(data.tcga$Precursor)
data.tcga$IsomiRExpNorm <- (data.tcga$IsomiRExp - min(data.tcga$IsomiRExp)) / (max(data.tcga$IsomiRExp) - min(data.tcga$IsomiRExp))
data.tcga$RatioNorm <- (data.tcga$Ratio - min(data.tcga$Ratio)) / (max(data.tcga$Ratio) - min(data.tcga$Ratio))
data.tcga$PreExpNorm <- (data.tcga$PreExp - min(data.tcga$PreExp)) / (max(data.tcga$PreExp) - min(data.tcga$PreExp))

temp.logit.tcga <- lm(Ratio ~ Sample + Precursor + Arm + Offset + IsomiRExp + PreExp, data.tcga)
#temp.logit.tcga <- lm(IsomiRExp ~ Sample + Precursor + Arm + Offset + Ratio + CanoExp, temp.logit.tcga)
summary(temp.logit.tcga)$coefficients

#nrow(data.tcga)
temp.logit.tcga <- lm(IsomiRExpNorm ~ RatioNorm + PreExpNorm, data.tcga)
summary(temp.logit.tcga)$coefficients

data.tsinghua$Sample <- as.numeric(data.tsinghua$Sample)
data.tsinghua$Precursor <- as.numeric(data.tsinghua$Precursor)
data.tsinghua$IsomiRExpNorm <- (data.tsinghua$IsomiRExp - min(data.tsinghua$IsomiRExp)) / (max(data.tsinghua$IsomiRExp) - min(data.tsinghua$IsomiRExp))
data.tsinghua$RatioNorm <- (data.tsinghua$Ratio - min(data.tsinghua$Ratio)) / (max(data.tsinghua$Ratio) - min(data.tsinghua$Ratio))
data.tsinghua$PreExpNorm <- (data.tsinghua$PreExp - min(data.tsinghua$PreExp)) / (max(data.tsinghua$PreExp) - min(data.tsinghua$PreExp))

temp.logit.tsinghua <- lm(Ratio ~ Sample + Precursor + Arm + Offset + IsomiRExp + PreExp, data.tsinghua)
summary(temp.logit.tsinghua)$coefficients

temp.logit.tsinghua <- lm(IsomiRExpNorm ~ RatioNorm + PreExpNorm, data.tsinghua)
summary(temp.logit.tsinghua)$coefficients
 