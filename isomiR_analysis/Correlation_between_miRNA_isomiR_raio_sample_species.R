library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)
library(data.table)

dataset <- "TCGA_rev1" # Catholic TCGA_rev1 Tsinghua

miRNA <- read.table(paste("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA//Total/Frequency/AlleleFrequencyDistAll_Mature_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep='\t')
isomir <- read.table(paste("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1Minus1_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep="\t")
fold <- read.table(paste("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_All.txt",sep=""),header=T,row.names = 1,check.names=F,quote="",sep="\t")

miRNA.fp <- miRNA[rownames(miRNA) %like% "FP",]
miRNA.tp <- miRNA[rownames(miRNA) %like% "TP",]
isomir.fp <- isomir[rownames(isomir) %like% "FP",]
isomir.tp <- isomir[rownames(isomir) %like% "TP",]
fold.fp <- fold[rownames(fold) %like% "FP",]
fold.tp <- fold[rownames(fold) %like% "TP",]

mature <- log10(miRNA.fp+0.1)
iso <- log10(isomir.fp+0.1)
ratio <- fold.fp
set <- intersect(rownames(iso),rownames(ratio))
iso.filtered <- iso[rownames(iso) %in% set,]
ratio.filtered <- ratio[rownames(ratio) %in% set,]
mature.filtered <- mature[rownames(mature) %in% set,]

total.frame <- data.frame()
#head(iso.filtered)
iso.filtered.melt <- melt(t(iso.filtered))
ratio.filtered.melt <- melt(t(ratio.filtered))
mature.filtered.melt <- melt(t(mature.filtered))

#sum(iso.filtered.melt$Var1 != ratio.filtered.melt$Var1)
#sum(iso.filtered.melt$Var1 != mature.filtered.melt$Var1)
total.data <- iso.filtered.melt
total.data$ratio <- ratio.filtered.melt$value
total.data$mature <- mature.filtered.melt$value

colnames(total.data) <- c("Sample","miRNA","IsomiR","Ratio","Mature")
total.data.cat <- total.data

miRNA <- read.table(paste("Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA//Total/Frequency/AlleleFrequencyDistAll_Mature_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep='\t')
isomir <- read.table(paste("Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1Minus1_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep="\t")
fold <- read.table(paste("Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_All.txt",sep=""),header=T,row.names = 1,check.names=F,quote="",sep="\t")

miRNA.fp <- miRNA[rownames(miRNA) %like% "FP",]
miRNA.tp <- miRNA[rownames(miRNA) %like% "TP",]
isomir.fp <- isomir[rownames(isomir) %like% "FP",]
isomir.tp <- isomir[rownames(isomir) %like% "TP",]
fold.fp <- fold[rownames(fold) %like% "FP",]
fold.tp <- fold[rownames(fold) %like% "TP",]

mature <- log10(miRNA.fp+0.1)
iso <- log10(isomir.fp+0.1)
ratio <- fold.fp
set <- intersect(rownames(iso),rownames(ratio))
iso.filtered <- iso[rownames(iso) %in% set,]
ratio.filtered <- ratio[rownames(ratio) %in% set,]
mature.filtered <- mature[rownames(mature) %in% set,]

total.frame <- data.frame()
#head(iso.filtered)
iso.filtered.melt <- melt(t(iso.filtered))
ratio.filtered.melt <- melt(t(ratio.filtered))
mature.filtered.melt <- melt(t(mature.filtered))

#sum(iso.filtered.melt$Var1 != ratio.filtered.melt$Var1)
#sum(iso.filtered.melt$Var1 != mature.filtered.melt$Var1)
total.data <- iso.filtered.melt
total.data$ratio <- ratio.filtered.melt$value
total.data$mature <- mature.filtered.melt$value

colnames(total.data) <- c("Sample","miRNA","IsomiR","Ratio","Mature")
total.data.tcga <- total.data

miRNA <- read.table(paste("Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA//Total/Frequency/AlleleFrequencyDistAll_Mature_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep='\t')
isomir <- read.table(paste("Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlus1Minus1_RPM.txt",sep=""),header=T,row.names=1,check.names=F,quote="",sep="\t")
fold <- read.table(paste("Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_All.txt",sep=""),header=T,row.names = 1,check.names=F,quote="",sep="\t")

miRNA.fp <- miRNA[rownames(miRNA) %like% "FP",]
miRNA.tp <- miRNA[rownames(miRNA) %like% "TP",]
isomir.fp <- isomir[rownames(isomir) %like% "FP",]
isomir.tp <- isomir[rownames(isomir) %like% "TP",]
fold.fp <- fold[rownames(fold) %like% "FP",]
fold.tp <- fold[rownames(fold) %like% "TP",]

mature <- log10(miRNA.fp+0.1)
iso <- log10(isomir.fp+0.1)
ratio <- fold.fp
set <- intersect(rownames(iso),rownames(ratio))
iso.filtered <- iso[rownames(iso) %in% set,]
ratio.filtered <- ratio[rownames(ratio) %in% set,]
mature.filtered <- mature[rownames(mature) %in% set,]

total.frame <- data.frame()
#head(iso.filtered)
iso.filtered.melt <- melt(t(iso.filtered))
ratio.filtered.melt <- melt(t(ratio.filtered))
mature.filtered.melt <- melt(t(mature.filtered))

#sum(iso.filtered.melt$Var1 != ratio.filtered.melt$Var1)
#sum(iso.filtered.melt$Var1 != mature.filtered.melt$Var1)
total.data <- iso.filtered.melt
total.data$ratio <- ratio.filtered.melt$value
total.data$mature <- mature.filtered.melt$value

colnames(total.data) <- c("Sample","miRNA","IsomiR","Ratio","Mature")
total.data.tsinghua <- total.data

set <- intersect(intersect(total.data.cat$miRNA, total.data.tcga$miRNA), total.data.tsinghua$miRNA)

total.data.cat.filtered <- total.data.cat[total.data.cat$miRNA %in% set,]
total.data.tcga.filtered <- total.data.tcga[total.data.tcga$miRNA %in% set,]
total.data.tsinghua.filtered <- total.data.tsinghua[total.data.tsinghua$miRNA %in% set,]

total.data.cat.filtered$Source <- "Catholic"
total.data.tcga.filtered$Source <- "TCGA"
total.data.tsinghua.filtered$Source <- "Tsinghua"

total.data <- rbind(rbind(total.data.cat.filtered, total.data.tcga.filtered), total.data.tsinghua.filtered)
total.data <- rbind(total.data.cat.filtered, total.data.tsinghua.filtered)
total.data$Source <- factor(total.data$Source)
total.data$Sample <- as.numeric(total.data$Sample)
total.data$miRNA <- as.numeric(total.data$miRNA)
total.data$Source <- as.numeric(total.data$Source)
head(total.data)
temp.logit <- lm(Ratio ~ Sample + miRNA + IsomiR + Mature + Source, total.data)
summary(temp.logit)

temp.logit <- lm(IsomiR ~ Sample + miRNA + Ratio + Mature + Source, total.data)
summary(temp.logit)

temp.logit <- lm(IsomiR ~ Ratio + Mature + Source, total.data)
summary(temp.logit)

total.data.cat$Sample <- as.numeric(total.data.cat$Sample)
total.data.cat$miRNA <- as.numeric(total.data.cat$miRNA)
temp.logit.cat <- lm(IsomiR ~ Sample + miRNA + Ratio + Mature, total.data.cat)
summary(temp.logit.cat)

temp.logit.cat <- lm(IsomiR ~ Ratio + Mature, total.data.cat)
summary(temp.logit.cat)

total.data.tcga$Sample <- as.numeric(total.data.tcga$Sample)
total.data.tcga$miRNA <- as.numeric(total.data.tcga$miRNA)
temp.logit.tcga <- lm(IsomiR ~ Sample + miRNA + Ratio + Mature, total.data.tcga)
summary(temp.logit.tcga)

temp.logit.tcga <- lm(IsomiR ~ Ratio + Mature, total.data.tcga)
summary(temp.logit.tcga)

total.data.tsinghua$Sample <- as.numeric(total.data.tsinghua$Sample)
total.data.tsinghua$miRNA <- as.numeric(total.data.tsinghua$miRNA)
temp.logit.tsinghua <- lm(IsomiR ~ Sample + miRNA + Ratio + Mature, total.data.tsinghua)
summary(temp.logit.tsinghua)

temp.logit.tsinghua <- lm(IsomiR ~ Ratio + Mature, total.data.tsinghua)
summary(temp.logit.tsinghua)
