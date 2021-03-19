library(rsgcc)
library(dendsort)
library(reshape2)
library(ggplot2)

data.cat <- read.table("New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_MatchedMature.txt",header=T,row.names=1,check.names=F,sep='\t',quote="")
output.cat <- "New/LiverCancer/Catholic/1.miRNA/Clustering/MatchedMature.pdf"
info.cat <- read.table("Project/LiverCancer/reference/Catholic_info_v3.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")

label.cat1 <- c("Normal","FL","FH","CS","DL","DH")
label.cat2 <- c("G1","G2","G3")
sample.order.cat <- c()
breaks.cat <- c()
for (label in label.cat1){
  temp.data <- info.cat[info.cat$Grade == label,]
  sorted.samples <- sort(as.character(temp.data$sRNAseq_ID))
  sample.order.cat <- c(sample.order.cat, sorted.samples)
  breaks.cat <- c(breaks.cat,sorted.samples[length(sorted.samples)])
}
#class(label)
info.cat.filtered <- info.cat[! info.cat$Grade %in% label.cat1,]
for (label in label.cat2){
  temp.data <- info.cat.filtered[as.character(info.cat.filtered$neoplasm_histologic_grade) == label,]
  sorted.samples <- sort(as.character(temp.data$sRNAseq_ID))
  sample.order.cat <- c(sample.order.cat, sorted.samples)
  breaks.cat <- c(breaks.cat,sorted.samples[length(sorted.samples)])
}
temp.data <- info.cat.filtered[is.na(info.cat.filtered$neoplasm_histologic_grade),]
sorted.samples <- sort(as.character(temp.data$sRNAseq_ID))
sample.order.cat <- c(sample.order.cat, sorted.samples)
breaks.cat <- c(breaks.cat,sorted.samples[length(sorted.samples)])
label.cat <- c(label.cat1,label.cat2,"NA")

order.file <- read.table("New/LiverCancer/Catholic/1.miRNA/Clustering/Cluster_info.txt",header=F,row.names=NULL,sep='\t',check.names=F,quote="")
order.file <- subset(order.file, V2 != 0)

# adjust cluster order
levels(order.file$V2)
#cluster.level = c("Cluster_4","Cluster_5","Cluster_1","Cluster_2","Cluster_3","Cluster_6")
cluster.level = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5")
cluster.level = c("Cluster_5","Cluster_4","Cluster_1","Cluster_3","Cluster_2")
order.file <- order.file[order.file$V2 %in% cluster.level,]
order.file$V2 <- factor(order.file$V2, levels = cluster.level)
order.file1 <- order.file[order(order.file$V2),]
gene.order <- as.character(order.file1$V1)

data.cat1 <- subset(data.cat, row.names(data.cat) %in% gene.order)
data.cat2 <- scale(t(data.cat1))
data.cat3 <- melt(data.cat2)
colnames(data.cat3) <- c("Sample","miRNA","Zscore")
data.cat3$miRNA <- factor(data.cat3$miRNA, levels = gene.order)
data.cat3$Sample <- factor(data.cat3$Sample, levels = sample.order.cat)
data.cat3$Zscore[data.cat3$Zscore < -2] <- -2
data.cat3$Zscore[data.cat3$Zscore > 2] <- 2
data.cat3[is.na(data.cat3$Zscore),]$Zscore <- 0
plot.cat <- ggplot(data.cat3,aes(x=Sample,y=miRNA)) + geom_tile(aes(fill=Zscore)) + scale_fill_gradient2(low="darkgoldenrod3",high="blue",mid="white",midpoint = 0) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=12), axis.ticks.y=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank()) + scale_x_discrete(breaks = breaks.cat, label = label.cat) + coord_fixed()
ggsave(output.cat,plot.cat,dpi=300,width=10,height=10)


data.tcga <- read.table("New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_MatchedIsomiR.txt",header=T,row.names=1,check.names=F,sep='\t',quote="")
output.tcga <- "New/LiverCancer/TCGA/1.miRNA/Clustering/MatchedIsomiR.pdf"

info.tcga <- read.table("Project/LiverCancer/src_rev1/SurvivalAnalysis/renew2/TCGA_New.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
label.tcga <- c("G1","G2","G3+")
sample.order.tcga <- c()
breaks.tcga <- c()
for (label in label.tcga){
  temp.data <- info.tcga[info.tcga$Grade == label,]
  sorted.samples <- sort(as.character(temp.data$RNAseq_ID))
  sample.order.tcga <- c(sample.order.tcga,sorted.samples)
  breaks.tcga <- c(breaks.tcga,sorted.samples[length(sorted.samples)])
}
sorted.samples <- sort(colnames(data.tcga)[! colnames(data.tcga) %in% info.tcga$RNAseq_ID])
sample.order.tcga <- c(sorted.samples,sample.order.tcga)
breaks.tcga <- c(sorted.samples[length(sorted.samples)],breaks.tcga)
label.tcga <- c("Normal",label.tcga)

order.file <- read.table("New/LiverCancer/TCGA/1.miRNA/Clustering/Cluster_info.txt",header=F,row.names=NULL,sep='\t',check.names=F,quote="")
order.file <- subset(order.file, V2 != 0)

# adjust cluster order
levels(order.file$V2)
#cluster.level = c("Cluster_4","Cluster_5","Cluster_1","Cluster_2","Cluster_3","Cluster_6")
cluster.level = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5")
cluster.level = c("Cluster_5","Cluster_3","Cluster_2","Cluster_4","Cluster_1")
order.file <- order.file[order.file$V2 %in% cluster.level,]
order.file$V2 <- factor(order.file$V2, levels = cluster.level)
order.file1 <- order.file[order(order.file$V2),]
gene.order <- as.character(order.file1$V1)

data.tcga1 <- subset(data.tcga, rownames(data.tcga) %in% gene.order)
data.tcga2 <- scale(t(data.tcga1))
data.tcga3 <- melt(data.tcga2)
colnames(data.tcga3) <- c("Sample","miRNA","Zscore")
data.tcga3$miRNA <- factor(data.tcga3$miRNA, levels = gene.order)
data.tcga3$Sample <- factor(data.tcga3$Sample, levels = sample.order.tcga)
data.tcga3$Zscore[data.tcga3$Zscore < -2] <- -2
data.tcga3$Zscore[data.tcga3$Zscore > 2] <- 2
data.tcga3[is.na(data.tcga3$Zscore),]$Zscore <- 0
plot.tcga <- ggplot(data.tcga3,aes(x=Sample,y=miRNA)) + geom_tile(aes(fill=Zscore)) + scale_fill_gradient2(low="darkgoldenrod3",high="blue",mid="white",midpoint = 0) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=12), axis.ticks.y=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank()) + scale_x_discrete(breaks = breaks.tcga, label = label.tcga) + coord_fixed()
ggsave(output.tcga,plot.tcga,dpi=300,width=10,height=10)


data.tsinghwa <- read.table("New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered_MatchedMature.txt", header=T,row.names=1,check.names=F,sep='\t',quote="")
output.tsinghwa <- "New/LiverCancer/Tsinghua/1.miRNA/MatchedMature.pdf"
head(data.tsinghwa)

info.tsinghwa <- read.table("Project/LiverCancer/Tsinghua/info/Match_GSM_SRR_Sample.txt",header=T,row.names=NULL,check.names=F,quote="",sep='\t')
label.tsinghwa <- c("Normal","Cancer")
sample.order.tsinghwa <- c()
breaks.tsinghwa <- c()
for (label in label.tsinghwa){
  temp.data <- info.tsinghwa[info.tsinghwa$Grade == label,]
  sorted.samples <- sort(as.character(temp.data$run_accession))
  sample.order.tsinghwa <- c(sample.order.tsinghwa,sorted.samples)
  breaks.tsinghwa <- c(breaks.tsinghwa,sorted.samples[length(sorted.samples)])
}
order.file <- read.table("New/LiverCancer/Tsinghua/1.miRNA/Clustering/Cluster_info.txt",header=F,row.names=NULL,sep='\t',check.names=F,quote="")
order.file <- subset(order.file, V2 != 0)

# adjust cluster order
levels(order.file$V2)
#cluster.level = c("Cluster_4","Cluster_5","Cluster_1","Cluster_2","Cluster_3","Cluster_6")
cluster.level = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4")
cluster.level = c("Cluster_4","Cluster_3","Cluster_2","Cluster_1")
order.file <- order.file[order.file$V2 %in% cluster.level,]
order.file$V2 <- factor(order.file$V2, levels = cluster.level)
order.file1 <- order.file[order(order.file$V2),]
gene.order <- as.character(order.file1$V1)

data.tsinghwa1 <- subset(data.tsinghwa, rownames(data.tsinghwa) %in% gene.order)
data.tsinghwa2 <- scale(t(data.tsinghwa1))
data.tsinghwa3 <- melt(data.tsinghwa2)
colnames(data.tsinghwa3) <- c("Sample","miRNA","Zscore")
data.tsinghwa3$miRNA <- factor(data.tsinghwa3$miRNA, levels = gene.order)
data.tsinghwa3$Sample <- factor(data.tsinghwa3$Sample, levels = sample.order.tsinghwa)
data.tsinghwa3$Zscore[data.tsinghwa3$Zscore < -2] <- -2
data.tsinghwa3$Zscore[data.tsinghwa3$Zscore > 2] <- 2

data.tsinghwa3[is.na(data.tsinghwa3$Zscore),]$Zscore <- 0
plot.tsinghwa <- ggplot(data.tsinghwa3,aes(x=Sample,y=miRNA)) + geom_tile(aes(fill=Zscore)) + scale_fill_gradient2(low="darkgoldenrod3",high="blue",mid="white",midpoint = 0) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=12), axis.ticks.y=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_blank()) + scale_x_discrete(breaks = breaks.tsinghwa, label = label.tsinghwa) + coord_fixed()
ggsave(output.tsinghwa,plot.tsinghwa,dpi=300,width=10,height=10)
