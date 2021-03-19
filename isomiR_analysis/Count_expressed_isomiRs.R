library(ggplot2)
library(reshape2)
library(scales)

data.cat <- read.table("New/LiverCancer/Catholic/1.miRNA/IsomiR_count_v2.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
data.tcga <- read.table("New/LiverCancer/TCGA/1.miRNA/IsomiR_count_v2.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
data.tsinghua <- read.table("New/LiverCancer/Tsinghua/1.miRNA/IsomiR_count_v2.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")

head(data.cat)
data.cat$Data <- "Catholic"
data.tcga$Data <- "TCGA"
data.tsinghua$Data <- "Tsinghua"

data <- rbind(data.cat, rbind(data.tcga, data.tsinghua))

data.melt <- melt(data, id.vars = c("Arm","Data","Order"))

head(data.melt)
data.melt$variable <- factor(data.melt$variable, levels=seq(-5,5,by=1))
plot.count <- ggplot(data.melt,aes(x=variable,y=value)) + geom_bar(stat="identity") +
  scale_y_continuous(breaks=c(0,70,140,210,280,350)) +
  facet_grid(Arm~Data) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(), axis.text=element_text(size=12,color="black"), axis.title = element_text(size=12,color="black")) + 
  xlab("Offset") + ylab("Count")

plot.count <- ggplot(data.melt[data.melt$Order == "All",],aes(x=variable,y=value)) + geom_bar(stat="identity") +
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600,700),limits=c(0,700)) +
  facet_grid(~Data) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(), axis.text=element_text(size=12,color="black"), axis.title = element_text(size=12,color="black")) + 
  xlab("Offset") + ylab("Count")

data.melt$Order <- factor(data.melt$Order, levels=c("All","Cano","Major_1","Major_2","Others"))
plot.count <- ggplot(data.melt[data.melt$Order != "All",],aes(x=variable,y=value,fill=Order)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#fed8ce","#8f3a39","#5d7db0","gray"),name="") + 
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600,700),limits=c(0,700)) +
  facet_grid(~Data) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(), axis.text=element_text(size=12,color="black"), axis.title = element_text(size=12,color="black")) + 
  xlab("Offset") + ylab("Count")

ggsave(plot.count, file="New/LiverCancer/Catholic/1.miRNA/Count_expressed_isomiRs_total_v2.pdf",dpi=300,width=10,height=10)

data <- read.table("New/LiverCancer/Tsinghua/1.miRNA/IsomiR_position_information_v2.txt",header=T,row.names=NULL,check.names=F,quote="",sep="\t")
head(data)

data$miRNA <- factor(data$miRNA, levels=data$miRNA[rev(order(data$MedianExp))])
data$Order <- factor(data.melt$Order, levels=c("Cano","Major_1","Major_2","Others"))

data$Group <- "miRNAs"
data[data$Position == -1 | data$Position == 1,]$Group <- "isomiRs"
data[! data$Position %in% c(-1,0,1),]$Group <- "Others"

data$Group <- factor(data$Group, levels=c("miRNAs","isomiRs","Others"))
#data <- droplevels(data[data$Group != "Others",])
plot.adundance <- ggplot(data,aes(x=miRNA, y=log10(MedianExp+0.1)+1, fill=Order)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#fed8ce","#8f3a39","#5d7db0","gray"),name="") + 
  scale_y_continuous(limits=c(0,6),breaks=seq(0,6,1),labels=seq(0,6,1)-1) +
  theme_bw() + 
  theme(axis.title=element_text(size=12,color="black"),axis.text.y=element_text(size=12, color="black"), axis.text.x=element_blank(), panel.grid = element_blank(), legend.text = element_text(size=12, color="black"), axis.ticks.x = element_blank()) + 
  xlab("miRNAs") + ylab("RPM (log10)")
nrow(data)
data.num <- data.frame(table(data$Order))

data.num$Freq <- data.num$Freq / sum(data.num$Freq) * 100

data.num$Var1 <- factor(data.num$Var1, levels=c("Cano","Major_1","Major_2","Others"))
plot.pie <- ggplot(data.num,aes(x="",y=Freq,fill=Var1)) + 
  geom_bar(width=1,stat="identity") + 
  coord_polar("y",start=0) + 
  scale_fill_manual(values=c("#fed8ce","#8f3a39","#5d7db0","gray"),guide=F) +
  theme_minimal() +
  theme(axis.text=element_blank(), axis.title=element_blank(), panel.grid = element_blank(), panel.border=element_blank(), axis.ticks = element_blank()) +
  geom_text(aes(y = Freq + c(0, cumsum(Freq)[-length(Freq)]), 
                label = Freq/100), size=5)

ggsave(plot.adundance, file="New/LiverCancer/Tsinghua/1.miRNA/Abundance_expressed_isomiRs_pmall_v2.pdf",dpi=300,width=10,height=10)
ggsave(plot.pie, file="New/LiverCancer/Tsinghua/1.miRNA/Abundance_expressed_isomiRs_pmall_pie_v2.pdf",dpi=300,width=10,height=10)