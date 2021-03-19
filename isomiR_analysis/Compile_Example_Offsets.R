library(ggplot2)
library(reshape2)
library(ggpubr)

mir.148b.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",10,20,14,90,80,85),ncol=3))
mir.148b.tp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",31,43,26,63,29,69),ncol=3))
mir.145.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",4,34,2,86,54,91),ncol=3))
mir.145.tp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",28,58,22,54,38,57),ncol=3))
mir.146b.tp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",25,8,27,69,90,66),ncol=3))
mir.146b.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",31,28,38,69,71,61),ncol=3))
mir.181c.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",0,2,0,100,98,100),ncol=3))
mir.181c.tp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",20,38,20,80,52,80),ncol=3))
mir.425.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",19,20,21,77,76,75),ncol=3))
mir.425.tp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",69,35,19,18,59,77),ncol=3))

mir.122.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",99.94,99.81,99.98,0.02,0.17,0.01),ncol=3))
mir.192.fp <- data.frame(matrix(c("Catholic","TCGA","Tsinghua",68.11,58.38,70.85,31.82,41.54,29.11),ncol=3))

draw.fig <- function(data, title){
  colnames(data) <- c("Source","red","blue")
  data.melt <- melt(data, id.vars = "Source")
  data.melt$Source <- factor(data.melt$Source, levels=c("Tsinghua","TCGA","Catholic"))
  data.melt$variable <- factor(data.melt$variable, levels=c("blue","red"))
  data.melt$value <- as.numeric(as.character(data.melt$value))
  
  ggplot(data.melt, aes(x=Source, y=value, fill=variable)) + 
    geom_bar(stat="identity") + 
    scale_fill_manual(values=c("blue","red"), name = "", labels=c("converted","original")) +
    scale_y_continuous(limits=c(0,100), name = "Read distribution (%)") + 
    theme_classic() + 
    theme(axis.text = element_text(size=7,color="black"), axis.title.x = element_blank(), axis.text.y = element_text(size=7, color="black"),panel.grid=element_blank(), legend.title = element_text(size=7, color="black"), legend.text = element_text(size=7,color="black"),aspect.ratio = 1, title = element_text(size=7,color="black")) + coord_flip() + ggtitle(title)
}

mir.148b.fp.fig <- draw.fig(mir.148b.fp,"mir.148b.fp")
mir.148b.tp.fig <- draw.fig(mir.148b.tp,"mir.148b.tp")
mir.145.fp.fig <- draw.fig(mir.145.fp,"mir.145.fp")
mir.145.tp.fig <- draw.fig(mir.145.tp,"mir.145.tp")
mir.146b.tp.fig <- draw.fig(mir.146b.tp,"mir.146b.tp")
mir.146b.fp.fig <- draw.fig(mir.146b.fp,"mir.146b.fp")
mir.181c.fp.fig <- draw.fig(mir.181c.fp,"mir.181c.fp")
mir.181c.tp.fig <- draw.fig(mir.181c.tp,"mir.181c.tp")
mir.425.fp.fig <- draw.fig(mir.425.fp,"mir.425.fp")
mir.425.tp.fig <- draw.fig(mir.425.tp,"mir.425.tp")

mir.122.fp.fig <- draw.fig(mir.122.fp,"mir.122.fp")
mir.192.fp.fig <- draw.fig(mir.192.fp,"mir.192.fp")

total <- ggarrange(mir.148b.fp.fig,mir.148b.tp.fig,mir.145.fp.fig,mir.145.tp.fig,mir.146b.tp.fig,mir.146b.fp.fig,mir.181c.fp.fig,mir.181c.tp.fig,mir.425.fp.fig,mir.425.tp.fig,mir.122.fp.fig,mir.192.fp.fig, nrow=6, ncol=2, common.legend = T)

ggsave(total,file="New/LiverCancer/Catholic/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/Example_isomiR_expression.pdf",dpi=300,width=5,height=10)
