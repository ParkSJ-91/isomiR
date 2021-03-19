library(ggplot2)
library(ggrepel)

shift_axis <- function(p, y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(y=y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())
}

data <- read.table("Project/LiverCancer/Catholic/3.TargetAnalysis/9.TargetAnalysis_noPVTT/TotalCorrelation_AllTarget_GeneLevel_include6mer/GO_Cancer_CS-0.1_IsomiR.txt",header=T, row.names=NULL, sep='\t',check.names=F)
inputD <- "New/LiverCancer/Catholic/3.TargetingAnalysis/PlotTargetSites/"
outputD <- "New/LiverCancer/Catholic/3.TargetingAnalysis/PlotTargetSites/"
geneSymbol <- "GHR"
miRNA <- "hsa-miR-21-5p_GCUUAUC_1"
for (i in 1:nrow(data)){
  geneSymbol <- data[i,]$GeneSymbol
  miRNA <- data[i,]$miRNA
  airD = read.table(paste(inputD,geneSymbol,"_",miRNA,"_air.txt",sep=""),header=T,row.names=NULL,check.names=F,sep='\t')
  siteD = read.table(paste(inputD,geneSymbol,"_",miRNA,"_site.txt",sep=""),header=T,row.names=NULL,check.names=F,sep='\t')
  newList <- c()
  for (type in factor(siteD$type)){
    if (type == "7mer-1a"){
      newList <- c(newList, "green")
      #siteD[i,]$color = "#009640"
    } else if (type == "7mer-m8"){
      newList <- c(newList, "cyan")
      #siteD[i,]$color = "#009ee3"
    } else if (type == "6mer"){
      newList <- c(newList, "orange")
      #siteD[i,]$color = "#ffed00"
    } else if (type == "8mer-1a"){
      newList <- c(newList, "blue")
      #siteD[i,]$color = "#302683"
    }
  }
  g <- ggplot() + geom_step(data=airD, aes(x=x,y=y),color="red",size=2)+ theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(size=25), axis.ticks.x = element_line(size=2),axis.ticks.length = unit(0.3,"cm")) + coord_cartesian(xlim=c(0,airD[nrow(airD),]$x)) + scale_x_continuous(expand=c(0,0)) + ylim(-1,1) 
  g <- g + geom_path(aes(x=c(0,airD[nrow(airD),]$x),y=0),color="grey")
  g <- g + geom_path(aes(x=c(0,airD[nrow(airD),]$x),y=0.25),color="grey")
  g <- g + geom_path(aes(x=c(0,airD[nrow(airD),]$x),y=0.5),color="grey")
  g <- g + geom_path(aes(x=c(0,airD[nrow(airD),]$x),y=0.75),color="grey")
  g <- g + geom_path(aes(x=c(0,airD[nrow(airD),]$x),y=1),color="grey")
  g <- g + geom_segment(data=siteD, aes(x=start,xend=end, y = -.7, yend = -.7, color=type), size=8) + scale_color_manual(values=newList) + geom_text_repel(data=siteD, aes(x = (start + end)/2 , y = -0.9,label = miRNA), nudge_y=0,point.padding=NA,size=7)
  g <- g + theme(legend.position ="none")
  newG <- shift_axis(g,-.2)
  ggsave(file=paste(outputD,geneSymbol,"_",miRNA,".pdf",sep=""),newG, dpi=300,width=30,height=3,unit="in")
  #ggsave(file=paste(outputD,geneSymbol,"_",miRNA,".png",sep=""),newG, dpi=300,width=30,height=3,unit="in")
}
