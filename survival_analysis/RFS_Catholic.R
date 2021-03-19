library(survival)
library(survminer)

inputFile = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Qnorm_RPKMFiltered_Intersect.txt'
outputD = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis_v2/'

InfoTable <- read.table('/home/seokju/Project/LiverCancer/src_rev1/SurvivalAnalysis/renew2/Catholic_New.txt',header=T,row.names=3,sep='\t',check.names=F)
covariates_ori <- c("neoplasm_histologic_grade","AFP","Age","Gender","VascularInvasion","Alcohol","HBV","HCV","T_stage","OP_method")
InfoTable.filtered <- InfoTable[,colnames(InfoTable) %in% c(covariates_ori,"RFS","Recur")]
InfoTable.filtered <- InfoTable.filtered[InfoTable.filtered$RFS >= 0 & ! is.na(InfoTable.filtered$RFS) & InfoTable.filtered$RFS <= 60,]
InfoTable.uni.filtered <- InfoTable[,colnames(InfoTable) %in% c("RFS","Recur")]
InfoTable.uni.filtered <- InfoTable.uni.filtered[InfoTable.uni.filtered$RFS >= 0 & ! is.na(InfoTable.uni.filtered$RFS) InfoTable.uni.filtered$RFS <= 60,]
InfoTable.uni.filtered <- InfoTable.uni.filtered[complete.cases(InfoTable.uni.filtered),]

expTable <- read.table(inputFile,header=T,row.names=1,sep='\t',check.names = F)
for (geneSymbol in rownames(expTable)){
  covariates <- covariates_ori
  print(geneSymbol)
  print(covariates)
  temp.exp <- expTable[geneSymbol,]
  temp.exp.t.filtered <- t(temp.exp)
  rownames(temp.exp.t.filtered)
  
  temp.exp.t.filtered.filtered.uni <- data.frame(subset(temp.exp.t.filtered, rownames(temp.exp.t.filtered) %in% rownames(InfoTable.uni.filtered)),check.names=F)
  temp.exp.t.filtered.filtered.uni$Exp <- ""
  
  if (median(temp.exp.t.filtered.filtered.uni[,1])==max(temp.exp.t.filtered.filtered.uni[,1])){
    temp.exp.t.filtered.filtered.uni$Exp[temp.exp.t.filtered.filtered.uni[,1] >= median(temp.exp.t.filtered.filtered.uni[,1])] <- "Upper"
    temp.exp.t.filtered.filtered.uni$Exp[temp.exp.t.filtered.filtered.uni[,1] < median(temp.exp.t.filtered.filtered.uni[,1])] <- "Lower"
  } else if (median(temp.exp.t.filtered.filtered.uni[,1]) == min(temp.exp.t.filtered.filtered.uni[,1])){
    temp.exp.t.filtered.filtered.uni$Exp[temp.exp.t.filtered.filtered.uni[,1] > median(temp.exp.t.filtered.filtered.uni[,1])] <- "Upper"
    temp.exp.t.filtered.filtered.uni$Exp[temp.exp.t.filtered.filtered.uni[,1] <= median(temp.exp.t.filtered.filtered.uni[,1])] <- "Lower"
  } else {
    temp.exp.t.filtered.filtered.uni$Exp[temp.exp.t.filtered.filtered.uni[,1] >= median(temp.exp.t.filtered.filtered.uni[,1])] <- "Upper"
    temp.exp.t.filtered.filtered.uni$Exp[temp.exp.t.filtered.filtered.uni[,1] < median(temp.exp.t.filtered.filtered.uni[,1])] <- "Lower"
  }
  InfoTable.result.uni <- cbind(InfoTable.uni.filtered, temp.exp.t.filtered.filtered.uni[,"Exp"][match(rownames(InfoTable.uni.filtered), rownames(temp.exp.t.filtered.filtered.uni))])
  colnames(InfoTable.result.uni) <- c(colnames(InfoTable.uni.filtered),"Exp")
  InfoTable.result.uni <- droplevels(InfoTable.result.uni)
  
  fitCPH.uni <- coxph(Surv(RFS,Recur) ~ Exp, data=InfoTable.result.uni)
  fitCPHsumm.uni <- summary(fitCPH.uni)
  summaryTable.uni <- data.frame(fitCPHsumm.uni$coefficients, check.names = F)
  validation.uni <- cox.zph(fitCPH.uni)
  validationTable.uni <- validation.uni$table
  p.value.uni <- fitCPHsumm.uni$sctest["pvalue"]
  summaryTable.uni["P(LogRank)"] <- p.value.uni
  
  
  summary <- ""
  zph <- ""
  for (i in 1:length(covariates)){
    InfoTable.filtered <- InfoTable.filtered[complete.cases(InfoTable.filtered),]
    temp.exp.t.filtered.filtered <- data.frame(subset(temp.exp.t.filtered, rownames(temp.exp.t.filtered) %in% rownames(InfoTable.filtered)),check.names=F)
    temp.exp.t.filtered.filtered$Exp <- ""
    
    if (median(temp.exp.t.filtered.filtered[,1]) == max(temp.exp.t.filtered.filtered[,1])){
      temp.exp.t.filtered.filtered$Exp[temp.exp.t.filtered.filtered[,1] >= median(temp.exp.t.filtered.filtered[,1])] <- "Upper"
      temp.exp.t.filtered.filtered$Exp[temp.exp.t.filtered.filtered[,1] < median(temp.exp.t.filtered.filtered[,1])] <- "Lower"
    } else if (median(temp.exp.t.filtered.filtered[,1]) == min(temp.exp.t.filtered.filtered[,1])){
      temp.exp.t.filtered.filtered$Exp[temp.exp.t.filtered.filtered[,1] > median(temp.exp.t.filtered.filtered[,1])] <- "Upper"
      temp.exp.t.filtered.filtered$Exp[temp.exp.t.filtered.filtered[,1] <= median(temp.exp.t.filtered.filtered[,1])] <- "Lower"
    } else {
      temp.exp.t.filtered.filtered$Exp[temp.exp.t.filtered.filtered[,1] >= median(temp.exp.t.filtered.filtered[,1])] <- "Upper"
      temp.exp.t.filtered.filtered$Exp[temp.exp.t.filtered.filtered[,1] < median(temp.exp.t.filtered.filtered[,1])] <- "Lower"
    }
    InfoTable.result <- cbind(InfoTable.filtered, temp.exp.t.filtered.filtered[,"Exp"][match(rownames(InfoTable.filtered), rownames(temp.exp.t.filtered.filtered))])
    colnames(InfoTable.result) <- c(colnames(InfoTable.filtered),"Exp")
    InfoTable.result <- droplevels(InfoTable.result)
    
    temp.formular <- as.formula(paste("Surv(RFS,Recur) ~ Exp + ",paste(covariates,collapse=" + "),sep=""))
    temp.model <- coxph(temp.formular,data=InfoTable.result)
    temp.zph <- cox.zph(temp.model)
    temp.zph.frame <- data.frame(temp.zph$table)
    dangers <- rownames(temp.zph.frame[temp.zph.frame$p <= 0.05,])
    danger.covariates <- c()
    for (temp.covariate in covariates){
      if (TRUE %in% startsWith(dangers, temp.covariate)){
        danger.covariates <- c(danger.covariates, temp.covariate)
      }
    }
    if (length(danger.covariates) > 0){
      covariates <- covariates[!covariates %in% danger.covariates]
      covariates <- c(covariates, as.character(sapply(danger.covariates, function(x) paste("strata(",x,")",sep=""))))
    } else{
      summary <- summary(temp.model)
      zph <- cox.zph(temp.model)
    }
  }
  
  summary <- data.frame(summary$coefficients, check.names = F)
  validation <- zph
  validationTable <- validation$table
  figureCPH <- survfit(Surv(RFS,Recur) ~ Exp, data = InfoTable.result.uni,conf.type="log-log")
  p.value <- summary["ExpUpper",]$`Pr(>|z|)`
  if (is.na(p.value)){next}
  write.table(summaryTable.uni,file=paste(outputD,"/RFS_uni/summary/",geneSymbol,"_sig.txt",sep=""),row.names=T,col.names=T)
  write.table(validationTable.uni,file=paste(outputD,"/RFS_uni/validation/",geneSymbol,"_sig.txt",sep=""),row.names=T,col.names=T)
  
  label = paste("P value","\n","Univariate   : ", round(p.value.uni,5),"\n","Multivariate : ", round(p.value,5),sep="")
  #plt <- ggsurvplot(figureCPH, xlab="Time(months)",risk.table = F,conf.int=F,font.y=c(25,"bold","black"), font.x=c(25,"bold","black"),font.tickslab = c(20,"plain","black"), legend.title = "Exp", legend.labs = c("Low","High"),palette=c("red","blue"),font.legend=c(15,"plain","black")) 
  #plt$plot <- plt$plot + ggplot2::geom_text(x=0,y=0.2,label= label,size=7,hjust=0) #+ scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  if (p.value <= 0.05){
    #ggsave(paste(outputD,"/RFS/sig/",geneSymbol,".pdf",sep=""),plt$plot,dpi=300,width=6,height=6)
    #ggsave(paste(outputD,"/RFS/sig/",geneSymbol,".png",sep=""),plt$plot,dpi=300,width=6,height=6)
    write.table(summary,file=paste(outputD,"/RFS//summary/",geneSymbol,"_sig.txt",sep=""),row.names=T,col.names=T)
    write.table(validationTable,file=paste(outputD,"/RFS//validation/",geneSymbol,"_sig.txt",sep=""),row.names=T,col.names=T)
  } else {
    #ggsave(paste(outputD,"/RFS//plot/",geneSymbol,".pdf",sep=""),plt$plot,dpi=300,width=6,height=6)
    #ggsave(paste(outputD,"/RFS//plot/",geneSymbol,".png",sep=""),plt$plot,dpi=300,width=6,height=6)
    write.table(summary,file=paste(outputD,"/RFS//summary/",geneSymbol,".txt",sep=""),row.names=T,col.names=T)
    write.table(validationTable,file=paste(outputD,"/RFS//validation/",geneSymbol,".txt",sep=""),row.names=T,col.names=T)
  }
}
warnings()
