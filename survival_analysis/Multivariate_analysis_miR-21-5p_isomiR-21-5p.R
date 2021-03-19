library(survival)
library(survminer)

inputFile = 'New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'

InfoTable <- read.table('Project/LiverCancer/src_rev1/SurvivalAnalysis/renew2/TCGA_New.txt',header=T,row.names=3,sep='\t',check.names=F)
InfoTable$OS <- as.numeric(as.character(InfoTable$OS))

covariates_ori <- c("neoplasm_histologic_grade","AFP","Age","Gender","VascularInvasion","Alcohol","HBV","HCV","T_stage")

InfoTable.uni.filtered <- InfoTable[,colnames(InfoTable) %in% c("OS","Survival")]
InfoTable.uni.filtered <- InfoTable.uni.filtered[InfoTable.uni.filtered$OS >= 0 & ! is.na(InfoTable.uni.filtered$OS) & InfoTable.uni.filtered$OS <= 1825,]
InfoTable.uni.filtered <- InfoTable.uni.filtered[complete.cases(InfoTable.uni.filtered),]

expTable <- read.table(inputFile,header=T,row.names=1,sep='\t',check.names = F)

geneSymbols <- c("hsa-miR-21-5p","hsa-miR-21-5p_GCUUAUC_1","hsa-miR-21-5p_UAGCUUA_-1","both")

InfoTable.result.uni <- InfoTable.uni.filtered
for (geneSymbol in geneSymbols){
  print(geneSymbol)
  covariates <- covariates_ori
  print(covariates)
    if (geneSymbol == "both"){
    temp.exp <- expTable["hsa-miR-21-5p_GCUUAUC_1",] + expTable["hsa-miR-21-5p_UAGCUUA_-1",]
    temp.exp.t.filtered <- t(temp.exp)
    rownames(temp.exp.t.filtered)
    } else{
      temp.exp <- expTable[geneSymbol,]
      temp.exp.t.filtered <- t(temp.exp)
      rownames(temp.exp.t.filtered)
    }
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
  InfoTable.result.uni <- cbind(InfoTable.result.uni, temp.exp.t.filtered.filtered.uni[,"Exp"][match(rownames(InfoTable.uni.filtered), rownames(temp.exp.t.filtered.filtered.uni))])
}
colnames(InfoTable.result.uni) <- c(colnames(InfoTable.uni.filtered),geneSymbols)
InfoTable.result.uni <- droplevels(InfoTable.result.uni)

survival.uni.result <- data.frame()
for (geneSymbol in geneSymbols){
  
  fitCPH.uni <- coxph(as.formula(paste("Surv(OS,Survival) ~ `",geneSymbol,"`",sep="")), data=InfoTable.result.uni)
  fitCPHsumm.uni <- summary(fitCPH.uni)
  print(geneSymbol)
  print(fitCPHsumm.uni)
  summaryTable.uni <- data.frame(fitCPHsumm.uni$coefficients, check.names = F)
  validation.uni <- cox.zph(fitCPH.uni)
  validationTable.uni <- validation.uni$table
  p.value.uni <- fitCPHsumm.uni$sctest["pvalue"]
  summaryTable.uni["P(LogRank)"] <- p.value.uni
  survival.uni.result <- rbind(survival.uni.result, summaryTable.uni)
}

model <- coxph(Surv(OS,Survival) ~ `hsa-miR-21-5p` + `hsa-miR-21-5p_GCUUAUC_1` + `hsa-miR-21-5p_UAGCUUA_-1`, data=InfoTable.result.uni)
model <- coxph(Surv(OS,Survival) ~ `hsa-miR-21-5p` + `both`, data=InfoTable.result.uni)
summary <- summary(model)
zph <- cox.zph(model)
summary <- data.frame(summary$coefficients, check.names = F)

# RFS


InfoTable <- read.table('Project/LiverCancer/src_rev1/SurvivalAnalysis/renew2/TCGA_New.txt',header=T,row.names=3,sep='\t',check.names=F)
InfoTable$RFS <- as.numeric(as.character(InfoTable$RFS))

covariates_ori <- c("neoplasm_histologic_grade","AFP","Age","Gender","VascularInvasion","Alcohol","HBV","HCV","T_stage")

InfoTable.uni.filtered <- InfoTable[,colnames(InfoTable) %in% c("RFS","Recur")]
InfoTable.uni.filtered <- InfoTable.uni.filtered[InfoTable.uni.filtered$RFS >= 0 & ! is.na(InfoTable.uni.filtered$RFS) & InfoTable.uni.filtered$RFS <= 1825,]
InfoTable.uni.filtered <- InfoTable.uni.filtered[complete.cases(InfoTable.uni.filtered),]

expTable <- read.table(inputFile,header=T,row.names=1,sep='\t',check.names = F)

geneSymbols <- c("hsa-miR-21-5p","hsa-miR-21-5p_GCUUAUC_1","hsa-miR-21-5p_UAGCUUA_-1","both")

InfoTable.result.uni <- InfoTable.uni.filtered
for (geneSymbol in geneSymbols){
  print(geneSymbol)
  covariates <- covariates_ori
  print(covariates)
  if (geneSymbol == "both"){
    temp.exp <- expTable["hsa-miR-21-5p_GCUUAUC_1",] + expTable["hsa-miR-21-5p_UAGCUUA_-1",]
    temp.exp.t.filtered <- t(temp.exp)
    rownames(temp.exp.t.filtered)
  } else{
    temp.exp <- expTable[geneSymbol,]
    temp.exp.t.filtered <- t(temp.exp)
    rownames(temp.exp.t.filtered)
  }
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
  InfoTable.result.uni <- cbind(InfoTable.result.uni, temp.exp.t.filtered.filtered.uni[,"Exp"][match(rownames(InfoTable.uni.filtered), rownames(temp.exp.t.filtered.filtered.uni))])
}
colnames(InfoTable.result.uni) <- c(colnames(InfoTable.uni.filtered),geneSymbols)
InfoTable.result.uni <- droplevels(InfoTable.result.uni)

survival.uni.result <- data.frame()
for (geneSymbol in geneSymbols){
  
  fitCPH.uni <- coxph(as.formula(paste("Surv(RFS,Recur) ~ `",geneSymbol,"`",sep="")), data=InfoTable.result.uni)
  fitCPHsumm.uni <- summary(fitCPH.uni)
  print(geneSymbol)
  print(fitCPHsumm.uni)
  summaryTable.uni <- data.frame(fitCPHsumm.uni$coefficients, check.names = F)
  validation.uni <- cox.zph(fitCPH.uni)
  validationTable.uni <- validation.uni$table
  p.value.uni <- fitCPHsumm.uni$sctest["pvalue"]
  summaryTable.uni["P(LogRank)"] <- p.value.uni
  survival.uni.result <- rbind(survival.uni.result, summaryTable.uni)
}

model <- coxph(Surv(RFS,Recur) ~ `hsa-miR-21-5p` + `hsa-miR-21-5p_GCUUAUC_1` + `hsa-miR-21-5p_UAGCUUA_-1`, data=InfoTable.result.uni)
model <- coxph(Surv(RFS,Recur) ~ `hsa-miR-21-5p` + `both`, data=InfoTable.result.uni)
summary <- summary(model)
zph <- cox.zph(model)
summary <- data.frame(summary$coefficients, check.names = F)
