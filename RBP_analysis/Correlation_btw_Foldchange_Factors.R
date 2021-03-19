
type <- "Tsinghua"
if (type == "Catholic"){
  data.gene <- read.table("Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Qnorm_RPKMFiltered.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.isoform <- read.table("Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Isoform_Qnorm.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.fold <- read.table("Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_FP.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.miRNA <- read.table("New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt",header=T,row.names=1,check.names=F,quote="",sep="\t")
} else if (type == "TCGA"){
  data.gene <- read.table("Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Qnorm_RPKMFiltered_Intersect.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.isoform <- read.table("Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Isoform_Qnorm.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.fold <- read.table("Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_FP.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.miRNA <- read.table("New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt",header=T,row.names=1,check.names=F,quote="",sep="\t")
} else if (type == "Tsinghua"){
  data.gene <- read.table("Project/LiverCancer/Tsinghua/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_All_Qnorm_RPKMFiltered_Intersect.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.isoform <- read.table("Project/LiverCancer/Tsinghua/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_All_Isoform_noPVTT_Qnorm.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.fold <- read.table("Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDistAll_OffsetPlusMinus1_Frequency_FP.txt",header=T,row.names=1,check.names=F,quote="",sep='\t')
  data.miRNA <- read.table("New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt",header=T,row.names=1,check.names=F,quote="",sep="\t")
  colnames(data.gene)
  colnames(data.isoform)
  ncol(data.fold)
  ncol(data.miRNA)
  
  #data.gene <- cbind(data.gene[,1:20], data.gene[41:60], data.gene[21:40])
  #data.isoform <- cbind(data.isoform[,1:20], data.isoform[41:60], data.isoform[21:40])
  colnames(data.fold) <- colnames(data.gene)
  colnames(data.miRNA) <- colnames(data.gene)
}
data.fold <- 2 ** data.fold

#colnames(cbind(data.isoform[,1:20],data.isoform[,41:60]))
colnames(data.fold) == colnames(data.miRNA)
colnames(data.gene) == colnames(data.isoform)
outputD <- "New/LiverCancer/Catholic/RBP_analysis/"

data <- rbind(data.fold[rownames(data.fold) == 'hsa-mir-21',],data.miRNA[rownames(data.miRNA) == 'hsa-miR-21-5p',], data.miRNA[rownames(data.miRNA) == 'hsa-miR-21-5p_GCUUAUC_1',] + data.miRNA[rownames(data.miRNA) == 'hsa-miR-21-5p_UAGCUUA_-1',],data.gene[row.names(data.gene) == "U2AF2",],data.isoform[rownames(data.isoform) == 'ENST00000592790.1',],data.gene[row.names(data.gene) == "HNRNPC",])
data.t <- t(data)
colnames(data.t) <- c("Fold","Exp_mature","Exp_isomiR","U2AF2","VMP1","HNRNPC")
data.t <- data.frame(data.t,check.names=F)

for (RBP in c("HNRNPC","U2AF2","both")){
  if (RBP == "both"){
    RBP = c("HNRNPC","U2AF2")
  }
  for (tempType in c("Fold","Exp")){
    univariate <- data.frame()
    #motif <- "BasalUG"
    if (tempType == "Fold"){
      y <- "Fold"
      motifs <- c("Exp_isomiR","Exp_mature",RBP,"VMP1")
    } else {
      y <- "Exp_isomiR"
      motifs <- c("Exp_mature","Fold",RBP,"VMP1")
    }
    for (motif in motifs){
      print(motif)
      temp.formula <- as.formula(paste(y," ~ ",motif,"",sep=""))
      tempLogit <- lm(temp.formula,data=data.t)
      a <- summary(tempLogit)
      temp.frame <- data.frame(a$coefficients,check.names=F)
      temp.frame$Name <- motif
      univariate <- rbind(univariate,temp.frame[2:nrow(temp.frame),])
    }
    univariate
    univariate[univariate$`Pr(>|z|)` > 0.05,]
    write.table(x=univariate,file=paste(outputD,"LM_miR21_",tempType,"_",paste(RBP,collapse="_"),"_univariateGLM_",type,"_notLog.txt",sep=""),quote = F,sep = '\t',row.names = F,col.names = T)
    
    multivariate <- data.frame()
    
    temp.formula <- as.formula(paste(y," ~ ",paste(motifs,collapse="+"),sep=""))
    tempLogit <- lm(temp.formula,data=data.t)
    a <- summary(tempLogit)
    b <- data.frame(a$coefficients,check.names=F)[-1,]
    b$Name <- motifs
    multivariate <- rbind(multivariate,b)
    
    multivariate
    write.table(x=multivariate,file=paste(outputD,"LM_miR21_",tempType,"_",paste(RBP,collapse="_"),"_multivariateGLM_",type,"_notLog.txt",sep=""),quote = F,sep = '\t',row.names = F,col.names = T)
  }
}