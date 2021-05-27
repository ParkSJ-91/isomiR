# isomiR

-Summary-

This project was launched to study the biogenesis and function of isomiR in cancer.

Over the course of the project, we confirmed that most of the isomiR is controlled by cis-factors, including sequence and structure. 
However, some miRNAs consistently exhibited abnormal isomiR expression in liver cancer samples.

The isomiR-21-5ps, which are 1nt offsets of 5' end of miR-21-5p, were the most representative dysregulated isomiRs and were also related to overal survival of HCC patients.
Those targets GHR irrespective of miR-21-5p, and those proved to be related to cancer progression _in vitro_ and _in vivo_.

Since no mutations on pri-mir-21, we suspected a trans-factors that regulates the expression of isomiR-21-5ps.
Through eCLIP-seq data, we found that two RNA binding proteins, hnRNPC and U2AF2, are related to the dysregulation of isomiR-21-5ps. 
Because hnRNPC showed a stronger effect on isomiR-21-5ps biogenesis, we were able to demonstrate the association of the hnRNPC-isomiR-21-5ps-GHR axis wich cancer progression _in vitro_.

-Code usage-
In the "mapping" folder,
  There are a series of codes for preprocessing, mapping, and calculating expression levels of miRNA-seq and RNA-seq data.
  We used miRDeep2 for miRNA-seq and BitSeq for RNA-seq.
  
In the "isomiR_analysis" folder,
  There are codes for calculating miRNA/isomiR expression leves and their ratio and finding associated cis-factors.
 
In the "targeting_analysis" folder,
  There are codes for priprotizing isomiRs and their targets using DEG analysis and TargetScan. 

In the "Survival_analysis" folder,
  There are codes for survival analysis.
  
In the "GHR_analysis" folder,
  There are codes for GO analysis of GHR correlated genes.
  
In the "RBP_analysis" folder,
  There are codes for finding trans-factors that regulates isomiR-21-5ps biogenesis.
