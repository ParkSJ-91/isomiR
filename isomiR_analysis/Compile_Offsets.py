# output form
# miRNA Arm offset 5p (percentage) 3p (percentage, calculated from main 5p)

import sys

def compile_ReadCounts(inF):
    fpDic = dict(); tpDic = dict()
    # fpDic = {precursor: {arm: [percentage of offset fp from -5 to +5 of each sample], ...}, ...}
    # tpDic = {precursor: {arm: {fp offset: [percentage of offset tp from -5 to +5 of each sample], ...}, ...}, ...}
    f = open(inF); lines = f.readlines(); f.close()

    tempReads = []
    for line in lines[1:]:
        line = line.strip().split('\t')
        precursor = line[0]
        arm = line[1]
        sample = line[2]
        offset_FP = int(line[3])
        offset_TP = int(line[4])
        readcount = int(line[5])
    
        if len(tempReads) == 0:
            tempReads.append(
        


def main():
    
    source = sys.argv[1]

    if source == "Catholic":
        inF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/DoubleOffset_Distribution.txt'
    elif source == "TCGA":
        inF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/DoubleOffset_Distribution.txt'
    elif source == "Tsinghua":
        ifF = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/DoubleOffset_Distribution.txt'
    
    outD = '/home/seokju/New/LiverCancer/Catholic/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/Total/Frequency/'

    compile_ReadCounts(inF)


if __name__=='__main__':
    main()
