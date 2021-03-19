def get_expTable(expF, source):
    expDic = dict(); 
    normalMinimum = []; cancerMinimum = []
    f = open(expF); lines = f.readlines(); f.close()
    header = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        line = line.strip().split('\t')
        name = line[0]
        exps = map(float,line[1:])
        expDic[name] = exps
        if source == 'Catholic':
            normalMinimum += filter(lambda x: x != 0, exps[:15])
            cancerMinimum += filter(lambda x: x != 0, exps[15:])
        elif source == 'TCGA':
            normalMinimum += filter(lambda x: x != 0, exps[:50])
            cancerMinimum += filter(lambda x: x != 0, exps[50:])
        elif source == 'Tsinghua':            
            normalMinimum += filter(lambda x: x != 0, exps[:20])
            cancerMinimum += filter(lambda x: x != 0, exps[20:])
    return expDic, header, min(normalMinimum), min(cancerMinimum)

def calc_DEG(expDic, header, normalMinimum, cancerMinimum):
    pGreaterDic = dict(); pLessDic = dict(); foldchangeDic = dict()
    for name in expDic.keys():
        normalExp = []
        cancerExp = []
        if source == 'Catholic':
            normalExp += expDic[name][:15]
            cancerExp += expDic[name][15:]
        elif source == 'TCGA':
            normalExp += expDic[name][:50]
            cancerExp += expDic[name][50:]
        elif source == 'Tsinghua':
            normalExp += expDic[name][:20]
            cancerExp += expDic[name][20:]
        
        normalMedian = np.median(normalExp)
        if normalMedian == 0:
            normalMedian = normalMinimum
        cancerMedian = np.median(cancerExp)
        if cancerMedian == 0:
            cancerMedian = cancerMinimum
        foldchangeDic[name] = cancerMedian / normalMedian

        sGreater, pGreater = stats.mannwhitneyu(normalExp,cancerExp,alternative="greater")
        sLess, pLess = stats.mannwhitneyu(normalExp,cancerExp,alternative="less")

        pGreaterDic[name] = pGreater
        pLessDic[name] = pLess

    fdrGreaterDic = adjust_pvalues(pGreaterDic)
    fdrLessDic = adjust_pvalues(pLessDic)
    return fdrGreaterDic, fdrLessDic, pGreaterDic, pLessDic, foldchangeDic

def adjust_pvalues(pvalueDic, method="BH"):
    stats = importr('stats')
    adjustDic = dict()
    temp_keys = []; temp_values = []
    for key in pvalueDic.keys():
        temp_keys.append(key)
        temp_values.append(pvalueDic[key])

    adjust_values = stats.p_adjust(FloatVector(temp_values),method=method)

    for i in xrange(len(temp_keys)):
        adjustDic[temp_keys[i]] = adjust_values[i]

    return adjustDic

def write_DEG(expDic,fdrDic,header, outputF, pDic, foldchangeDic):
    writer = open(outputF,'w')
    writer.write('Gene' + '\t' + '\t'.join(header) + '\n')
    writer2 = open('./testtest.txt','w')
    for name in expDic.keys():
        fdr = fdrDic[name]
        foldchange = foldchangeDic[name]
        if fdr <= 0.05 and abs(math.log2(foldchange)) >= 2:
            writer.write(name + '\t' + '\t'.join(map(str,expDic[name]) ) + '\n')
            writer2.write(name + '\t' + str(pDic[name]) + '\n')
    writer.close()
    writer2.close()

def write_DEG2(expDic, fdrGreaterDic, fdrLessDic, header, inputF, pGreaterDic, pLessDic, foldchangeDic):
    writer_Greater = open(inputF.split('.txt')[0] + '_GreaterDEG.txt','w')
    writer_Less = open(inputF.split('.txt')[0] + '_LessDEG.txt','w')
    writer_Total = open(inputF.split('.txt')[0] + '_TotalDEG.txt','w')
    writer_Greater.write('Gene' + '\t' + '\t'.join(header) + '\n')
    writer_Less.write('Gene' + '\t' + '\t'.join(header) + '\n')
    #writer_Total.write('Gene' + '\t' + '\t'.join(header) + '\n')
    writer_Total.write('Gene' + '\t' + 'Foldchange' + '\t' + 'FDR' + '\n')
    for name in expDic.keys():
        fdrGreater = fdrGreaterDic[name]
        beWritten  = name + '\t' + '\t'.join(map(str,expDic[name])) + '\n'
        foldchange = foldchangeDic[name]
        #if abs(math.log(foldchange,2)) < math.log(1.5,2): continue
        #is it possible that miRNAs can regulate expression of target genes more than 1.5 fold ?
        #the miRNAs even are neither overexpressed nor knock-out !
        #it should be thought more.
        if abs(math.log(foldchange,2)) >= math.log(1.5,2) and fdrGreater <= 0.05:
            writer_Greater.write(beWritten)
            #writer_Total.write(beWritten)

        fdrLess = fdrLessDic[name]
        if abs(math.log(foldchange,2)) >= math.log(1.5,2) and fdrLess <= 0.05:
            writer_Less.write(beWritten)
            #writer_Total.write(beWritten)
        writer_Total.write(name + '\t' + str(math.log(foldchange,2)) + '\t' + str(min(fdrGreater,fdrLess)) + '\n')
    writer_Greater.close()
    writer_Less.close()
    writer_Total.close()


if __name__=='__main__':
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import FloatVector
    from scipy import stats
    import numpy as np
    import math
    import sys
    if len(sys.argv) != 2:
        print "Usage : FiltereDEG.py Source"
        exit()
    #inputD = sys.argv[1]
    source = sys.argv[1]

    if source == 'Catholic':
        inputF = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
    elif source == 'TCGA':
        inputF = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
    elif source == 'Tsinghua':
        inputF = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt'

    expDic, header, normalMinimum, cancerMinimum = get_expTable(inputF,source)
    fdrGreaterDic, fdrLessDic, pGreaterDic, pLessDic, foldchangeDic = calc_DEG(expDic, header, normalMinimum, cancerMinimum)
    write_DEG2(expDic, fdrGreaterDic, fdrLessDic, header, inputF, pGreaterDic, pLessDic, foldchangeDic)
