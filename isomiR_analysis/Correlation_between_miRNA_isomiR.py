import sys
import numpy as np 
from scipy import stats

def get_Exp(expF):
    expDic = dict()
    f = open(expF); lines = f.readlines(); f.close()
    totalsamples = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        line = line.strip().split('\t')
        mi = line[0]
        if mi.count('_') != 0:
            mi = mi.split('_')[0] + '_' + mi.split('_')[2]
        exps = map(float, line[1:])
        expDic[mi] = exps
    return expDic, totalsamples

def get_Reannotation(reAnnoF):
    reAnnoDic = dict()
    f = open(reAnnoF); lines = f.readlines(); f.close()
    for line in lines[1:]:
        line = line.strip().split('\t')
        mi = line[0]
        arm = line[1]
        offset = int(line[7])
        if not reAnnoDic.has_key(mi): reAnnoDic[mi] = dict()
        reAnnoDic[mi][arm] = offset
    return reAnnoDic

def get_Sequence(seqF):
    seqDic = dict()
    f = open(seqF); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        lines = chunk.split('\n')
        priName = lines[0].split(' ')[0]
        seq = ''.join(lines[1:])
        seqDic[priName] = seq 
    return seqDic

def get_AnnoArm(pairDic, matureSeqDic, hairpinDic):
    armDic = dict()
    for pre in pairDic.keys():
        preSeq = hairpinDic[pre]
        for mature in pairDic[pre]:
            matureSeq = matureSeqDic[mature]
            f = 0; t = 0 
            if mature.endswith('-5p'):
                f = 1 
            elif mature.endswith('-3p'):
                t = 1 
            elif preSeq.index(matureSeq) < len(preSeq) - (len(matureSeq) + preSeq.index(matureSeq)):
                f = 1 
            else:
                t = 1 
            if f == 1 and t == 0:
                armDic[mature + '_' + pre] = 'FP' 
            elif f == 0 and t == 1:
    
                armDic[mature + '_' + pre] = 'TP' 
            else:
                continue
    return armDic

def get_Pair(pairF):
    pairDic = dict(); pairDic2 = dict()
    f = open(pairF); lines = f.readlines(); f.close()
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        name = line[8].split(';')[2].split('=')[1]
        if line[2] == 'miRNA_primary_transcript':
            pri = name
            pairDic2[pri] = []
        else:
            pairDic[name] = pri
            pairDic2[pri].append(name)
    return pairDic, pairDic2

def get_MatchedSet(expDic, pairDic, reAnnoDic, armDic):
    matchDic = dict()

    for mi in expDic.keys():
        if mi.count('_') == 0: continue
        #pre = pairDic[mi]

        mi_ori = mi.split('_')[0]
        pre = pairDic[mi_ori]
        offset = int(mi.split('_')[1])
        arm = armDic[mi_ori+'_'+pre]

        if reAnnoDic.has_key(mi_ori):
            offset_balancer = reAnnoDic[mi_ori][arm]
        else:
            offset_balancer = 0

        if offset - offset_balancer == 1 or offset - offset_balancer == -1:

            if offset_balancer == 0:
                mi_new = mi_ori
            else:
                mi_new = mi_ori + '_' + str(offset_balancer * -1)
            if not matchDic.has_key(mi_new): matchDic[mi_new] = []
            matchDic[mi_new].append(mi)
    return matchDic

def write_result(matchDic, expDic, outPrefix, pairDic, armDic, totalsamples):
    #writer_mature = open(outPrefix + '_MatchedMature.txt','w')
    #writer_isomiR = open(outPrefix + '_MatchedIsomiR.txt','w')
    #writer_corr = open(outPrefix + '_MatchedCorrelation.txt','w')
    #writer_mature.write('Name' + '\t' + '\t'.join(totalsamples) + '\n')
    #writer_isomiR.write('Name' + '\t' + '\t'.join(totalsamples) + '\n')
    #writer_corr.write('miRNA' + '\t' + 'Precursor' + '\t' + 'Arm' + '\t' + 'PearsonR' + '\t' + 'PearsonP' + '\n')

    for mi_new in matchDic.keys():
        if not expDic.has_key(mi_new):
            print mi_new, matchDic[mi_new]
            continue
            #exit()
        if len(matchDic[mi_new]) == 0: continue
        #writer_mature.write(mi_new + '\t' + '\t'.join(map(str,expDic[mi_new])) + '\n')
        print mi_new, matchDic[mi_new]
        if len(matchDic[mi_new]) == 1:
            isomiR_exp = expDic[matchDic[mi_new][0]]
        else:
            isomiR_exp = np.asarray(expDic[matchDic[mi_new][0]]) + np.asarray(expDic[matchDic[mi_new][1]])
        #writer_isomiR.write(mi_new + '\t' + '\t'.join(map(str,isomiR_exp)) + '\n')

        pre = pairDic[mi_new.split('_')[0]]
        arm = armDic[mi_new.split('_')[0]+'_'+pre]

        pearsonR, pearsonP = stats.pearsonr(expDic[mi_new], isomiR_exp)
        #writer_corr.write(mi_new + '\t' + pre + '\t' + arm + '\t' + str(pearsonR) + '\t' + str(pearsonP) + '\n')

    #writer_mature.close()
    #writer_isomiR.close()
    #writer_corr.close()

def main():

    dataset = sys.argv[1] # Catholic / TCGA / Tsinghua
    if dataset == 'Catholic':
        expF = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
    elif dataset == 'TCGA':
        expF = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
    elif dataset == 'Tsinghua':
        expF = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt'
    outPrefix = expF.split('.txt')[0]

    reAnnoF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_v2_cutoff10_DD/Total/Frequency/AlleleFrequencyDist_OffsetPlusMiNum1_Frequency_IntersectSet.txt'
    matureF = '/home/seokju/New/LiverCancer/data/miRNA.fa'
    preF = '/home/seokju/New/LiverCancer/data/hairpin_hsa.fa'
    pairF = '/home/seokju/New/LiverCancer/data/hsa.hg19.gff3'

    expDic, totalsamples = get_Exp(expF)
    matureSeqDic = get_Sequence(matureF)
    preSeqDic = get_Sequence(preF)
    pairDic, pairDic2 = get_Pair(pairF)
    reAnnoDic = get_Reannotation(reAnnoF)

    armDic = get_AnnoArm(pairDic2, matureSeqDic, preSeqDic)

    matchDic = get_MatchedSet(expDic, pairDic, reAnnoDic, armDic)
    write_result(matchDic, expDic, outPrefix, pairDic, armDic, totalsamples)


if __name__=="__main__":

    main()
