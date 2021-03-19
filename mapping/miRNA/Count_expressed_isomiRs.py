def get_Annotation(annoF):
    annoDic = dict(); annoDic2 = dict()
    f = open(annoF); lines = f.readlines(); f.close()
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        lineType = line[2]
        name = line[8].split(';')[2].split('=')[1]

        if lineType == 'miRNA_primary_transcript':
            pre = name
            annoDic[pre] = []
        else:
            annoDic[pre].append(name)
            annoDic2[name] = pre
    return annoDic, annoDic2

def get_reAnnotation(reAnnoF, annoDic):
    reAnnoDic = dict()
    f = open(reAnnoF); lines = f.readlines(); f.close()
    for line in lines:
        line = line.strip().split(' ')
        if line[0] != "5'": continue
        
        pre = line[2]
        arm = line[3]
        offset = int(line[4])

        if len(annoDic[pre]) == 2:
            if arm == 'FP':
                mature = filter(lambda x: x.endswith('-5p'), annoDic[pre])[0]
            else:
                mature = filter(lambda x: x.endswith('-3p'), annoDic[pre])[0]
            reAnnoDic[mature] = offset
        elif len(annoDic[pre]) == 1:
            reAnnoDic[annoDic[pre][0]] = offset
        else:
            print pre, arm, offset
            exit()
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

def get_Exp(expF):
    expDic = dict()
    f = open(expF); lines = f.readlines(); f.close()
    totalsamples = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        line = line.strip().split('\t')
        mature = line[0]
        exps = map(float, line[1:])
        expDic[mature] = exps
    return expDic, totalsamples

def calc_expressedIsomiRs(expDic, armDic, annoDic2, reAnnoDic):
    expIsoDic = {'FP':[0,0,0,0,0,0,0,0,0,0,0], 'TP':[0,0,0,0,0,0,0,0,0,0,0]}
    eachDic = dict()
    for mi in expDic.keys():
        pre = annoDic2[mi.split('_')[0]]
        arm = armDic[mi.split('_')[0] + '_' + pre]
        
        if mi.count('_') > 0:
            offset = int(mi.split('_')[-1])
        else:
            offset = 0

        if reAnnoDic.has_key(mi.split('_')[0]):
            dominantOffset = reAnnoDic[mi.split('_')[0]]
            offset = offset - dominantOffset
        if offset < -5 or offset > 5: continue
        expIsoDic[arm][offset+5] += 1
        eachDic[mi + '\t' + arm] = offset
    return expIsoDic, eachDic
        
def write_result(expDic, expIsoDic, eachDic, outputD):
    writer = open(outputD + '/IsomiR_count.txt','w'); writer2 = open(outputD + '/IsomiR_position_information.txt','w')
    writer.write('Arm' + '\t' + '\t'.join(map(str, range(-5,6))) + '\n')
    writer2.write('miRNA' + '\t' + 'Arm' + '\t' + 'Position' + '\t' + 'MedianExp' + '\n')
    for arm in expIsoDic.keys():
        writer.write(arm + '\t' + '\t'.join(map(str, expIsoDic[arm])) + '\n')
    for miInfo in eachDic.keys():
        writer2.write(miInfo + '\t' + str(eachDic[miInfo]) + '\t' + str(np.median(expDic[miInfo.split('\t')[0]])) + '\n')
    writer.close(); writer2.close()

def main(matureF, preF, annoF, reAnnoF, expF):

    matureSeqDic = get_Sequence(matureF)
    preSeqDic = get_Sequence(preF)

    annoDic, annoDic2 = get_Annotation(annoF)
    reAnnoDic = get_reAnnotation(reAnnoF, annoDic)

    armDic = get_AnnoArm(annoDic, matureSeqDic, preSeqDic)

    expDic, totalsamples = get_Exp(expF)

    expIsoDic, eachDic = calc_expressedIsomiRs(expDic, armDic, annoDic2, reAnnoDic)
    
    write_result(expDic, expIsoDic, eachDic, '/'.join(expF.split('/')[:-1]))

if __name__=='__main__':
    import sys
    import numpy as np

    matureF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/miRNA.fa'
    preF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hairpin_hsa.fa'
    annoF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hsa.hg19.gff3'
    reAnnoF = '/home/seokju/Project/LiverCancer/src_rev1/miRDeep2/Isomir/Stratification/renew/Re-annotation_DrCS_DiCS.txt'
    if data == 'Catholic':
        expF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_Intersect.txt'
    elif data == 'TCGA':
        expF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_Intersect.txt'
    elif data == 'Tsinghua':
        expF = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_All_Qnorm_RPMFiltered_Intersect.txt'

    main(matureF, preF, annoF, reAnnoF, expF)
