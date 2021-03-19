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
    smallest = []
    f = open(expF); lines = f.readlines(); f.close()
    totalsamples = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        line = line.strip().split('\t')
        mature = line[0]
        exps = map(float, line[1:])
        smallest += filter(lambda x: x != 0, exps)
        expDic[mature] = exps
    return expDic, totalsamples, min(smallest)

def get_DEMi(inF):  
    DEs = []
    f = open(inF); lines = f.readlines(); f.close()
    for line in lines[1:]:
        DEs.append(line.strip().split('\t')[0])
    return DEs

def get_ExpOrder(expDic, oriExpDic, smallest, reAnnoDic, source):
    pairExpDic = dict()
    orderDic = dict()
    corrDic = dict() # {Cano:[median,sd,[median,sd,corr,order(Major_1)], ...]}
    anovaDic = dict()
    ranksumUpDic = dict(); ranksumDownDic = dict()

    canoMiIDs = set(list(map(lambda x: x.split('_')[0], expDic.keys())))
    for canoMiID in canoMiIDs:
        if reAnnoDic.has_key(canoMiID):
            offset = reAnnoDic[canoMiID]
        else:
            offset = 0
        if offset == 0:
            canoMi = canoMiID
        else:
            canoMi = filter(lambda x: x.split('_')[0] == canoMiID and int(x.split('_')[-1]) == offset, filter(lambda y: y.count('_') > 0, expDic.keys()))
            if len(canoMi) != 1:
                print canoMiID, "has no canonical miRNA ?"
            else:
                canoMi = canoMi[0]
        if expDic.has_key(canoMi):
            #print(oriExpDic[canoMi])
            corrDic[canoMi] = [np.median(list(map(lambda x: math.log(x+smallest,10), oriExpDic[canoMi]))), np.std(list(map(lambda x: math.log(x+smallest,10),oriExpDic[canoMi])))]
        allMiExps = []
        tempDic = dict()
        for tempMi in filter(lambda x: x.split('_')[0] == canoMiID, expDic.keys()):
            #corrR, corrP = stats.pearsonr(expDic[tempMi], expDic[canoMi])
            if tempMi.count('_') > 0:
                offset = int(tempMi.split('_')[-1])
            else:
                offset = 0
            if reAnnoDic.has_key(tempMi.split('_')[0]):
                dominantOffset = reAnnoDic[tempMi.split('_')[0]]
                offset = offset - dominantOffset
            if offset == 0: continue # pass canonical miRNA
            allMiExps.append([tempMi, np.median(oriExpDic[tempMi])])

            if expDic.has_key(canoMi):
                # Calc isomiR ratio based on RPM
                ratio = [math.log((oriExpDic[tempMi][x]+0.1)/(oriExpDic[canoMi][x]+0.1),10) for x in range(len(oriExpDic[canoMi]))]
                #ratio = [(oriExpDic[tempMi][x]+0.1)/(oriExpDic[canoMi][x]+0.1) for x in range(len(oriExpDic[canoMi]))]
                pairExpDic[tempMi] = [expDic[canoMi], expDic[tempMi], ratio]
                if source == "Catholic":
                    ratio_N = ratio[:15]
                    ratio_C = ratio[15:]
                elif source == "TCGA":
                    ratio_N = ratio[:50]
                    ratio_C = ratio[50:]
                elif source == "Tsinghua":
                    ratio_N = ratio[:20]
                    ratio_C = ratio[20:]
                #if canoMi == 'hsa-miR-21-5p':
                #    print tempMi
                #    print ratio
                anovaDic[tempMi] = stats.f_oneway(ratio_N,ratio_C)[1]
                ranksumUpDic[tempMi] = stats.mannwhitneyu(ratio_N,ratio_C,alternative="greater")[1]
                ranksumDownDic[tempMi] = stats.mannwhitneyu(ratio_N,ratio_C,alternative="less")[1]

                #if canoMi == 'hsa-miR-21-5p' or canoMi == 'hsa-miR-10a-5p':
                #    print tempMi
                #    print np.median(list(map(lambda x: math.log(x+smallest,10), oriExpDic[tempMi])))
                #    print ratio
                #    print "Anova:",anovaDic[tempMi] 
                #    print "Ranksum:",ranksumUpDic[tempMi] ,ranksumDownDic[tempMi]
                #if expDic.has_key(canoMi):
                corrR, corrP = stats.pearsonr(oriExpDic[tempMi], oriExpDic[canoMi])
                dist = distance.euclidean(oriExpDic[tempMi], oriExpDic[canoMi])
                #tempDic[tempMi] = [np.median(list(map(lambda x: math.log(x,2), expDic[tempMi]))), np.std(list(map(lambda x: math.log(x,2), expDic[tempMi]))), dist, '_']
                tempDic[tempMi] = [tempMi, np.median(list(map(lambda x: math.log(x+smallest,10), oriExpDic[tempMi]))), np.std(list(map(lambda x: math.log(x+smallest,10), oriExpDic[tempMi]))), corrR, corrP, np.median(ratio), np.std(ratio), np.median(ratio_N), np.std(ratio_N), np.median(ratio_C), np.std(ratio_C)]
        allMiExps.sort(key=lambda row:row[1])

        for i in range(len(allMiExps)):
            isomiR, exp = allMiExps[::-1][i]
            orderDic[isomiR] = i
            if expDic.has_key(canoMi):
                corrDic[canoMi].append(tempDic[isomiR])
    #print corrDic['hsa-miR-21-5p']
    #print anovaDic['hsa-miR-21-5p_GCUUAUC_1']
    #print anovaDic['hsa-miR-21-5p_UAGCUUA_-1']
    #exit()
    anovaFdrDic = adjust_pvalues(anovaDic)
    ranksumUpFdrDic = adjust_pvalues(ranksumUpDic)
    ranksumDownFdrDic = adjust_pvalues(ranksumDownDic)

    ranksumFdrDic = dict()
    for id in ranksumUpFdrDic.keys():
        ranksumFdrDic[id] = min(ranksumUpFdrDic[id], ranksumDownFdrDic[id])

    return pairExpDic, orderDic, corrDic, anovaFdrDic, ranksumFdrDic

def adjust_pvalues(pvalueDic, method="BH"):
    tempStats = importr('stats')
    adjustDic = dict()
    temp_keys = []; temp_values = []
    for key in pvalueDic.keys():
        temp_keys.append(key)
        temp_values.append(pvalueDic[key])
    adjust_values = list(tempStats.p_adjust(FloatVector(temp_values),method=method))

    for i in xrange(len(temp_keys)):
        key = temp_keys.pop(0)
        FDR = adjust_values.pop(0)
        adjustDic[key] = FDR

    return adjustDic

def calc_expressedIsomiRs(expDic, armDic, annoDic2, reAnnoDic, orderDic):
    expIsoDic = {'All':{'FP':[0,0,0,0,0,0,0,0,0,0,0], 'TP':[0,0,0,0,0,0,0,0,0,0,0]},\
                 'Cano':{'FP':[0,0,0,0,0,0,0,0,0,0,0], 'TP':[0,0,0,0,0,0,0,0,0,0,0]},\
                 'Major_1':{'FP':[0,0,0,0,0,0,0,0,0,0,0], 'TP':[0,0,0,0,0,0,0,0,0,0,0]},\
                 'Major_2':{'FP':[0,0,0,0,0,0,0,0,0,0,0], 'TP':[0,0,0,0,0,0,0,0,0,0,0]},\
                 'Others':{'FP':[0,0,0,0,0,0,0,0,0,0,0], 'TP':[0,0,0,0,0,0,0,0,0,0,0]}}

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
        expIsoDic['All'][arm][offset+5] += 1
        if not orderDic.has_key(mi):
            expIsoDic['Cano'][arm][offset+5] += 1
            eachDic[mi+'\t'+arm] = [offset,'Cano']
        else:
            if orderDic[mi] == 0:
                expIsoDic['Major_1'][arm][offset+5] += 1
                eachDic[mi+'\t'+arm] = [offset,'Major_1']
            elif orderDic[mi] == 1:
                expIsoDic['Major_2'][arm][offset+5] += 1
                eachDic[mi+'\t'+arm] = [offset,'Major_2']
            else:
                expIsoDic['Others'][arm][offset+5] += 1
                eachDic[mi+'\t'+arm] = [offset,'Others']
        #eachDic[mi + '\t' + arm] = offset
    return expIsoDic, eachDic
        
def write_result(expDic, expIsoDic, eachDic, pairExpDic, orderDic, corrDic, anovaDic, ranksumDic, armDic, annoDic2, outputD, source):
    writer = open(outputD + '/IsomiR_count_v2.txt','w'); writer2 = open(outputD + '/IsomiR_position_information_v2.txt','w')
    writer.write('Arm' + '\t' + '\t'.join(map(str, range(-5,6))) + '\t' + 'Order' + '\n')
    writer2.write('miRNA' + '\t' + 'Arm' + '\t' + 'Position' + '\t' + 'MedianExp' + '\t' + 'Order' + '\n')
    for order in expIsoDic.keys():
        for arm in expIsoDic[order].keys():
            writer.write(arm + '\t' + '\t'.join(map(str, expIsoDic[order][arm])) + '\t' + order + '\n')
    for miInfo in eachDic.keys():
        writer2.write(miInfo + '\t' + str(eachDic[miInfo][0]) + '\t' + str(np.median(expDic[miInfo.split('\t')[0]])) + '\t' + eachDic[miInfo][1] + '\n')
    writer.close(); writer2.close()

    writer = open(outputD + '/Correlation_Mature_AllIsomiRs_AllSamples.txt','w')
    writer.write('IsomiR' + '\t' + 'Sample' + '\t' + 'Precursor' + '\t' + 'Arm' + '\t' + 'Offset' + '\t' + 'CanoExp' + '\t' + 'IsomiRExp' + '\t' + 'PreExp' + '\t' + 'Ratio' + '\t' + 'Source' + '\n')
    for isomiR in pairExpDic.keys():

        pre = annoDic2[isomiR.split('_')[0]]
        arm = armDic[isomiR.split('_')[0] + '_' + pre]
        offset = eachDic[isomiR+'\t'+arm][0]

        canoMi = filter(lambda x: x.split('_')[0] == isomiR.split('_')[0], corrDic.keys())
        if len(canoMi) != 1:
            print isomiR, "doesn't have correct canonical miRNA"
            exit()
        canoMi = canoMi[0]
        allIsomiRs = map(lambda x: x[0], corrDic[canoMi][2:])
        for i in xrange(len(pairExpDic[isomiR][0])):
            Sample = source + '_' + str(i+1)
            CanoExp = pairExpDic[isomiR][0][i]
            IsomiRExp = pairExpDic[isomiR][1][i]
            Ratio = pairExpDic[isomiR][2][i]
            AllIsomiRExp = sum(map(lambda x: expDic[x][i], allIsomiRs))

            writer.write('\t'.join(map(str,[isomiR, Sample,pre,arm,offset,CanoExp,IsomiRExp,CanoExp+AllIsomiRExp,Ratio,source])) + '\n')

    writer.close()

    writer = open(outputD + '/Correlation_Mature_AllIsomiRs.txt','w')
    writer.write('isomiR' + '\t' + 'Precursor' + '\t' + 'Arm' + '\t' + 'Median_Cano' + '\t' + 'SD_Cano' + '\t' + 'CorrR_IsomiR' + '\t' + 'CorrP_IsomiR' + '\t' + 'Median_IsomiR' + '\t' + 'SD_IsomiR' + '\t' + 'Median_Ratio' + '\t' + 'SD_Ratio' + '\t' + 'Median_Ratio_Normal' + '\t' + 'SD_Ratio_Normal' + '\t' + 'Median_Ratio_Cancer' + '\t' + 'SD_Ratio_Cancer' + '\t' + 'Anova' + '\t' + 'Ranksum' + '\n')
    for mi in corrDic.keys():

        pre = annoDic2[mi.split('_')[0]]
        arm = armDic[mi.split('_')[0] + '_' + pre]

        median_cano = corrDic[mi][0]
        sd_cano = corrDic[mi][1]

        for isomiR_info in corrDic[mi][2:]:
            isomiR, median, sd, corrR, corrP, median_ratio, sd_ratio, median_ratio_normal, sd_ratio_normal, median_ratio_cancer, sd_ratio_cancer = isomiR_info
            writer.write('\t'.join(map(str,[isomiR,pre,arm,median_cano,sd_cano,corrR,corrP,median,sd,median_ratio,sd_ratio,median_ratio_normal,sd_ratio_normal,median_ratio_cancer,sd_ratio_cancer,anovaDic[isomiR],ranksumDic[isomiR]])) + '\n')

    writer.close()

    writer = open(outputD + '/Correlation_Mature_MajorIsomiR.txt','w')
    writer.write('Precursor' + '\t' + 'Arm' + '\t' + 'Median_Cano' + '\t' + 'SD_Cano' + '\t' + 'Median_Major_1' + '\t' + 'SD_Major_1' + '\t' + 'CorrR_Major_1' + '\t' + 'CorrP_Major_1' + '\t' + 'Median_Ratio_Major_1' + '\t' + 'SD_Ratio_Major_1' + '\t' + 'Median_Ratio_Normal_Major_1' + '\t' + 'SD_Ratio_Normal_Major_1' + '\t' + 'Median_Ratio_Cancer_Major_1' + '\t' + 'SD_Ratio_Cancer_Major_1' + '\t' + 'Anova_Major_1' + '\t' + 'Ransum_Major_1' + '\t' + 'Median_Major_2' + '\t' + 'SD_Major_2' + '\t' + 'CorrR_Major_2' + '\t' + 'CorrP_Major_2' + '\t' + 'Median_Ratio_Major_2' + '\t' + 'SD_Ratio_Major_2' + '\t' + 'Median_Ratio_Normal_Major_2' + '\t' + 'SD_Ratio_Normal_Major_2' + '\t' + 'Median_Ratio_Cancer_Major_2' + '\t' + 'SD_Ratio_Cancer_Major_2' + '\t' + 'Anova_Major_2' + '\t' + 'Ranksum_Major_2' + '\n')
    for mi in corrDic.keys():

        pre = annoDic2[mi.split('_')[0]]
        arm = armDic[mi.split('_')[0] + '_' + pre]

        median_cano = corrDic[mi][0]
        sd_cano = corrDic[mi][1]

        if len(corrDic[mi]) > 2:
            isomiR_major_1, median_major_1, sd_major_1, corrR_major_1, corrP_major_1, median_ratio_major_1, sd_ratio_major_1, median_ratio_normal_major_1, sd_ratio_normal_major_1, median_ratio_cancer_major_1, sd_ratio_cancer_major_1 = corrDic[mi][2]
            anovaP_major_1 = anovaDic[isomiR_major_1]
            ranksumP_major_1 = ranksumDic[isomiR_major_1]
        else:
            median_major_1 = 'NA'; sd_major_1 = 'NA'; corrR_major_1 = 'NA'; corrP_major_1 = 'NA'; median_ratio_major_1 = 'NA'; sd_ratio_major_1 = 'NA'; median_ratio_normal_major_1 = 'NA'; sd_ratio_normal_major_1 = 'NA'; median_ratio_cancer_major_1 = 'NA'; sd_ratio_cancer_major_1 = 'NA'
            anovaP_major_1 = 'NA'
            ranksumP_major_1 = 'NA'
        if len(corrDic[mi]) > 3:
            isomiR_major_2, median_major_2, sd_major_2, corrR_major_2, corrP_major_2, median_ratio_major_2, sd_ratio_major_2, median_ratio_normal_major_2, sd_ratio_normal_major_2, median_ratio_cancer_major_2, sd_ratio_cancer_major_2 = corrDic[mi][3]
            anovaP_major_2 = anovaDic[isomiR_major_2]
            ranksumP_major_2 = ranksumDic[isomiR_major_2]
        else:
            median_major_2 = 'NA'; sd_major_2 = 'NA'; corrR_major_2 = 'NA'; corrP_major_2 = 'NA'; corrP_major_2 = 'NA'; median_ratio_major_2 = 'NA'; sd_ratio_major_2 = 'NA'; median_ratio_normal_major_2 = 'NA'; sd_ratio_normal_major_2 = 'NA'; median_ratio_cancer_major_2 = 'NA'; sd_ratio_cancer_major_2 = 'NA'
            anovaP_major_2 = 'NA'
            ranksumP_major_2 = 'NA'
            
        writer.write('\t'.join(map(str,[pre, arm, median_cano, sd_cano, median_major_1, sd_major_1, corrR_major_1, corrP_major_1, median_ratio_major_1, sd_ratio_major_1, median_ratio_normal_major_1, sd_ratio_normal_major_1, median_ratio_cancer_major_1, sd_ratio_cancer_major_1, anovaP_major_1, ranksumP_major_1, median_major_2, sd_major_2, corrR_major_2, corrP_major_2, median_ratio_major_2, sd_ratio_major_2, median_ratio_normal_major_2, sd_ratio_normal_major_2, median_ratio_cancer_major_2, sd_ratio_cancer_major_2, anovaP_major_2, ranksumP_major_2])) + '\n')
    writer.close()

    writer = open(outputD + '/DE_isomiR_ratio.txt','w')
    writer.write('isomiR' + '\t' + 'order' + '\t' + 'Anova' + '\t' + 'Ranksum' + '\n')
    for mi in anovaDic.keys():
        writer.write('\t'.join(map(str,[mi, orderDic[mi], anovaDic[mi], ranksumDic[mi]])) + '\n')
    writer.close()

def main(matureF, preF, annoF, reAnnoF, expF, oriExpF, DEMiF, data):

    matureSeqDic = get_Sequence(matureF)
    preSeqDic = get_Sequence(preF)

    annoDic, annoDic2 = get_Annotation(annoF)
    reAnnoDic = get_reAnnotation(reAnnoF, annoDic)

    DEMis = get_DEMi(DEMiF)
    armDic = get_AnnoArm(annoDic, matureSeqDic, preSeqDic)

    expDic, totalsamples, _ = get_Exp(expF)
    oriExpDic, _, smallest = get_Exp(oriExpF)

    pairExpDic, orderDic, corrDic, anovaDic, ranksumDic = get_ExpOrder(expDic, oriExpDic, smallest, reAnnoDic, data)
    expIsoDic, eachDic = calc_expressedIsomiRs(expDic, armDic, annoDic2, reAnnoDic, orderDic)
    
    for DEMi in DEMis:
        #if DEMi.count('_') == 0: continue
        if anovaDic.has_key(DEMi):
            print DEMi, orderDic[DEMi], anovaDic[DEMi], ranksumDic[DEMi]
    #exit()
    write_result(expDic, expIsoDic, eachDic, pairExpDic, orderDic, corrDic, anovaDic, ranksumDic, armDic, annoDic2, '/'.join(expF.split('/')[:-1]), data)

if __name__=='__main__':
    import sys
    import math
    import numpy as np
    from scipy import stats
    from scipy.spatial import distance
    from rpy2.robjects.vectors import FloatVector
    from rpy2.robjects.packages import importr

    data = sys.argv[1]
    matureF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/miRNA.fa'
    preF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hairpin_hsa.fa'
    annoF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hsa.hg19.gff3'
    reAnnoF = '/home/seokju/Project/LiverCancer/src_rev1/miRDeep2/Isomir/Stratification/renew/Re-annotation_DrCS_DiCS.txt'
    DEMiF = '/home/seokju/Project/LiverCancer/Catholic/3.TargetAnalysis/9.TargetAnalysis_noPVTT/TotalCorrelation_AllTarget_GeneLevel_include6mer/DEMi_Cancer.txt'
    if data == 'Catholic':
        expF = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
        oriExpF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_withoutSNV3/miRDeep2_Expression_neoplasm.txt'
    elif data == 'TCGA':
        expF = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
        oriExpF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_withoutSNV3/miRDeep2_Expression_neoplasm.txt'
    elif data == 'Tsinghua':
        expF = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt'
        oriExpF = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_All.txt'
    main(matureF, preF, annoF, reAnnoF, expF, oriExpF, DEMiF, data)
