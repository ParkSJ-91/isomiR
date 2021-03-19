def get_hairpin(hairpinF):
    hairpinDic = dict()
    f = open(hairpinF); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        lines = chunk.split('\n')
        priName = lines[0].split(' ')[0]
        seq = ''.join(lines[1:])
        hairpinDic[priName] = seq
    return hairpinDic

def get_pairDic(pairF):
    pairDic = dict(); coordDic = dict()
    f = open(pairF); lines = f.readlines(); f.close()
    pre = '' 
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        chr = line[0]
        mirType = line[2]
        start = line[3]
        end = line[4]
        strand = line[6]
        info = line[8].split(';')
        id = info[2].split('=')[1]
        if mirType == 'miRNA_primary_transcript':
            pre = id
            coordDic[pre] = [chr,start,end,strand]
        else:
            if not pairDic.has_key(pre): pairDic[pre] = []
            pairDic[pre].append(id)
    pairDic['hsa-mir-1302-2'] = ['hsa-miR-1302']
    pairDic['hsa-mir-939'] = ['hsa-miR-939-5p','hsa-miR-939-3p']
    return pairDic, coordDic

def get_DroshaDependentSet(droshaDependencyF):
    droshaDependentSet = []; dicerDependentSet = []
    droshaFoldchangeDic_rep1 = dict(); droshaNormalizedDic_rep1 = dict()
    droshaFoldchangeDic_rep2 = dict(); droshaNormalizedDic_rep2 = dict()
    dicerFoldchangeDic_rep1 = dict(); dicerNormalizedDic_rep1 = dict()
    dicerFoldchangeDic_rep2 = dict(); dicerNormalizedDic_rep2 = dict()
    f = open(droshaDependencyF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        name = line[0]
        Drosha_rep1 = float(line[3])
        Drosha_rep2 = float(line[4])
        Dicer_rep1 = float(line[1])
        Dicer_rep2 = float(line[2])
        WT_rep1 = float(line[8])
        WT_rep2 = float(line[9])
        WT_Dicer = float(line[7])

        droshaFoldchangeDic_rep1[name] = Drosha_rep1 / WT_rep1
        droshaFoldchangeDic_rep2[name] = Drosha_rep2 / WT_rep2
        
        dicerFoldchangeDic_rep1[name] = Dicer_rep1 / WT_Dicer
        dicerFoldchangeDic_rep2[name] = Dicer_rep2 / WT_Dicer

    for name in droshaFoldchangeDic_rep1.keys():
        droshaNormalizedDic_rep1[name] = math.log(droshaFoldchangeDic_rep1[name] / droshaFoldchangeDic_rep1['hsa-miR-320a'],10)
        droshaNormalizedDic_rep2[name] = math.log(droshaFoldchangeDic_rep2[name] / droshaFoldchangeDic_rep2['hsa-miR-320a'],10)
        
        dicerNormalizedDic_rep1[name] = math.log(dicerFoldchangeDic_rep1[name] / dicerFoldchangeDic_rep1['hsa-miR-451a'],10)
        dicerNormalizedDic_rep2[name] = math.log(dicerFoldchangeDic_rep2[name] / dicerFoldchangeDic_rep2['hsa-miR-451a'],10)

    for name in droshaNormalizedDic_rep1.keys():
        if droshaNormalizedDic_rep1[name] < -2 or droshaNormalizedDic_rep2[name] < -2:
            droshaDependentSet.append(name)
        if dicerNormalizedDic_rep1[name] < -2 or dicerNormalizedDic_rep2[name] < -2:
            dicerDependentSet.append(name)
    return droshaDependentSet, dicerDependentSet

def get_Exp(expF):
    expDic = dict()
    f = open(expF); lines = f.readlines(); f.close()
    totalsamples = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        line = line.strip().split('\t')
        name = line[0]
        exps = map(float,line[1:])
        expDic[name] = exps
    return expDic, totalsamples

def get_AnnoArm(pairDic, matureSeqDic, hairpinDic, mirtrons, droshaDependentSet, dicerDependentSet):
    armDic = dict(); armDic2 = dict()
    for pre in pairDic.keys():
        if pre in mirtrons: continue
        preSeq = hairpinDic[pre]
        armDic2[pre] = dict()
        for mature in pairDic[pre]:
            if not mature in droshaDependentSet: continue
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
                armDic2[pre]['FP'] = mature
            elif f == 0 and t == 1:
                
                armDic[mature + '_' + pre] = 'TP'
                armDic2[pre]['TP'] = mature
            else:
                continue
    return armDic, armDic2

def get_InfoDic(InfoF, gradeType):
    if gradeType == 'Total':
        grades = ['Normal','G1','G2','G3','G4','HPC','Cancer']
    elif gradeType == 'Normal':
        grades = ['Normal']
    elif gradeType == 'Cancer':
        grades = ['G1','G2','G3','G4','HPC','Cancer']
    totalsamples = []
    InfoDic = dict(); InfoDic2 = dict()
    f = open(InfoF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        if len(line) > 6:
            grade = line[1]
            sampleN = line[2]
        else:
            grade = line[1]
            sampleN = line[3]
        if sampleN in ['Sample_25','Sample_68','Sample_85','Sample_86','Sample_87','Sample_88']: continue
        if not grade in grades: continue
        if not InfoDic.has_key(grade): InfoDic[grade] = []
        InfoDic[grade].append(sampleN)
        InfoDic2[sampleN] = grade
    for grade in grades:
        if not InfoDic.has_key(grade): continue
        samples = InfoDic[grade]
        samples.sort()
        totalsamples += samples
    return InfoDic, InfoDic2, totalsamples

def get_ReadCount(IsomiRD, totalsamples, armDic, infoDic):
    rcDic = dict(); medianRCDic = dict(); dOffsetDic = dict()

    for i in xrange(len(totalsamples)):
        sample = totalsamples[i]
        targetFile = filter(lambda x: x.startswith(sample) and x.endswith('_OffsetPerfect.txt'), os.listdir(IsomiRD))
        if len(targetFile) == 0:
            continue
        elif len(targetFile) == 1: pass
        else:
            print targetFile
            exit()
        f = open(IsomiRD + targetFile[0]); lines = f.readlines(); f.close()
        tempDic = dict()
        for line in lines:
            line = line.strip().split('\t')
            mature = line[0]
            pre = line[1]
            offset = line[3]
            IsomirAllele = int(line[6])
            MatureAllele = int(line[7])
            if not armDic.has_key(mature + '_' + pre):
                continue
            arm = armDic[mature +'_'+ pre]
            if not rcDic.has_key(pre): rcDic[pre] = dict()
            if not rcDic[pre].has_key(arm):
                rcDic[pre][arm] = []
                for n in xrange(11):
                    rcDic[pre][arm].append([0 for y in xrange(len(totalsamples))])

            rcDic[pre][arm][int(offset) + 5][i] += IsomirAllele# / float(len(totalsamples))
            if not tempDic.has_key(pre + '_' + arm):
                rcDic[pre][arm][5][i] += MatureAllele # / float(len(totalsamples))
                tempDic[pre+'_'+arm] = 'ok'

        # process double offset files
        dOffsetF = filter(lambda x: x.startswith(sample) and x.endswith("_DoubleOffset.txt"), os.listdir(IsomiRD))
        if len(dOffsetF) == 0:
            continue
        elif len(dOffsetF) == 1: pass
        else:
            print dOffsetF
            exit()
        f = open(IsomiRD + dOffsetF[0]); lines = f.readlines(); f.close()
        for line in lines:
            line = line.strip().split('\t')
            mature = line[0]
            fpOffset = int(line[1].split('_')[1])
            pre = line[2]
            if not armDic.has_key(mature+'_'+pre): continue
            arm = armDic[mature + '_' + pre]
            tempReadCounts = map(int, line[3:]) # tpOffset -5 ~ 5
            if not dOffsetDic.has_key(pre): dOffsetDic[pre] = dict()
            if not dOffsetDic[pre].has_key(arm):
                dOffsetDic[pre][arm] = []
                for n in xrange(11):
                    dOffsetDic[pre][arm].append([[0 for y in xrange(11)] for z in xrange(len(totalsamples))])
            for tpOffset in xrange(len(tempReadCounts)):
                dOffsetDic[pre][arm][fpOffset+5][i][tpOffset] += tempReadCounts[tpOffset]
            
    for pre in rcDic.keys():
        medianRCDic[pre] = dict()
        for arm in rcDic[pre].keys():
            medianRCDic[pre][arm] = []
            for offset in range(len(rcDic[pre][arm])):
                temp = []
                if len(rcDic[pre][arm][offset]) == 0:
                    temp += [0,0,0]
                else:
                    if len(rcDic[pre][arm][offset]) != len(totalsamples):
                        print pre, arm, offset
                        exit()
                    normal = []; cancer = []
                    for i in xrange(len(totalsamples)):
                        sample = totalsamples[i]
                        if infoDic[sample] == 'Normal':
                            normal.append(rcDic[pre][arm][offset][i])
                        else:
                            cancer.append(rcDic[pre][arm][offset][i])
                    temp += [np.mean(normal+cancer), np.mean(normal), np.mean(cancer)]
                medianRCDic[pre][arm].append(temp)
    return rcDic, medianRCDic, dOffsetDic

def adjust_offset(rcDic_cat, rcDic_tcga, rcDic_tsinghwa, medianRcDic_cat, medianRcDic_tcga, medianRcDic_tsinghwa, dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa):
    changedDic = dict()

    for pre in rcDic_cat.keys():
        if not rcDic_tcga.has_key(pre): continue
        if not rcDic_tsinghwa.has_key(pre): continue
        for arm in rcDic_cat[pre].keys():
            if not rcDic_tcga[pre].has_key(arm): continue
            if not rcDic_tsinghwa[pre].has_key(arm): continue
            zeroCat = 0; zeroTCGA = 0; zeroTsinghwa = 0
            maxIndices_cat = []
            for i in xrange(len(rcDic_cat[pre][arm][0])):
                exps = map(lambda x: x[i], rcDic_cat[pre][arm])
                if max(exps) == 0:
                    zeroCat += 1
                    continue
                maxIndices_cat += [index for index,value in enumerate(exps) if value == max(exps)]
            maxIndices_tcga = []
            for i in xrange(len(rcDic_tcga[pre][arm][0])):
                exps = map(lambda x: x[i], rcDic_tcga[pre][arm])
                if max(exps) == 0:
                    zeroTCGA += 1
                    continue
                maxIndices_tcga += [index for index,value in enumerate(exps) if value == max(exps)]
            maxIndices_tsinghwa = []
            for i in xrange(len(rcDic_tsinghwa[pre][arm][0])):
                exps = map(lambda x: x[i], rcDic_tsinghwa[pre][arm])
                if max(exps) == 0:
                    zeroTsinghwa += 1
                    continue
                maxIndices_tsinghwa += [index for index,value in enumerate(exps) if value == max(exps)]

            HitIndices_cat = []
            if len(rcDic_cat[pre][arm][0])-zeroCat != 0:
                for i in xrange(len(rcDic_cat[pre][arm])):
                    if maxIndices_cat.count(i) >= (len(rcDic_cat[pre][arm][0])-zeroCat)*0.9:
                        HitIndices_cat.append(i)

            HitIndices_tcga = []
            if len(rcDic_tcga[pre][arm][0])-zeroTCGA != 0:
                for i in xrange(len(rcDic_tcga[pre][arm])):
                    if maxIndices_tcga.count(i) >= (len(rcDic_tcga[pre][arm][0])-zeroTCGA)*0.9:
                        HitIndices_tcga.append(i)

            HitIndices_tsinghwa = []
            if len(rcDic_tsinghwa[pre][arm][0])-zeroTsinghwa != 0:
                for i in xrange(len(rcDic_tsinghwa[pre][arm])):
                    if maxIndices_tsinghwa.count(i) >= (len(rcDic_tsinghwa[pre][arm][0])-zeroTsinghwa)*0.9:
                        HitIndices_tsinghwa.append(i)

            Hit_Cat = 0; Hit_TCGA = 0; Hit_Tsinghwa = 0
            hits = []
            if len(HitIndices_cat) == 1 and not 5 in HitIndices_cat:
                Hit_Cat += 1
                hits += HitIndices_cat
            if len(HitIndices_tcga) == 1 and not 5 in HitIndices_tcga:
                Hit_TCGA += 1
                hits += HitIndices_tcga
            if len(HitIndices_tsinghwa) == 1 and not 5 in HitIndices_tsinghwa:
                Hit_Tsinghwa += 1
                hits += HitIndices_tsinghwa
            if Hit_Cat + Hit_TCGA + Hit_Tsinghwa >= 2:
                finalHit = filter(lambda x: hits.count(x) >= 2, set(hits))
                if len(finalHit) != 1:
                    changedDic[pre + '_' + arm] = 0
                    continue
                print "5' end", pre, arm, finalHit[0] - 5
                print '\t', "Support by Cat:", Hit_Cat
                print '\t', "Support by TCGA:", Hit_TCGA
                print '\t', "Support by Tsinghua", Hit_Tsinghwa

                HitIndex_cat = HitIndex_tcga = HitIndex_tsinghwa = finalHit[0]
                changedDic[pre + '_' + arm] = finalHit[0] - 5
                newExps_cat = []; newExps_tcga = []; newExps_tsinghwa = []
                newMedian_cat = []; newMedian_tcga = []; newMedian_tsinghwa = []
                newDOffset_cat = []; newDOffset_tcga = []; newDOffset_tsinghwa = []
                adjustFactor = 5 - HitIndex_cat
                for i in xrange(11):
                    if i - adjustFactor < 0 or i - adjustFactor > 10:
                        newExps_cat.append([0 for y in xrange(len(rcDic_cat[pre][arm][0]))])
                        newExps_tcga.append([0 for y in xrange(len(rcDic_tcga[pre][arm][0]))])
                        newExps_tsinghwa.append([0 for y in xrange(len(rcDic_tsinghwa[pre][arm][0]))])
                        newMedian_cat.append([0,0,0])
                        newMedian_tcga.append([0,0,0])
                        newMedian_tsinghwa.append([0,0,0])
                        newDOffset_cat.append([[0 for y in range(11)] for z in xrange(len(rcDic_cat[pre][arm][0]))])
                        newDOffset_tcga.append([[0 for y in range(11)] for z in xrange(len(rcDic_tcga[pre][arm][0]))])
                        newDOffset_tsinghwa.append([[0 for y in range(11)] for z in xrange(len(rcDic_tsinghwa[pre][arm][0]))])
                    else:
                        newExps_cat.append(rcDic_cat[pre][arm][i-adjustFactor])
                        newExps_tcga.append(rcDic_tcga[pre][arm][i-adjustFactor])
                        newExps_tsinghwa.append(rcDic_tsinghwa[pre][arm][i-adjustFactor])
                        newMedian_cat.append(medianRcDic_cat[pre][arm][i-adjustFactor])
                        newMedian_tcga.append(medianRcDic_tcga[pre][arm][i-adjustFactor])
                        newMedian_tsinghwa.append(medianRcDic_tsinghwa[pre][arm][i-adjustFactor])
                        newDOffset_cat.append(dOffsetDic_cat[pre][arm][i-adjustFactor])
                        newDOffset_tcga.append(dOffsetDic_tcga[pre][arm][i-adjustFactor])
                        newDOffset_tsinghwa.append(dOffsetDic_tsinghwa[pre][arm][i-adjustFactor])
                rcDic_cat[pre][arm] = newExps_cat
                rcDic_tcga[pre][arm] = newExps_tcga
                rcDic_tsinghwa[pre][arm] = newExps_tsinghwa
                medianRcDic_cat[pre][arm] = newMedian_cat
                medianRcDic_tcga[pre][arm] = newMedian_tcga
                medianRcDic_tsinghwa[pre][arm] = newMedian_tsinghwa
                dOffsetDic_cat[pre][arm] = newDOffset_cat
                dOffsetDic_tcga[pre][arm] = newDOffset_tcga
                dOffsetDic_tsinghwa[pre][arm] = newDOffset_tsinghwa
            else:
                changedDic[pre + '_' + arm] = 0
                continue
    return rcDic_cat, rcDic_tcga, rcDic_tsinghwa,medianRcDic_cat, medianRcDic_tcga, medianRcDic_tsinghwa, dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa, changedDic

def adjust_tpOffset(dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa, fpChangedDic):
    changedDic = dict()

    for pre in dOffsetDic_cat.keys():
        if not dOffsetDic_tcga.has_key(pre): continue
        if not dOffsetDic_tsinghwa.has_key(pre): continue
        for arm in dOffsetDic_cat[pre].keys():
            if not dOffsetDic_tcga[pre].has_key(arm): continue
            if not dOffsetDic_tsinghwa[pre].has_key(arm): continue
            zeroCat = 0; zeroTCGA = 0; zeroTsinghwa = 0
            maxIndices_cat = []
            for i in xrange(len(dOffsetDic_cat[pre][arm][0])):
                exps = dOffsetDic_cat[pre][arm][5][i]
                if max(exps) == 0:
                    zeroCat += 1
                    continue
                maxIndices_cat += [index for index, value in enumerate(exps) if value == max(exps)]
            maxIndices_tcga = []
            for i in xrange(len(dOffsetDic_tcga[pre][arm][0])):
                exps = dOffsetDic_tcga[pre][arm][5][i]
                if max(exps) == 0:
                    zeroTCGA += 1
                    continue
                maxIndices_tcga += [index for index, value in enumerate(exps) if value == max(exps)]
            maxIndices_tsinghwa = []
            for i in xrange(len(dOffsetDic_tsinghwa[pre][arm][0])):
                exps = dOffsetDic_tsinghwa[pre][arm][5][i]
                if max(exps) == 0:
                    zeroTsinghwa += 1
                    continue
                maxIndices_tsinghwa += [index for index, value in enumerate(exps) if value == max(exps)]

            HitIndices_cat = []
            if len(dOffsetDic_cat[pre][arm][0]) - zeroCat != 0:
                for i in xrange(len(dOffsetDic_cat[pre][arm][0][0])):
                    if maxIndices_cat.count(i) >= (len(dOffsetDic_cat[pre][arm][0]) - zeroCat) * 0.9:
                        HitIndices_cat.append(i)
            HitIndices_tcga = []
            if len(dOffsetDic_tcga[pre][arm][0]) - zeroTCGA != 0:
                for i in xrange(len(dOffsetDic_tcga[pre][arm][0][0])):
                    if maxIndices_tcga.count(i) >= (len(dOffsetDic_tcga[pre][arm][0]) - zeroTCGA) * 0.9:
                        HitIndices_tcga.append(i)
            HitIndices_tsinghwa = []
            if len(dOffsetDic_tsinghwa[pre][arm][0]) - zeroTsinghwa != 0:
                for i in xrange(len(dOffsetDic_tsinghwa[pre][arm][0][0])):
                    if maxIndices_tsinghwa.count(i) >= (len(dOffsetDic_tsinghwa[pre][arm][0]) - zeroTCGA) * 0.9:
                        HitIndices_tsinghwa.append(i)

            Hit_Cat = 0; Hit_TCGA = 0; Hit_Tsinghwa = 0
            hits = []
            if len(HitIndices_cat) == 1 and not 5 in HitIndices_cat:
                Hit_Cat += 1
                hits += HitIndices_cat
            if len(HitIndices_tcga) == 1 and not 5 in HitIndices_tcga:
                Hit_TCGA += 1
                hits += HitIndices_tcga
            if len(HitIndices_tsinghwa) == 1 and not 5 in HitIndices_tsinghwa:
                Hit_Tsinghwa += 1
                hits += HitIndices_tsinghwa
            if Hit_Cat + Hit_TCGA + Hit_Tsinghwa >= 2:
                finalHit = filter(lambda x: hits.count(x) >= 2, set(hits))
                if len(finalHit) != 1:
                    changedDic[pre + '_' + arm] = 0
                    continue
                print "3' end", pre, arm, finalHit[0] - 5
                print '\t', "Support by Cat:", Hit_Cat
                print '\t', "Support by TCGA:", Hit_TCGA
                print '\t', "Support by Tsinghua", Hit_Tsinghwa
                HitIndex_cat = HitIndex_tcga = HitIndex_tsinghwa = finalHit[0]
                changedDic[pre + '_' + arm] = finalHit[0] - 5
                newDOffset_cat = []; newDOffset_tcga = []; newDOffset_tsinghwa = []
                adjustFactor = 5 - HitIndex_cat
                for fpOffset in xrange(11):
                    tpExpsAllSamples_cat = []
                    for sampleN in xrange(len(dOffsetDic_cat[pre][arm][0])):
                        tpExps = []
                        for tpOffset in xrange(11):
                            if tpOffset - adjustFactor < 0 or tpOffset - adjustFactor > 10:
                                tpExps.append(0)
                            else:
                                tpExps.append(dOffsetDic_cat[pre][arm][fpOffset][sampleN][tpOffset - adjustFactor])
                        tpExpsAllSamples_cat.append(tpExps)
                    newDOffset_cat.append(tpExpsAllSamples_cat)
                    
                    tpExpsAllSamples_tcga = []
                    for sampleN in xrange(len(dOffsetDic_tcga[pre][arm][0])):
                        tpExps = []
                        for tpOffset in xrange(11):
                            if tpOffset - adjustFactor < 0 or tpOffset - adjustFactor > 10:
                                tpExps.append(0)
                            else:
                                tpExps.append(dOffsetDic_tcga[pre][arm][fpOffset][sampleN][tpOffset - adjustFactor])
                        tpExpsAllSamples_tcga.append(tpExps)
                    newDOffset_tcga.append(tpExpsAllSamples_tcga)

                    tpExpsAllSamples_tsinghwa = []
                    for sampleN in xrange(len(dOffsetDic_tsinghwa[pre][arm][0])):
                        tpExps = []
                        for tpOffset in xrange(11):
                            if tpOffset - adjustFactor < 0 or tpOffset - adjustFactor > 10:
                                tpExps.append(0)
                            else:
                                tpExps.append(dOffsetDic_tsinghwa[pre][arm][fpOffset][sampleN][tpOffset - adjustFactor])
                        tpExpsAllSamples_tsinghwa.append(tpExps)
                    newDOffset_tsinghwa.append(tpExpsAllSamples_tsinghwa)

                dOffsetDic_cat[pre][arm] = newDOffset_cat
                dOffsetDic_tcga[pre][arm] = newDOffset_tcga
                dOffsetDic_tsinghwa[pre][arm] = newDOffset_tsinghwa
            else:
                    changedDic[pre + '_' + arm] = 0
                    continue
    return dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa, changedDic

def get_MajorIsoformIndex(medianReads, isoformType):
    # get major isomiR
    tempList = []
    for i in range(11):
        if i == 5: continue
        tempList.append([i,medianReads[i][0]])
    tempList.sort(key=lambda row:row[1])
    if isoformType == 'first':
        majorIndex = tempList[-1][0]
    else:
        majorIndex = tempList[-2][0]
    return majorIndex

def get_AlleleFrequency(rcDic, rcDic2, cutoff, dataset, isoformType):
    afDic_for_compare_offset = dict(); afDic_for_compare_arm = dict()
    afDic_for_compare_frequency = dict(); afDic_normal = dict(); afDic_cancer = dict()
    afDic_all_sd = dict(); afDic_normal_sd = dict(); afDic_cancer_sd = dict()
    afDic_for_compare_frequency_up = dict(); afDic_for_compare_frequency_down = dict()
    anovaDic = dict(); ranksumUpDic = dict(); ranksumDownDic = dict()
    for pre in rcDic.keys():
        for arm in rcDic[pre].keys():
            reads = rcDic[pre][arm] # median
            majorIndex = get_MajorIsoformIndex(reads, isoformType) 
            reads2 = rcDic2[pre][arm] # all reads {pre:{arm:[[exps of offsets], ...]}}
            
            #if reads[4][1] + reads[5][1] + reads[6][1] >= cutoff or reads[4][2] + reads[5][2] + reads[6][2] >= cutoff:
            if reads[majorIndex][1] + reads[5][1] >= cutoff or reads[majorIndex][2] + reads[5][2] >= cutoff:
                if isoformType != 'first':
                    if reads[majorIndex][1] < 10 and reads[majorIndex][2] < 10: continue
                afDic_for_compare_offset[pre+'\t'+arm+'\t'+'-1'] = (reads[4][0]) / float(reads[4][0]+reads[5][0]+reads[6][0])
                afDic_for_compare_offset[pre+'\t'+arm+'\t'+'+1'] = (reads[6][0]) / float(reads[4][0]+reads[5][0]+reads[6][0])
                #temp_frequencies = map(lambda x: math.log((reads2[4][x] + reads2[6][x] + 0.2) / float(reads2[4][x] + reads2[5][x] + reads2[6][x] + 0.3),10), range(len(reads2[0])))
                temp_frequencies = map(lambda x: math.log((reads2[majorIndex][x] + 0.1) / float(reads2[5][x] + 0.1), 10), range(len(reads2[0])))
                if dataset == 'Catholic':
                    temp_frequencies_N = temp_frequencies[:15]
                    temp_frequencies_C = temp_frequencies[15:]
                elif dataset == "TCGA":
                    temp_frequencies_N = temp_frequencies[:50]
                    temp_frequencies_C = temp_frequencies[50:]
                elif dataset == 'Tsinghua':
                    temp_frequencies_N = temp_frequencies[:20]
                    temp_frequencies_C = temp_frequencies[20:]

                afDic_for_compare_frequency[pre+'\t'+arm] = np.mean(temp_frequencies)
                afDic_normal[pre+'\t'+arm] = np.mean(temp_frequencies_N)
                afDic_cancer[pre+'\t'+arm] = np.mean(temp_frequencies_C)

                afDic_all_sd[pre+'\t'+arm] = np.std(temp_frequencies)
                afDic_normal_sd[pre+'\t'+arm] = np.std(temp_frequencies_N)
                afDic_cancer_sd[pre+'\t'+arm] = np.std(temp_frequencies_C)

                anovaDic[pre+'\t'+arm] = stats.f_oneway(temp_frequencies_N,temp_frequencies_C)[1]
                ranksumUpDic[pre+'\t'+arm] = stats.mannwhitneyu(temp_frequencies_N,temp_frequencies_C,alternative="greater")[1]
                ranksumDownDic[pre+'\t'+arm] = stats.mannwhitneyu(temp_frequencies_N,temp_frequencies_C,alternative="less")[1]
                #afDic_for_compare_frequency_up[pre+'\t'+arm] = (reads[4][0]+0.1) / float(reads[5][0]+0.1)
                #afDic_for_compare_frequency_down[pre+'\t'+arm] = (reads[6][0]+0.1) / float(reads[5][0]+0.1)
                if majorIndex < 5:
                    afDic_for_compare_frequency_up[pre+'\t'+arm] = (reads[majorIndex][0]+0.1) / float(reads[5][0]+0.1)
                elif majorIndex > 5:
                    afDic_for_compare_frequency_down[pre+'\t'+arm] = (reads[majorIndex][0]+0.1) / float(reads[5][0]+0.1)
        if rcDic[pre].has_key('FP') and rcDic[pre].has_key('TP'):
            reads_FP = rcDic[pre]['FP']
            reads_TP = rcDic[pre]['TP']
    anovaFdrDic = adjust_pvalues(anovaDic)
    ranksumUpFdrDic = adjust_pvalues(ranksumUpDic)
    ranksumDownFdrDic = adjust_pvalues(ranksumDownDic)

    ranksumFdrDic = dict()
    for id in ranksumUpFdrDic.keys():
        ranksumFdrDic[id] = min(ranksumUpFdrDic[id], ranksumDownFdrDic[id])

    return afDic_for_compare_offset, afDic_for_compare_arm, afDic_for_compare_frequency, afDic_for_compare_frequency_up, afDic_for_compare_frequency_down, afDic_normal, afDic_cancer, afDic_all_sd, afDic_normal_sd, afDic_cancer_sd, anovaFdrDic, ranksumFdrDic

def intersect_Frequency(catAfDic_for_compare_frequency, tcgaAfDic_for_compare_frequency, tsinghwaAfDic_for_compare_frequency, cat_corr_A, tcga_corr_A, tsinghwa_corr_A):
    intersectFrequencyDic = {'FP':{'Low':[],'Mid':[],'High':[],'other':[]},'TP':{'Low':[],'Mid':[],'High':[],'other':[]}}

    # calculate median values among highly correlated miRNAs in both Catholic and Tsinghwa
    #cat_FP = np.mean(map(lambda x: catAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'FP' and y in cat_corr_A and y in tsinghwa_corr_A, catAfDic_for_compare_frequency.keys())))#; rv_cat_FP = fw_cat_FP
    cat_FP = np.mean(map(lambda x: catAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'FP', catAfDic_for_compare_frequency.keys())))
    #cat_TP = np.mean(map(lambda x: catAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'TP' and y in cat_corr_A and y in tsinghwa_corr_A, catAfDic_for_compare_frequency.keys())))#; rv_cat_TP = fw_cat_TP
    cat_TP = np.mean(map(lambda x: catAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'TP', catAfDic_for_compare_frequency.keys())))
    #tcga_FP = np.mean(map(lambda x: tcgaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'FP' and y in tcga_corr_A, tcgaAfDic_for_compare_frequency.keys())))#; rv_tcga_FP = fw_tcga_FP
    tcga_FP = np.mean(map(lambda x: tcgaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'FP', tcgaAfDic_for_compare_frequency.keys())))
    #tcga_TP = np.mean(map(lambda x: tcgaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'TP' and y in tcga_corr_A, tcgaAfDic_for_compare_frequency.keys())))#; rv_tcga_TP = fw_tcga_TP
    tcga_TP = np.mean(map(lambda x: tcgaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'TP', tcgaAfDic_for_compare_frequency.keys())))
    #tsinghwa_FP = np.mean(map(lambda x: tsinghwaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'FP' and y in tsinghwa_corr_A and y in cat_corr_A, tsinghwaAfDic_for_compare_frequency.keys())))#; rv_tsinghwa_FP = fw_tsinghwa_FP
    tsinghwa_FP = np.mean(map(lambda x: tsinghwaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'FP', tsinghwaAfDic_for_compare_frequency.keys())))
    #tsinghwa_TP = np.mean(map(lambda x: tsinghwaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'TP' and y in tsinghwa_corr_A and y in cat_corr_A, tsinghwaAfDic_for_compare_frequency.keys())))#; rv_tsinghwa_TP = fw_tsinghwa_TP
    tsinghwa_TP = np.mean(map(lambda x: tsinghwaAfDic_for_compare_frequency[x], filter(lambda y: y.split('\t')[1] == 'TP', tsinghwaAfDic_for_compare_frequency.keys())))

    intersectedNames = set(catAfDic_for_compare_frequency.keys()) | set(tsinghwaAfDic_for_compare_frequency.keys())
    hit = 0
    for name in intersectedNames:
        pre, arm =  name.split('\t')
        if arm == 'FP':
            cat_Med = cat_FP
            tcga_Med = tcga_FP
            tsinghwa_Med = tsinghwa_FP
        else:
            cat_Med = cat_TP
            tcga_Med = tcga_TP
            tsinghwa_Med = tsinghwa_TP
        # compare Catholic and Tsinghwa because TCGA is not well correlated with other data
        # only correlated miRNAs in both Catholic and Tsinghwa ? 
        #if catAfDic_for_compare_frequency.has_key(name) and tsinghwaAfDic_for_compare_frequency.has_key(name) and name in cat_corr_A and name in tsinghwa_corr_A:
        if catAfDic_for_compare_frequency.has_key(name) and tsinghwaAfDic_for_compare_frequency.has_key(name):
            if catAfDic_for_compare_frequency[name] < cat_Med and tsinghwaAfDic_for_compare_frequency[name] < tsinghwa_Med:
                intersectFrequencyDic[arm]['Low'].append(pre)
                hit += 1
            elif cat_Med <= catAfDic_for_compare_frequency[name] < cat_Med and tsinghwa_Med <= tsinghwaAfDic_for_compare_frequency[name] < tsinghwa_Med:
                intersectFrequencyDic[arm]['Mid'].append(pre)
                hit += 1
            elif cat_Med <= catAfDic_for_compare_frequency[name] and tsinghwa_Med <= tsinghwaAfDic_for_compare_frequency[name]:
                intersectFrequencyDic[arm]['High'].append(pre)
                hit += 1
            else:
                intersectFrequencyDic[arm]['other'].append(pre)
                hit += 1
        else:
            intersectFrequencyDic[arm]['other'].append(pre)
            hit += 1
    for arm in intersectFrequencyDic.keys():
        print arm
        for trend in intersectFrequencyDic[arm].keys():
            print trend, len(intersectFrequencyDic[arm][trend])
    return intersectFrequencyDic

def Calc_Correlation(rcDic, rcDic2, totalsamples, afDic, isoformType):
    correlation = []
    for pre in rcDic2.keys():
        for arm in rcDic2[pre].keys():
            if not afDic.has_key(pre+'\t'+arm): continue
            medianReads = rcDic[pre][arm]
            majorIndex = get_MajorIsoformIndex(medianReads, isoformType)
            offsetList = []; canoList = []
            for i in xrange(len(totalsamples)):
                sample = totalsamples[i]
                mature = rcDic2[pre][arm][5][i]
                offsetList.append(rcDic2[pre][arm][majorIndex][i])
                canoList.append(rcDic2[pre][arm][5][i])
            corrR, corrP = stats.pearsonr(offsetList, canoList)
            if pre == 'hsa-mir-1307' and arm == 'FP':
                print corrR, corrP
            if corrP <= 0.05:
                correlation.append(pre + '\t' + arm)
    return correlation

def write_rcDic(expDic, rcDic, rcDic2, outputD, totalsamples, afDic, analysistype, afDic2, infoDic, armDic2, changedDic, isoformType):
    correlation_N = []; correlation_C = []; correlation_A = []; pDic = dict()
    writer = open(outputD + 'AlleleReadCountDist_OffsetPlustMinus1_'+analysistype+'.txt','w')
    #writer_P = open(outputD + 'AlleleReadCountDist_OffsetMinus1_'+analysistype+'.txt','w')
    #writer_M = open(outputD + 'AlleleReadCountDist_OffsetPlus1_'+analysistype+'.txt','w')
    writerFP = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_'+analysistype+'_FP.txt','w')
    writerTP = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_'+analysistype+'_TP.txt','w')
    writerFPExp = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_Exp_FP.txt','w')
    writerFPRPM = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_RPM_FP.txt','w')
    #writerFPRPM_P = open(outputD + 'AlleleFrequencyDistAll_OffsetPlus1_RPM_FP.txt','w')
    #writerFPRPM_M = open(outputD + 'AlleleFrequencyDistAll_OffsetMinus1_RPM_FP.txt','w')
    writerMatureFP = open(outputD + 'AlleleFrequencyDistAll_Mature_RPM_FP.txt','w')
    #writerFP_P = open(outputD + 'AlleleFrequencyDistAll_OffsetPlus1_'+analysistype+'_FP.txt','w')
    #writerFP_M = open(outputD + 'AlleleFrequencyDistAll_OffsetMinus1_'+analysistype+'_FP.txt','w')
    #writerTP_P = open(outputD + 'AlleleFrequencyDistAll_OffsetPlus1_' +analysistype+'_TP.txt','w')
    #writerTP_M = open(outputD + 'AlleleFrequencyDistAll_OffsetMinus1_' + analysistype+'_TP.txt','w')
    writerTPExp = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_Exp_TP.txt','w')
    writerTPRPM = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_RPM_TP.txt','w')
    #writerTPRPM_P = open(outputD + 'AlleleFrequencyDistAll_OffsetPlus1_RPM_TP.txt','w')
    #writerTPRPM_M = open(outputD + 'AlleleFrequencyDistAll_OffsetMinus1_RPM_TP.txt','w')
    writerMatureTP = open(outputD + 'AlleleFrequencyDistAll_Mature_RPM_TP.txt','w')
    writerMatureRPM = open(outputD + 'AlleleFrequencyDistAll_Mature_RPM.txt','w')
    writerOffsetRPM = open(outputD + 'AlleleFrequencyDistAll_OffsetPlus1Minus1_RPM.txt','w')
    writerMatureExp = open(outputD + 'AlleleFrequencyDistAll_Mature_Exp.txt','w')
    writerOffsetExp = open(outputD + 'AlleleFrequencyDistAll_OffsetPlus1Minus1_Exp.txt','w')
    writerAll = open(outputD + 'AlleleFrequencyDistAll_OffsetPlusMinus1_'+analysistype+'_All.txt','w')
    writerProportion = open(outputD + 'AlleleReadCountDist_TotalOffset_FP.txt','w')
    writerFP.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerFPExp.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerFPRPM.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerFPRPM_P.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerFPRPM_M.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerMatureFP.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerFP_P.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerFP_M.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerTP_P.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerTP_M.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerTP.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerTPExp.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerTPRPM.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerTPRPM_P.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    #writerTPRPM_M.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerMatureTP.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerMatureRPM.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerOffsetRPM.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerMatureExp.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerOffsetExp.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerAll.write('miRNA' + '\t' + '\t'.join(totalsamples) + '\n')
    writerProportion.write('miRNA' + '\t' + '\t'.join(map(str,range(-5,6))) + '\n')
    writer.write('Precursor' + '\t' + 'Arm' + '\t' + 'Sample' + '\t' + 'OffsetReads' + '\t' + 'MatureReads' + '\t' + 'CorrelationR' + '\t' + 'CorrelationP' + '\n')
    for arm in afDic.keys():
        for trend in afDic[arm].keys():
            for pre in afDic[arm][trend]:
                if not afDic2.has_key(pre+'\t'+arm): continue

                majorIndex = get_MajorIsoformIndex(rcDic[pre][arm], isoformType)
                # write proportion 
                tempSum = sum(map(lambda x: x[0], rcDic[pre][arm]))
                if arm == "FP":
                    writerProportion.write(pre + '\t' + '\t'.join(map(lambda x: str(x[0] / float(tempSum)), rcDic[pre][arm])) + '\n')
                matureName = armDic2[pre][arm]
                #if changedDic[pre + '_' + arm] == 0:
                    #offsetPName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == "1", filter(lambda x: x.count('_') > 0, expDic.keys()))
                    #offsetMName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == "-1", filter(lambda x: x.count('_') > 0, expDic.keys()))
                #    offsetName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(majorIndex), filter(lambda x: x.count('_') > 0, expDic.keys()))
                #else:
                #    if changedDic[pre + '_' + arm] - 1 == 0:
                        #offsetMName = [matureName]
                        #offsetPName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(1 + changedDic[pre + '_' + arm]), filter(lambda x: x.count('_') > 0, expDic.keys()))
                #        print "Here !!!!!!!!!!!!!!!!!!", pre, arm,offsetPName, changedDic[pre+'_'+arm]
                #    elif changedDic[pre + '_' + arm] + 1 == 0:
                #        offsetPName = [matureName]
                #        offsetMName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(-1 + changedDic[pre + '_' + arm]), filter(lambda x: x.count('_') > 0, expDic.keys()))
                #        print "Here !!!!!!!!!!!!!!!!!", pre,arm,offsetMName,changedDic[pre+'_'+arm]
                #    else:
                #        offsetPName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(1 + changedDic[pre + '_' + arm]), filter(lambda x: x.count('_') > 0, expDic.keys()))
                #        offsetMName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(-1 + changedDic[pre + '_' + arm]), filter(lambda x: x.count('_') > 0, expDic.keys()))
                #    matureName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(changedDic[pre + '_' + arm]), filter(lambda x: x.count('_') > 0, expDic.keys()))[0]
                if changedDic[pre+'_'+arm] + majorIndex -5 == 0:
                    offsetName = [matureName]
                else:
                    offsetName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(changedDic[pre+'_'+arm] + majorIndex - 5), filter(lambda x: x.count('_') > 0, expDic.keys()))
                if changedDic[pre+'_'+arm]!=0:
                    matureName = filter(lambda x: x.split('_')[0] == matureName and x.split('_')[2] == str(changedDic[pre+'_'+arm]), filter(lambda x: x.count('_') > 0, expDic.keys()))[0]
                frequencyList = []#; frequencyList_M = []; frequencyList_P = []
                expList_mature = []; expList_isomiR = []
                rpmList_mature = []; rpmList_isomiR = []#; rpmList_M = []; rpmList_P = []; rpmList_MP = []
                offsetList_N = []; canoList_N = []
                #offsetList_N_M = []#; canoList_N_M = []
                #offsetList_N_P = []#; canoList_N_P = []
                offsetList_C = []; canoList_C = []
                #offsetList_C_M = []#; canoList_C_M = []
                #offsetList_C_P = []#; canoList_C_P = []
                print(rcDic2[pre][arm][5])
                print(rcDic2[pre][arm][majorIndex])
                for i in xrange(len(totalsamples)):
                    sample = totalsamples[i]
                    #if len(rcDic2[pre][arm][4]) == 0:
                    #    m1 = 0
                    #else:
                    #    m1 = rcDic2[pre][arm][4][i]
                    #if len(rcDic2[pre][arm][6]) == 0:
                    #    p1 = 0
                    #else:
                    #    p1 = rcDic2[pre][arm][6][i]
                    isoform = rcDic2[pre][arm][majorIndex][i]
                    mature = rcDic2[pre][arm][5][i]
                    #frequencyList.append(math.log((m1+p1+0.1)/float(mature+0.1),10))
                    frequencyList.append(math.log((isoform+0.1)/float(mature+0.1),10))
                    #frequencyList_M.append(math.log((m1+0.1)/float(mature+0.1),10))
                    #frequencyList_P.append(math.log((p1+0.1)/float(mature+0.1),10))
                    expList_isomiR.append(isoform)
                    expList_mature.append(mature)
                    rpmList_mature.append(expDic[matureName][i])
                    #if len(offsetMName) > 0:
                    #    rpmList_M.append(expDic[offsetMName[0]][i])
                    #if len(offsetPName) > 0:
                    #    rpmList_P.append(expDic[offsetPName[0]][i])
                    #if len(offsetMName) > 0 and len(offsetPName) > 0:
                    #    rpmList_MP.append(expDic[offsetMName[0]][i] + expDic[offsetPName[0]][i])
                    #elif len(offsetMName) > 0:
                    #    rpmList_MP.append(expDic[offsetMName[0]][i])
                    #elif len(offsetMName) > 0:
                    #    rpmList_MP.append(expDic[offsetPName[0]][i])
                    if len(offsetName) == 0:
                        rpmList_isomiR.append(0)
                    else:
                        rpmList_isomiR.append(expDic[offsetName[0]][i])

                    if infoDic[sample] == 'Normal':
                        #offsetList_N.append(m1+p1)
                        #offsetList_N_M.append(m1)
                        #offsetList_N_P.append(p1)
                        offsetList_N.append(isoform)
                        canoList_N.append(mature)
                    else:
                        #offsetList_C.append(m1+p1)
                        #offsetList_C_M.append(m1)
                        #offsetList_C_P.append(p1)
                        offsetList_C.append(isoform)
                        canoList_C.append(mature)
                
                corrR_N, corrP_N = stats.pearsonr(offsetList_N, canoList_N)
                corrR_C, corrP_C = stats.pearsonr(offsetList_C, canoList_C)
                corrR_A, corrP_A = stats.pearsonr(offsetList_N+offsetList_C, canoList_N+canoList_C)

                if corrP_N <= 0.05:
                    correlation_N.append(pre + '\t' + arm)
                if corrP_C <= 0.05:
                    correlation_C.append(pre + '\t' + arm)
                if corrP_A <= 0.05:
                    correlation_A.append(pre + '\t' + arm)
                pDic[pre + '\t' + arm] = corrP_A
                if arm == 'FP':
                    writerFP.write(pre + '\t' + '\t'.join(map(str,frequencyList)) + '\n')
                    #writerFP_M.write(pre + '\t' + '\t'.join(map(str,frequencyList_M)) + '\n')
                    #writerFP_P.write(pre + '\t' + '\t'.join(map(str,frequencyList_P)) + '\n')
                    writerFPExp.write(pre + '\t' + '\t'.join(map(str,expList_isomiR)) + '\n')
                    writerMatureFP.write(pre + '\t' + '\t'.join(map(str,rpmList_mature)) + '\n')
                    writerFPRPM.write(pre + '\t' + '\t'.join(map(str,rpmList_isomiR)) + '\n')
                    writerMatureRPM.write(pre + '_FP' + '\t' + '\t'.join(map(str,rpmList_mature)) + '\n')
                    writerMatureExp.write(pre + '_FP' + '\t' + '\t'.join(map(str,expList_mature)) + '\n')
                    writerOffsetExp.write(pre + '_FP' + '\t' + '\t'.join(map(str,expList_isomiR)) + '\n')
                    writerAll.write(pre + '_FP' + '\t' + '\t'.join(map(str,frequencyList)) + '\n')
                    if len(rpmList_isomiR) > 0:
                        writerOffsetRPM.write(pre + '_FP' + '\t' + '\t'.join(map(str,rpmList_isomiR)) + '\n')
                    #if len(rpmList_M) > 0 and sum(rpmList_M) > 0:
                    #    writerFPRPM_M.write(pre + '\t' + '\t'.join(map(str,rpmList_M)) + '\n')
                        
                    #if len(rpmList_P) > 0 and sum(rpmList_P) > 0:
                    #    writerFPRPM_P.write(pre + '\t' + '\t'.join(map(str,rpmList_P)) + '\n')
                else:
                    writerTP.write(pre + '\t' + '\t'.join(map(str,frequencyList)) + '\n')
                    #writerTP_M.write(pre + '\t' + '\t'.join(map(str,frequencyList_M)) + '\n')
                    #writerTP_P.write(pre + '\t' + '\t'.join(map(str,frequencyList_P)) + '\n')
                    writerTPExp.write(pre + '\t' + '\t'.join(map(str,expList_isomiR)) + '\n')
                    writerMatureTP.write(pre + '\t' + '\t'.join(map(str,rpmList_mature)) + '\n')
                    writerTPRPM.write(pre + '\t' + '\t'.join(map(str,rpmList_isomiR)) + '\n')
                    writerMatureRPM.write(pre + '_TP' + '\t' + '\t'.join(map(str,rpmList_mature)) + '\n')
                    writerMatureExp.write(pre + '_TP' + '\t' + '\t'.join(map(str,expList_mature)) + '\n')
                    writerOffsetExp.write(pre + '_TP' + '\t' + '\t'.join(map(str,expList_isomiR)) + '\n')
                    writerAll.write(pre + '_TP' + '\t' + '\t'.join(map(str,frequencyList)) + '\n')
                    if len(rpmList_isomiR) > 0:
                        writerOffsetRPM.write(pre + '_TP' + '\t' + '\t'.join(map(str,rpmList_isomiR)) + '\n')
                    #if len(rpmList_M) > 0 and sum(rpmList_M) > 0:
                    #    writerTPRPM_M.write(pre + '\t' + '\t'.join(map(str,rpmList_M)) + '\n')
                    #if len(rpmList_P) > 0 and sum(rpmList_P) > 0:
                    #    writerTPRPM_P.write(pre + '\t' + '\t'.join(map(str,rpmList_P)) + '\n')

                for i in xrange(len(totalsamples)):
                    sample = totalsamples[i]
                    offset = (offsetList_N+offsetList_C)[i]
                    #offset_M = (offsetList_N_M+offsetList_C_M)[i]
                    #offset_P = (offsetList_N_P+offsetList_C_P)[i]
                    canonical = (canoList_N+canoList_C)[i]
                    writer.write(pre + '\t' + arm + '\t' + sample + '\t' + str(offset) + '\t' + str(canonical) + '\t' + str(corrR_A) + '\t' + str(corrP_A) + '\n')
                    #writer_M.write(pre + '\t' + arm + '\t' + sample + '\t' + str(offset_M) + '\t' + str(canonical) + '\t' + str(corrR_A) + '\t' + str(corrP_A) + '\n')
                    #writer_P.write(pre + '\t' + arm + '\t' + sample + '\t' + str(offset_P) + '\t' + str(canonical) + '\t' + str(corrR_A) + '\t' + str(corrP_A) + '\n')
    writer.close()
    #writer_M.close()
    #writer_P.close()
    writerFP.close()
    writerFPExp.close()
    writerFPRPM.close()
    #writerFPRPM_P.close()
    #writerFPRPM_M.close()
    writerMatureFP.close()
    #writerFP_M.close()
    #writerFP_P.close()
    writerTP.close()
    writerTPExp.close()
    writerTPRPM.close()
    #writerTPRPM_P.close()
    #writerTPRPM_M.close()
    writerMatureTP.close()
    #writerTP_M.close()
    #writerTP_P.close()
    writerMatureRPM.close()
    writerOffsetRPM.close()
    writerMatureExp.close()
    writerOffsetExp.close()
    writerAll.close()
    writerProportion.close()
    return correlation_N, correlation_C, correlation_A, pDic

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

def write_afDic(catAfDic, catAfDic_up, catAfDic_down, tcgaAfDic, tcgaAfDic_up, tcgaAfDic_down, tsinghwaAfDic, tsinghwaAfDic_up, tsinghwaAfDic_down, outputD, analysistype, armDic2, matureSeqDic, hairpinDic, intersectFrequencyDic, changedDic, tpChangedDic, cat_corr_A, tcga_corr_A, tsinghwa_corr_A, medianRcDic_cat, medianRcDic_tcga, medianRcDic_tsinghwa, cat_pDic, tcga_pDic, tsinghwa_pDic, catAfDic_normal, catAfDic_cancer, catAfDic_all_sd, catAfDic_normal_sd, catAfDic_cancer_sd, tcgaAfDic_normal, tcgaAfDic_cancer, tcgaAfDic_all_sd, tcgaAfDic_normal_sd, tcgaAfDic_cancer_sd, tsinghwaAfDic_normal, tsinghwaAfDic_cancer, tsinghwaAfDic_all_sd, tsinghwaAfDic_normal_sd, tsinghwaAfDic_cancer_sd, cat_anovaDic, tcga_anovaDic, tsinghwa_anovaDic, cat_ranksumDic, tcga_ranksumDic, tsinghwa_ranksumDic):
    writer = open(outputD + 'AlleleFrequencyDist_OffsetPlusMiNum1_'+analysistype+'_IntersectSet_FDR.txt','w')
    writer.write('Precursor' + '\t' + 'Arm' + '\t' + 'Cat_Frequency' + '\t' + 'TCGA_Frequency' + '\t' + 'Tsinghwa_Frequency' + '\t' + 'CanonicalFPend' + '\t' + 'FrequencyTrend' + '\t' + 'OffsetFP_Changed' + '\t' + 'OffsetTP_Changed' + '\t' + 'Cat_Frequency_Up' + '\t' + 'Cat_Frequency_Down' + '\t' + 'TCGA_Frequency_Up' + '\t' + 'TCGA_Frequency_Down' + '\t' + 'Tsinghwa_Frequency_Up' + '\t' + 'Tsinghwa_Frequency_Down' + '\t' + 'Catholic_Correlation' + '\t' + 'TCGA_Correlation' + '\t' + 'Tsinghwa_Correlation' + '\t' + 'MedianRC_Cat' + '\t' + 'MedianRC_TCGA' + '\t' + 'MedianRC_Tsinghwa' + '\t' + 'P_cat' + '\t' + 'P_TCGA' + '\t' + 'P_Tsinghwa' + '\t' + 'Cat_Frequency_N' + '\t' + 'Cat_Frequency_C' + '\t' + 'TCGA_Frequency_N' + '\t' + 'TCGA_Frequency_C' + '\t' + 'Tsinghua_Frequency_N' + '\t' + 'Tsinghua_Frequency_C' + '\t' + 'Cat_Frequency_SD' + '\t' + 'Cat_Frequency_N_SD' + '\t' + 'Cat_Frequency_C_SD' + '\t' + 'TCGA_Frequency_SD' + '\t' + 'TCGA_Frequency_N_SD' + '\t' + 'TCGA_Frequency_C_SD' + '\t' + 'Tsinghua_Frequency_SD' + '\t' + 'Tsinghua_Frequency_N_SD' + '\t' + 'Tsinghua_Frequency_C_SD' + '\t' + 'OneWayAnova_Catholic' + '\t' + 'OneWayAnova_TCGA' + '\t' + 'OneWayAnova_Tsinghua' + '\t' + 'Ranksum_Catholic' + '\t' + 'Ranksum_TCGA' + '\t' + 'Ranksum_Tsinghua' + '\n')
    for arm in intersectFrequencyDic.keys():
        for trend in ['Low','High','other']:
            for pre in intersectFrequencyDic[arm][trend]:
                name = pre + '\t' + arm
                mature = armDic2[pre][arm]
                matureSeq = matureSeqDic[mature]

                canonicalFPend = matureSeq[0]
                if catAfDic.has_key(name):
                    cat = str(catAfDic[name])
                    cat_N = str(catAfDic_normal[name])
                    cat_C = str(catAfDic_cancer[name])
                    cat_SD = str(catAfDic_all_sd[name])
                    cat_N_SD = str(catAfDic_normal_sd[name])
                    cat_C_SD = str(catAfDic_cancer_sd[name])
                    if catAfDic_up.has_key(name):
                        cat_up = str(catAfDic_up[name])
                    else:
                        cat_up = 'NA'
                    if catAfDic_down.has_key(name):
                        cat_down = str(catAfDic_down[name])
                    else:
                        cat_down = 'NA'
                    cat_anova = str(cat_anovaDic[name])
                    cat_ranksum = str(cat_ranksumDic[name])
                else:
                    cat = 'NA'
                    cat_N = 'NA'
                    cat_C = 'NA'
                    cat_SD = 'NA'
                    cat_N_SD ='NA'
                    cat_C_SD = 'NA'
                    cat_up = 'NA'
                    cat_down = 'NA'
                    cat_anova = 'NA'
                    cat_ranksum = 'NA'
                if name in cat_corr_A:
                    cat_corr = 'T'
                else:
                    cat_corr = 'F'
                if tcgaAfDic.has_key(name):
                    tcga = str(tcgaAfDic[name])
                    tcga_N = str(tcgaAfDic_normal[name])
                    tcga_C = str(tcgaAfDic_cancer[name])
                    tcga_SD = str(tcgaAfDic_all_sd[name])
                    tcga_N_SD = str(tcgaAfDic_normal_sd[name])
                    tcga_C_SD = str(tcgaAfDic_cancer_sd[name])
                    if tcgaAfDic_up.has_key(name):
                        tcga_up = str(tcgaAfDic_up[name])
                    else:
                        tcga_up = 'NA'
                    if tcgaAfDic_down.has_key(name):
                        tcga_down = str(tcgaAfDic_down[name])
                    else:
                        tcga_down = 'NA'
                    tcga_anova = str(tcga_anovaDic[name])
                    tcga_ranksum = str(tcga_ranksumDic[name])
                else:
                    tcga = 'NA'
                    tcga_N = 'NA'
                    tcga_C = 'NA'
                    tcga_SD = 'NA'
                    tcga_N_SD = 'NA'
                    tcga_C_SD = 'NA'
                    tcga_up = 'NA'
                    tcga_down = 'NA'
                    tcga_anova = 'NA'
                    tcga_ranksum = 'NA'
                if name in tcga_corr_A:
                    tcga_corr = 'T'
                else:
                    tcga_corr = 'F'
                if tsinghwaAfDic.has_key(name):
                    tsinghwa = str(tsinghwaAfDic[name])
                    tsinghua_N = str(tsinghwaAfDic_normal[name])
                    tsinghua_C = str(tsinghwaAfDic_cancer[name])
                    tsinghua_SD = str(tsinghwaAfDic_all_sd[name])
                    tsinghua_N_SD = str(tsinghwaAfDic_normal_sd[name])
                    tsinghua_C_SD = str(tsinghwaAfDic_cancer_sd[name])
                    if tsinghwaAfDic_up.has_key(name):
                        tsinghwa_up = str(tsinghwaAfDic_up[name])
                    else:
                        tsinghwa_up = 'NA'
                    if tsinghwaAfDic_down.has_key(name):
                        tsinghwa_down = str(tsinghwaAfDic_down[name])
                    else:
                        tsinghwa_down = 'NA'
                    tsinghwa_anova = str(tsinghwa_anovaDic[name])
                    tsinghwa_ranksum = str(tsinghwa_ranksumDic[name])
                else:
                    tsinghwa = 'NA'
                    tsinghua_N = 'NA'
                    tsinghua_C = 'NA'
                    tsinghua_SD = 'NA'
                    tsinghua_N_SD = 'NA'
                    tsinghua_C_SD = 'NA'
                    tsinghwa_up = 'NA'
                    tsinghwa_down = 'NA'
                    tsinghwa_anova = 'NA'
                    tsinghwa_ranksum = 'NA'
                if name in tsinghwa_corr_A:
                    tsinghwa_corr = 'T'
                else:
                    tsinghwa_corr = 'F'
                if changedDic.has_key(pre+'_'+arm):
                    change = changedDic[pre+'_'+arm]
                else:
                    change = 0
                if tpChangedDic.has_key(pre + '_' + arm):
                    changedTP = tpChangedDic[pre+'_'+arm]
                else:
                    changedTP = 0
                if cat_pDic.has_key(pre+'\t'+arm):
                    cat_P = cat_pDic[pre+'\t'+arm]
                else:
                    cat_P = 'NA'
                if tcga_pDic.has_key(pre+'\t'+arm):
                    tcga_P = tcga_pDic[pre+'\t'+arm]
                else:
                    tcga_P = 'NA'
                if tsinghwa_pDic.has_key(pre+'\t'+arm):
                    tsinghwa_P = tsinghwa_pDic[pre+'\t'+arm]
                else:
                    tsinghwa_P = 'NA'
                writer.write(name + '\t' + cat + '\t' + tcga + '\t' + tsinghwa + '\t' + canonicalFPend + '\t' + trend + '\t' + str(change) + '\t' + str(changedTP) + '\t' + '\t'.join([cat_up,cat_down,tcga_up,tcga_down,tsinghwa_up,tsinghwa_down]) + '\t' + '\t'.join([cat_corr,tcga_corr, tsinghwa_corr]) + '\t' + '\t'.join(map(str,[medianRcDic_cat[pre][arm][4][2] + medianRcDic_cat[pre][arm][6][2], medianRcDic_tcga[pre][arm][4][2] + medianRcDic_tcga[pre][arm][6][2], medianRcDic_tsinghwa[pre][arm][4][2] + medianRcDic_tsinghwa[pre][arm][6][2]])) + '\t' + '\t'.join(map(str,[cat_P, tcga_P,tsinghwa_P])) + '\t' + '\t'.join(map(str, [cat_N,cat_C,tcga_N,tcga_C,tsinghua_N, tsinghua_C, cat_SD, cat_N_SD, cat_C_SD, tcga_SD, tcga_N_SD, tcga_C_SD, tsinghua_SD, tsinghua_N_SD, tsinghua_C_SD, cat_anova, tcga_anova, tsinghwa_anova, cat_ranksum, tcga_ranksum, tsinghwa_ranksum])) + '\n')
    writer.close()

def write_dOffset(afDic_for_compare_frequency, corr_A, dOffsetDic, totalsamples, infoDic, outputD):
    writer_sampleFrequency = open(outputD + '/Correlation_Bw_5p_3p_SampleFrequency.txt','w')
    writer_positionalExp = open(outputD + '/Correlation_Bw_5p_3p_PositionalRelativeExpression.txt','w')
    writer_positionalExp_FP = open(outputD + '/Correlation_Bw_5p_PositionalRelativeExpression.txt','w')
    writer_positionalExp_norm = open(outputD + '/Correlation_Bw_5p_3p_PositionalNormalizedRelativeExpression.txt','w')
    writer_DoubleOffset = open(outputD + '/DoubleOffset_Distribution.txt','w')

    writer_sampleFrequency.write('Precursor' + '\t' + 'Arm' + '\t' + 'Offset_FP' + '\t' + '\t'.join(map(str,range(-5,6))) + '\n')
    writer_positionalExp.write('Precursor' + '\t' + 'Arm' + '\t' + 'Offset_FP' + '\t' + 'Offset_TP' + '\t' + 'Median' + '\t' + 'STD' + '\n')
    writer_positionalExp_FP.write('Precursor' + '\t' + 'Arm' + '\t' + 'Offset_FP' + '\t' + 'Median' + '\t' + 'STD' + '\n')
    writer_positionalExp_norm.write('Precursor' + '\t' + 'Arm' + '\t' + 'Offset_FP' + '\t' + 'Offset_TP' + '\t' + 'Median' + '\t' + 'STD' + '\n')
    writer_DoubleOffset.write('Precursor' + '\t' + 'Arm' + '\t' + 'Sample' + '\t' + 'Offset_FP' + '\t' + 'Offset_TP' + '\t' + 'ReadCount' + '\n')
    for pre in dOffsetDic.keys():
        for arm in dOffsetDic[pre].keys():
            if not afDic_for_compare_frequency.has_key(pre + '\t' + arm): continue
            # normalization factor
            tpOffsetExp_norm = dOffsetDic[pre][arm][5]
            relativeMedianDic_FP = dict()

            sampleFrequency_FP = [0 for y in range(11)]
            for i in xrange(len(totalsamples)):
                fpOffsetExp = map(lambda x: sum(dOffsetDic[pre][arm][x][i]), range(11))
                for fpOffset in range(11):
                    if not relativeMedianDic_FP.has_key(fpOffset): relativeMedianDic_FP[fpOffset] = []
                    relativeMedianDic_FP[fpOffset].append((fpOffsetExp[fpOffset]+0.1)/(fpOffsetExp[5]+0.1))
                    
            for tempFpOffset in range(len(dOffsetDic[pre][arm])):
                fpOffset = tempFpOffset - 5
                if not pre + '\t' + arm in corr_A: continue

                medianDic = dict()#; relativeMedianDic = dict()
                normalDic = dict(); cancerDic = dict()
                for i in xrange(len(totalsamples)):
                    sample = totalsamples[i]
                    for tpOffset in range(11):
                        if not medianDic.has_key(tpOffset): medianDic[tpOffset] = []
                        if not normalDic.has_key(tpOffset): normalDic[tpOffset] = []
                        if not cancerDic.has_key(tpOffset): cancerDic[tpOffset] = []
                        if len(dOffsetDic[pre][arm][int(fpOffset) + 5]) == 0:
                            medianDic[tpOffset].append(0)
                            if infoDic[sample] == 'Normal':
                                normalDic[tpOffset].append(0)
                            else:
                                cancerDic[tpOffset].append(0)

                        else:
                            medianDic[tpOffset].append(dOffsetDic[pre][arm][int(fpOffset) + 5][i][tpOffset])
                            if infoDic[sample] == 'Normal':
                                normalDic[tpOffset].append(dOffsetDic[pre][arm][int(fpOffset) + 5][i][tpOffset])
                            else:
                                cancerDic[tpOffset].append(dOffsetDic[pre][arm][int(fpOffset) + 5][i][tpOffset])
                medianSum = 0; normalSum = 0; cancerSum = 0 
                for tpOffset in medianDic.keys():
                    medianSum += np.mean(medianDic[tpOffset])
                    normalSum += np.mean(normalDic[tpOffset])
                    cancerSum += np.mean(cancerDic[tpOffset])

                if normalSum < 10 and cancerSum < 10: continue
                
                sampleFrequency = [0 for y in range(11)]
                relativeMedianDic = dict(); normalizedRelativeMedianDic = dict()
                for i in xrange(len(dOffsetDic[pre][arm][int(fpOffset) + 5])):
                    sample = totalsamples[i]
                    tpOffsetExp = dOffsetDic[pre][arm][int(fpOffset) + 5][i]
                    maxExp = max(tpOffsetExp)
                    maxNumber = tpOffsetExp.count(maxExp)
                    for tpOffset in range(11):
                        writer_DoubleOffset.write(pre + '\t' + arm + '\t' + sample + '\t' + str(fpOffset) + '\t' + str(int(tpOffset) - 5) + '\t' + str(tpOffsetExp[tpOffset]) + '\n')
                        if tpOffsetExp[tpOffset] == maxExp:
                            sampleFrequency[tpOffset] += 1.0 / maxNumber

                        if not relativeMedianDic.has_key(tpOffset): relativeMedianDic[tpOffset] = []
                        # add pseudocount=1
                        relativeMedianDic[tpOffset].append((tpOffsetExp[tpOffset]+0.1) / (float(tpOffsetExp[5])+0.1))
                        if not normalizedRelativeMedianDic.has_key(tpOffset): normalizedRelativeMedianDic[tpOffset] = []
                        normalizedRelativeMedianDic[tpOffset].append(((tpOffsetExp[tpOffset]+0.1) / (float(tpOffsetExp[5])+0.1)) / ((tpOffsetExp_norm[i][tpOffset]+0.1) / (float(tpOffsetExp_norm[i][5])+0.1)))
                writer_sampleFrequency.write(pre + '\t' + arm + '\t' + str(fpOffset) + '\t' + '\t'.join(map(str,sampleFrequency)) + '\n')

                for tpOffset in relativeMedianDic.keys():
                    median = np.mean(relativeMedianDic[tpOffset])
                    std = np.std(relativeMedianDic[tpOffset]) / math.sqrt(len(relativeMedianDic[tpOffset]))
                    writer_positionalExp.write(pre + '\t' + arm + '\t' + str(fpOffset) + '\t' + str(int(tpOffset) - 5) + '\t' + str(median) + '\t' + str(std) + '\n')

                    normalizedMedian = np.mean(normalizedRelativeMedianDic[tpOffset])
                    normalizedStd = np.std(normalizedRelativeMedianDic[tpOffset]) / math.sqrt(len(normalizedRelativeMedianDic[tpOffset]))
                    writer_positionalExp_norm.write(pre + '\t' + arm + '\t' + str(fpOffset) + '\t' + str(int(tpOffset) - 5) + '\t' + str(normalizedMedian) + '\t' + str(normalizedStd) + '\n')
            for fpOffset in relativeMedianDic_FP.keys():
                median = np.mean(relativeMedianDic_FP[fpOffset])
                std = np.std(relativeMedianDic_FP[fpOffset]) / math.sqrt(len(relativeMedianDic_FP[fpOffset]))
                writer_positionalExp_FP.write(pre + '\t' + arm + '\t' + str(fpOffset) + '\t' + str(int(fpOffset) - 5) + '\t' + str(median) + '\t' + str(std) + '\n')
    writer_sampleFrequency.close()
    writer_positionalExp.close()
    writer_positionalExp_FP.close()
    writer_positionalExp_norm.close()
    writer_DoubleOffset.close()


if __name__=='__main__':
    from rpy2.robjects.vectors import FloatVector
    from rpy2.robjects.packages import importr
    from scipy import stats
    import numpy as np
    import math 
    import sys
    import os
    import re
    cutoff = 10 # cutoff read counts ( sum of offset +1, offset -1, and canonical )
    isoformType = sys.argv[1] # first or second (major isomiR)
    hairpinF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hairpin_hsa.fa'
    matureF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/miRNA.fa'
    pairF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hsa.hg19.gff3'
    mirtronF = '/home/seokju/Project/LiverCancer/reference/human_mirtrons.txt'
    droshaDependencyF = '/home/seokju/Project/miRNAediting/3.ExpressionProfiling/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt'
    f = open(mirtronF); mirtrons = map(lambda x: x.strip(), f.readlines()); f.close()
    hairpinDic = get_hairpin(hairpinF)
    matureSeqDic = get_hairpin(matureF)
    pairDic, coordDic = get_pairDic(pairF)
    droshaDependentSet, dicerDependentSet = get_DroshaDependentSet(droshaDependencyF)
    armDic, armDic2 = get_AnnoArm(pairDic, matureSeqDic, hairpinDic, mirtrons, droshaDependentSet, dicerDependentSet)

    gradeType = "Total"
    catInfoF = '/home/seokju/Project/LiverCancer/reference/Catholic_info_v3.txt'
    catExpF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_neoplasm_Qnorm.txt'
    #catIsomiRD = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/1.miRDeep2/3.result_TotalSNV3/'
    catIsomiRD = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/3.result/'
    #catOutputD = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/' + gradeType + '/'
    catOutputD = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis/' + gradeType + '/' + isoformType + '/'

    tcgaInfoF = '/home/seokju/Project/LiverCancer/reference/TCGA_LIHC_Clinical_UCSC_v3_for_profiling.txt'
    #tcgaExpF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_neoplasm_Qnorm.txt'
    tcgaExpF = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm.txt'
    #tcgaIsomiRD = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/1.miRDeep2/4.result_TotalSNV_v3/'
    tcgaIsomiRD = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/3.result/'
    #tcgaOutputD = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/' + gradeType + '/'
    tcgaOutputD = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/5.IsomiR/MotifAnalysis/' + gradeType + '/' + isoformType + '/'

    tsinghwaInfoF = '/home/seokju/Project/LiverCancer/Tsinghua/info/Match_GSM_SRR_Sample.txt'
    tsinghwaExpF = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_All_Qnorm.txt'
    #tsinghwaIsomiRD = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/1.miRDeep2/3.result_TotalSNV2/'
    tsinghwaIsomiRD = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/3.result/'
    #tsinghwaOutputD = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/' + gradeType + '/'
    tsinghwaOutputD = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis/' + gradeType + '/'  + isoformType + '/'

    if not os.path.exists(catOutputD): os.makedirs(catOutputD)
    if not os.path.exists(tcgaOutputD): os.makedirs(tcgaOutputD)
    if not os.path.exists(tsinghwaOutputD): os.makedirs(tsinghwaOutputD)

    catInfoDic, catInfoDic2, catTotalsamples = get_InfoDic(catInfoF, gradeType)
    tcgaInfoDic, tcgaInfoDic2, tcgaTotalsamples = get_InfoDic(tcgaInfoF, gradeType)
    tsinghwaInfoDic, tsinghwaInfoDic2, tsinghwaTotalsamples = get_InfoDic(tsinghwaInfoF, gradeType)

    catExpDic, catTotalsamples2 = get_Exp(catExpF)
    tcgaExpDic, tcgaTotalsamples2 = get_Exp(tcgaExpF)
    tsinghwaExpDic, tsinghwaTotalsamples2 = get_Exp(tsinghwaExpF)

    if catTotalsamples != catTotalsamples2 or tcgaTotalsamples != tcgaTotalsamples2 or tsinghwaTotalsamples != tsinghwaTotalsamples2:
        print "Cat:", catTotalsamples != catTotalsamples2
        print "TCGA:", tcgaTotalsamples != tcgaTotalsamples2
        print "Tsingha:", tsinghwaTotalsamples != tsinghwaTotalsamples2
        print tsinghwaTotalsamples
        print tsinghwaTotalsamples2
        print "Different set between Exp and fold"
        exit()

    rcDic_cat, medianRcDic_cat, dOffsetDic_cat = get_ReadCount(catIsomiRD, catTotalsamples,armDic, catInfoDic2)
    rcDic_tcga, medianRcDic_tcga, dOffsetDic_tcga = get_ReadCount(tcgaIsomiRD, tcgaTotalsamples, armDic, tcgaInfoDic2)
    rcDic_tsinghwa, medianRcDic_tsinghwa, dOffsetDic_tsinghwa = get_ReadCount(tsinghwaIsomiRD, tsinghwaTotalsamples, armDic, tsinghwaInfoDic2)

    rcDic_cat, rcDic_tcga, rcDic_tsinghwa, medianRcDic_cat, medianRcDic_tcga, medianRcDic_tsinghwa, dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa, changedDic = adjust_offset(rcDic_cat, rcDic_tcga, rcDic_tsinghwa, medianRcDic_cat, medianRcDic_tcga, medianRcDic_tsinghwa, dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa)

    dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa, tpChangedDic = adjust_tpOffset(dOffsetDic_cat, dOffsetDic_tcga, dOffsetDic_tsinghwa, changedDic)

    catAfDic_for_compare_offset, catAfDic_for_compare_arm, catAfDic_for_compare_frequency, catAfDic_for_compare_frequency_up, catAfDic_for_compare_frequency_down, catAfDic_normal, catAfDic_cancer, catAfDic_all_sd, catAfDic_normal_sd, catAfDic_cancer_sd, cat_anovaDic, cat_ranksumDic = get_AlleleFrequency(medianRcDic_cat, rcDic_cat, cutoff, 'Catholic', isoformType)
    tcgaAfDic_for_compare_offset, tcgaAfDic_for_compare_arm, tcgaAfDic_for_compare_frequency, tcgaAfDic_for_compare_frequency_up, tcgaAfDic_for_compare_frequency_down, tcgaAfDic_normal, tcgaAfDic_cancer, tcgaAfDic_all_sd, tcgaAfDic_normal_sd, tcgaAfDic_cancer_sd, tcga_anovaDic, tcga_ranksumDic = get_AlleleFrequency(medianRcDic_tcga, rcDic_tcga, cutoff, 'TCGA', isoformType)
    tsinghwaAfDic_for_compare_offset, tsinghwaAfDic_for_compare_arm, tsinghwaAfDic_for_compare_frequency, tsinghwaAfDic_for_compare_frequency_up, tsinghwaAfDic_for_compare_frequency_down, tsinghwaAfDic_normal, tsinghwaAfDic_cancer, tsinghwaAfDic_all_sd, tsinghwaAfDic_normal_sd, tsinghwaAfDic_cancer_sd, tsinghwa_anovaDic, tsinghwa_ranksumDic= get_AlleleFrequency(medianRcDic_tsinghwa, rcDic_tsinghwa, cutoff, 'Tsinghua', isoformType)

    cat_corr = Calc_Correlation(medianRcDic_cat, rcDic_cat, catTotalsamples, catAfDic_for_compare_frequency, isoformType)
    tcga_corr = Calc_Correlation(medianRcDic_tcga, rcDic_tcga, tcgaTotalsamples, tcgaAfDic_for_compare_frequency, isoformType)
    tsinghwa_corr = Calc_Correlation(medianRcDic_tsinghwa, rcDic_tsinghwa, tsinghwaTotalsamples, tsinghwaAfDic_for_compare_frequency, isoformType)
    if not os.path.exists(catOutputD + '/Frequency/'): os.makedirs(catOutputD + '/Frequency/')
    if not os.path.exists(tcgaOutputD + '/Frequency/'): os.makedirs(tcgaOutputD + '/Frequency/')
    if not os.path.exists(tsinghwaOutputD + '/Frequency/'): os.makedirs(tsinghwaOutputD + '/Frequency/')

    intersectFrequencyDic = intersect_Frequency(catAfDic_for_compare_frequency, tcgaAfDic_for_compare_frequency, tsinghwaAfDic_for_compare_frequency, cat_corr, tcga_corr, tsinghwa_corr)

    cat_corr_N, cat_corr_C, cat_corr_A, cat_pDic= write_rcDic(catExpDic, medianRcDic_cat, rcDic_cat, catOutputD+ '/Frequency/', catTotalsamples, intersectFrequencyDic,'Frequency',catAfDic_for_compare_frequency, catInfoDic2, armDic2, changedDic, isoformType)
    tcga_corr_N, tcga_corr_C, tcga_corr_A, tcga_pDic = write_rcDic(tcgaExpDic, medianRcDic_tcga, rcDic_tcga, tcgaOutputD+ '/Frequency/', tcgaTotalsamples, intersectFrequencyDic,'Frequency',tcgaAfDic_for_compare_frequency, tcgaInfoDic2, armDic2, changedDic, isoformType)
    tsinghwa_corr_N, tsinghwa_corr_C, tsinghwa_corr_A, tsinghwa_pDic = write_rcDic(tsinghwaExpDic, medianRcDic_tsinghwa, rcDic_tsinghwa, tsinghwaOutputD + '/Frequency/', tsinghwaTotalsamples, intersectFrequencyDic,'Frequency',tsinghwaAfDic_for_compare_frequency, tsinghwaInfoDic2, armDic2, changedDic, isoformType)

    write_afDic(catAfDic_for_compare_frequency, catAfDic_for_compare_frequency_up, catAfDic_for_compare_frequency_down, tcgaAfDic_for_compare_frequency, tcgaAfDic_for_compare_frequency_up, tcgaAfDic_for_compare_frequency_down, tsinghwaAfDic_for_compare_frequency, tsinghwaAfDic_for_compare_frequency_up, tsinghwaAfDic_for_compare_frequency_down, catOutputD + '/Frequency/', 'Frequency', armDic2, matureSeqDic, hairpinDic, intersectFrequencyDic, changedDic, tpChangedDic, cat_corr_A, tcga_corr_A, tsinghwa_corr_A, medianRcDic_cat, medianRcDic_tcga, medianRcDic_tsinghwa, cat_pDic, tcga_pDic, tsinghwa_pDic, catAfDic_normal, catAfDic_cancer, catAfDic_all_sd, catAfDic_normal_sd, catAfDic_cancer_sd, tcgaAfDic_normal, tcgaAfDic_cancer, tcgaAfDic_all_sd, tcgaAfDic_normal_sd, tcgaAfDic_cancer_sd, tsinghwaAfDic_normal, tsinghwaAfDic_cancer, tsinghwaAfDic_all_sd, tsinghwaAfDic_normal_sd, tsinghwaAfDic_cancer_sd, cat_anovaDic, tcga_anovaDic, tsinghwa_anovaDic, cat_ranksumDic, tcga_ranksumDic, tsinghwa_ranksumDic)

    #write_dOffset(catAfDic_for_compare_frequency, cat_corr_A, dOffsetDic_cat, catTotalsamples, catInfoDic2, catOutputD + '/Frequency/')
    #write_dOffset(tcgaAfDic_for_compare_frequency, tcga_corr_A, dOffsetDic_tcga, tcgaTotalsamples, tcgaInfoDic2, tcgaOutputD + '/Frequency/')
    #write_dOffset(tsinghwaAfDic_for_compare_frequency, tsinghwa_corr_A, dOffsetDic_tsinghwa, tsinghwaTotalsamples, tsinghwaInfoDic2, tsinghwaOutputD + '/Frequency/')
