def get_expTable(expF, geneType, dataType):
    expDic = dict()
    f = open(expF); lines = f.readlines(); f.close()
    totalsamples = lines[0].strip().split('\t')[1:]
    index_83 = 0
    if 'Sample_83' in totalsamples:
        index_83 = totalsamples.index('Sample_83')
        totalsamples.pop(index_83)
    minimumCand = []
    for line in lines[1:]:
        line = line.strip().split('\t')
        id = line[0]
        exps = map(float,line[1:])
        if index_83:
            exps.pop(index_83)
        if geneType == 'miRNA' and max(exps) < 50: continue
        expDic[id] = exps

        if dataType == 'catholic':
            NormalExp = exps[:15]
        elif dataType == 'TCGA':
            NormalExp = exps[:50]
        else:
            NormalExp = exps[:20]
        minimumCand += filter(lambda x: x>0, NormalExp)
    return expDic, totalsamples, min(minimumCand)

def get_infoDic(infoF):
    infoDic = dict(); infoDic2 = dict()
    f = open(infoF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        grade = line[1]
        sample = line[2]
        if grade == 'Normal': newGrade = 'Normal'
        elif grade in ['G1','G2','G3','HPC','G4','Cancer','PVTT']: newGrade = 'Cancer'
        else: continue
        infoDic[sample] = newGrade
        if not infoDic2.has_key(newGrade): infoDic2[newGrade] = []
        infoDic2[newGrade].append(sample)
    
    return infoDic

def do_DEG(expDic, grade, direction, totalsamples, infoDic, NormalMin):
    DEG_P_Dic = dict(); foldchangeDic = dict()
    for id in expDic.keys():
        NormalExps = []
        for sample in totalsamples:
            if infoDic[sample] == 'Normal':
                NormalExps.append(expDic[id][totalsamples.index(sample)])
        GradeExps = []
        for sample in totalsamples:
            if infoDic[sample] == grade:
                GradeExps.append(expDic[id][totalsamples.index(sample)])
        
        NormalMedian = np.median(NormalExps)
        if NormalMedian == 0:
            NormalMedian = NormalMin
        foldchange = np.median(GradeExps) / NormalMedian

        foldchangeDic[id] = foldchange

        try:
            U, P = stats.mannwhitneyu(GradeExps, NormalExps, alternative=direction)
        except ValueError:
            print id, grade
            continue
        DEG_P_Dic[id] = P

    
    fdrDic = adjust_pvalues(DEG_P_Dic)
    return fdrDic, foldchangeDic
    
def get_DEG(expDic, totalsamples, infoDic, NormalMin):
    DEG_Dic = dict(); foldDic = dict()
    directionDic = dict()
    grades = list(set(infoDic.values()))
    directions = ['less','greater']

    threads = []
    for grade in grades:
        if grade == 'Normal': continue
        DEG_Dic[grade] = dict()
        for direction in directions:
            DEG_Dic[grade][direction] = dict()
            tempDEG_Dic, tempFoldchange_Dic = do_DEG(expDic, grade, direction, totalsamples, infoDic, NormalMin)
            DEG_Dic[grade][direction] = tempDEG_Dic
        foldDic[grade] = dict()
        foldDic[grade] = tempFoldchange_Dic
    for grade in grades:
        if grade == 'Normal': continue
        for direction in directions:
            for id in DEG_Dic[grade][direction].keys():
                if DEG_Dic[grade][direction][id] > 0.05: continue
                if not directionDic.has_key(id): directionDic[id] = dict()
                if not directionDic[id].has_key(grade): directionDic[id][grade] = direction
    return DEG_Dic, directionDic, foldDic

def intersect_DEGs(cat,tcga,tsinghwa):
    newCat = dict(); newTCGA = dict(); newTsinghwa = dict()
    for grade in cat.keys():
        if grade == 'Normal': continue
        for direction in ['less','greater']:
            for id in set(cat[grade][direction].keys()) & set(tcga[grade][direction].keys()) & set(tsinghwa[grade][direction].keys()):
                if not newCat.has_key(grade): newCat[grade] = dict()
                if not newTCGA.has_key(grade): newTCGA[grade] = dict()
                if not newTsinghwa.has_key(grade): newTsinghwa[grade] = dict()

                if not newCat[grade].has_key(direction): newCat[grade][direction] = dict()
                if not newTCGA[grade].has_key(direction): newTCGA[grade][direction] = dict()
                if not newTsinghwa[grade].has_key(direction): newTsinghwa[grade][direction] = dict()

                newCat[grade][direction][id] = cat[grade][direction][id]
                newTCGA[grade][direction][id] = tcga[grade][direction][id]
                newTsinghwa[grade][direction][id] = tsinghwa[grade][direction][id]
                if id in ["HNRNPC","U2AF2"]:
                    print grade, direction, id, cat[grade][direction][id], tcga[grade][direction][id], tsinghwa[grade][direction][id]
    return newCat, newTCGA, newTsinghwa

def adjust_pvalues(pvalueDic, method="BH"):
    stats = importr('stats')
    adjustDic = dict()
    temp_keys = []; temp_values = []
    for key in pvalueDic.keys():
        temp_keys.append(key)
        temp_values.append(pvalueDic[key])

    adjust_values = list(stats.p_adjust(FloatVector(temp_values),method=method))

    for i in xrange(len(temp_keys)):
        key = temp_keys.pop(0)
        FDR = adjust_values.pop(0)
        adjustDic[key] = FDR

    return adjustDic

def filter_further(miDEG_Dic, miDirectionDic, miFoldDic, miOriDEG_Dic, miOriDirectionDic):
    for grade in miDEG_Dic.keys():
        for direction in miDEG_Dic[grade].keys():
            for id in miDEG_Dic[grade][direction].keys():
                if not miOriDEG_Dic[grade][direction].has_key(id):
                    del miDEG_Dic[grade][direction][id]
                    continue
                if miOriDEG_Dic[grade][direction][id] > 0.05:
                    del miDEG_Dic[grade][direction][id]
                    if miDirectionDic.has_key(id):
                        if miDirectionDic[id].has_key(grade):
                            del miDirectionDic[id][grade]
    return miDEG_Dic, miDirectionDic, miFoldDic 

def get_matchDic(refFlatF):
    matchDic = dict(); matchDic2 = dict()
    f = open(refFlatF); lines = f.readlines(); f.close()
    for line in lines:
        line = line.strip().split('\t')
        geneSymbol = line[0]
        txnID = line[1]
        matchDic[txnID] = geneSymbol
        if not matchDic2.has_key(geneSymbol): matchDic2[geneSymbol] = []
        matchDic2[geneSymbol].append(txnID)
    return matchDic, matchDic2

def get_SameUTRIsoforms(refFlatF_ori, matchDic):
    # In TargetScan, they use unique isoform with unique 3' UTR start
    isoformDic = dict(); tempDic = dict()
    f = open(refFlatF_ori); lines = f.readlines(); f.close()
    for line in lines:
        line = line.strip().split('\t')
        geneSymbol = line[0]
        txnID = line[1]
        if txnID.startswith('ENSTR'): continue
        chr = line[2]
        strand = line[3]
        txStart = int(line[4])
        txEnd = int(line[5])
        cdStart = int(line[6])
        cdEnd = int(line[7])
        exStarts = map(float,line[9].split(',')[:-1])
        exEnds = map(float,line[10].split(',')[:-1])
        exons = []
        for i in xrange(len(exStarts)):
            exons.append([exStarts[i],exEnds[i]])
        fpUTRLength = 0
        if strand == '+':
            hit = 0
            for exon in exons:
                if hit == 0:
                    if exon[0] <= cdEnd <= exon[1]:
                        hit += 1
                        fpUTRLength += exon[1] - cdEnd
                    else:
                        continue
                else:
                    fpUTRLength += exon[1] - exon[0]
        else:
            hit = 0
            for exon in exons[::-1]:
                if hit == 0:
                    if exon[0] <= cdStart <= exon[1]:
                        hit += 1
                        fpUTRLength += cdStart - exon[0]
                    else:
                        continue
                else:
                    fpUTRLength += exon[1] - exon[0]

        if cdStart == cdEnd: continue
        
        if strand == '+':
            newName = geneSymbol + chr + strand + str(cdEnd)
        else:
            newName = geneSymbol + chr + strand + str(cdStart)

        if not tempDic.has_key(newName): tempDic[newName] = dict()
        tempDic[newName][txnID] = fpUTRLength

    for newName in tempDic.keys():
        txnIDs = tempDic[newName].keys()
        if len(set(map(lambda x: matchDic[x], txnIDs))) != 1:
            print txnIDs, "is from multi gene symbols !"
            exit()
        
        geneSymbol = map(lambda x: matchDic[x], txnIDs)[0]

        if not isoformDic.has_key(geneSymbol): isoformDic[geneSymbol] = []
        isoformDic[geneSymbol].append(txnIDs)
    return isoformDic

def get_MajorIsoform(matchDic2, catIsoformExpDic):
    majorDic = dict()
    for geneSymbol in matchDic2.keys():
        now = 0; major = ''
        for txnID in matchDic2[geneSymbol]:
            if not catIsoformExpDic.has_key(txnID): continue
            medianExp = np.median(catIsoformExpDic[txnID])
            if medianExp > now:
                now = medianExp
                major = txnID
            else:
                continue
        if major == '': continue
        majorDic[geneSymbol] = major#txnID
    return majorDic

def get_TargetRelationship(targetF, matchDic, isoformDic, majorDic):
    targetDic = dict()
    allTargetDic = dict()
    lineDic = dict()
    for line in open(targetF):
        if line.startswith('Gene'): continue
        line = line.strip().split('\t')
        repreTxnID = line[0]
        geneSymbol = matchDic[repreTxnID]
        miRNA = line[2]
        siteType = line[3]
        site = [line[4], line[5]]
        score = float(line[30])
        if score == 0: continue
        sameUTRIsoforms = filter(lambda x: repreTxnID in x, isoformDic[geneSymbol])
        if len(sameUTRIsoforms) != 1:
            print "Maybe multiple TxnIDs"
            print repreTxnID, geneSymbol, sameUTRIsoforms
            exit()
        for txnID in sameUTRIsoforms[0]:
            if majorDic[geneSymbol] != txnID: continue
            if not (score > -0.1 and siteType == '6mer'):
                if not targetDic.has_key(geneSymbol): targetDic[geneSymbol] = dict()
                if not targetDic[geneSymbol].has_key(miRNA): targetDic[geneSymbol][miRNA] = 0 
                targetDic[geneSymbol][miRNA] += score
            if not allTargetDic.has_key(geneSymbol): allTargetDic[geneSymbol] = dict()
            if not allTargetDic[geneSymbol].has_key(miRNA): allTargetDic[geneSymbol][miRNA] = 0 
            allTargetDic[geneSymbol][miRNA] += score
            if not lineDic.has_key(geneSymbol): lineDic[geneSymbol] = dict()
            if not lineDic[geneSymbol].has_key(miRNA): lineDic[geneSymbol][miRNA] = []
            lineDic[geneSymbol][miRNA].append(site)
    return targetDic, allTargetDic, lineDic

def get_GO(goF):
    goDic = dict()
    f = open(goF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        id = line[0]
        terms = []
        for term in line:
            terms += filter(lambda x: 'signaling pathway' in x.lower(), term.split(','))
        goDic[id] = terms
    return goDic

def get_Anno(gff3):
    annoDic = dict()
    f = open(gff3); lines = f.readlines(); f.close()
    pri = ''
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        maturationType = line[2]
        id = line[8].split(';')[2].split('=')[1]
        if line[2] == 'miRNA_primary_transcript':
            pri = id
        else:
            if not annoDic.has_key(pri): annoDic[pri] = []
            annoDic[pri].append(id)
    return annoDic

def get_representative(repreF):
    repreDic = dict()
    f = open(repreF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        repre = line[0]
        miRNAs = line[2].split('/')
        for miRNA in miRNAs:
            repreDic[miRNA] = repre
    return repreDic

def get_Known(knownF, annoDic, repreDic, matchDic2):
    knownGenes = []; knownMis = []
    f = open(knownF); lines = f.readlines(); f.close()
    for line in lines:
        if line.startswith('>'): continue
        line = line.strip().split(' ')
        miRNA = line[0].lower()
        geneSymbol = line[1]
        knownGenes.append(geneSymbol)
        if not annoDic.has_key(miRNA):
            print miRNA
            continue
        pairs = annoDic[miRNA]
        for pair in pairs:
            if not repreDic.has_key(pair): continue
            repre = repreDic[pair]
            knownMis.append(repre)
    return knownGenes, knownMis

def updateSurvival(survivalD, survivalDic, type, num):
    files = os.listdir(survivalD + '/summary/')
    for file in files:
        name = file.split('.txt')[0]
        if '_sig' in name:
            name = name.split('_sig')[0]
        f = open(survivalD + '/summary/'+file); lines = f.readlines(); f.close()
        expLine = filter(lambda x: x.startswith('"ExpUpper"'), lines)[0]
        coefficient = float(expLine.strip().split(' ')[1])
        pValue = float(expLine.strip().split(' ')[5])
        if (type == "multivariate" and len(lines) == 2) or pValue > 0.05:
            pValue =  1

        if not survivalDic.has_key(name): 
            survivalDic[name] = {'P':[1,1,1,1,1,1,1,1],'hazard':["","","","","","","",""]}
        survivalDic[name]['P'][num] = pValue

        if coefficient < 0:
            survivalDic[name]['hazard'][num] = 'Down'
        else:
            survivalDic[name]['hazard'][num] = 'Up'

def calc_Correlations(miExpDic, miNormExpDic, geneExpDic, geneNormExpDic, targetDic, corrDic, pValDic, geneSymbolList, results, i, miNormMin, geneNormMin):
    for geneSymbol in geneSymbolList:
        if not geneExpDic.has_key(geneSymbol): continue
        for miRNA in targetDic[geneSymbol].keys():
            if not miExpDic.has_key(miRNA): continue
            results[i] += 1
            corrR, corrP = stats.pearsonr(miNormExpDic[miRNA] + miExpDic[miRNA], geneNormExpDic[geneSymbol] + geneExpDic[geneSymbol])
            if np.median(miNormExpDic[miRNA]) == 0:
                miNorm = miNormMin
            else:
                miNorm = np.median(miNormExpDic[miRNA])
            miFoldChange = np.median(miExpDic[miRNA]) / miNorm

            if np.median(geneNormExpDic[geneSymbol]) == 0:
                geneNorm = geneNormMin
            else:
                geneNorm = np.median(geneNormExpDic[geneSymbol])
            geneFoldChange = np.median(geneExpDic[geneSymbol]) / geneNorm

            newName = geneSymbol + '\t' + miRNA
            if corrP <= 0.05 and ((miFoldChange > 1 and geneFoldChange < 1) or (miFoldChange < 1 and geneFoldChange > 1)):
                corrDic[newName] = corrR
                pValDic[newName] = corrP

def change_Grade(expDic, targetGrades, totalsamples, infoDic):
    newExpDic = dict()
    for id in expDic.keys():
        tempExps = []
        for i in xrange(len(totalsamples)):
            sample = totalsamples[i]
            if infoDic[sample] in targetGrades:
                tempExps.append(expDic[id][i])
        newExpDic[id] = tempExps
    return newExpDic
    

def get_Correlations(miExpDic, geneExpDic, miDirectionDic, geneDirectionDic, miFoldDic, geneFoldDic, totalsamples, infoDic, miDEG_Dic, geneDEG_Dic, targetDic, allTargetDic, thr, miNormalMin, geneNormalMin, survivalDic, goDic, outputD, knownGenes, knownMis, lineDic, matchDic):
    manager = mp.Manager()
   
    grades = list(set(infoDic.values()))
    grades.remove("Normal")
    if len(grades) > 1:
        grades = ['Fibrosis','Cirrhosis','Dysplastic','Cancer']
    else:
        grades = ['Cancer']
    geneSymbols = targetDic.keys()

    lists = []
    for i in xrange(thr):
        lists.append([])
    for i in xrange(len(geneSymbols)):
        geneSymbol = geneSymbols.pop(0)
        lists[i%thr].append(geneSymbol)

    corrDic = manager.dict(); pValDic = manager.dict(); fdrDic = dict()
    targetGrades = list(grades)
    miGradeExpDic = change_Grade(miExpDic, targetGrades, totalsamples, infoDic)
    miNormalExpDic = change_Grade(miExpDic, "Normal", totalsamples, infoDic)
    geneGradeExpDic = change_Grade(geneExpDic, targetGrades, totalsamples, infoDic)
    geneNormalExpDic = change_Grade(geneExpDic, "Normal", totalsamples, infoDic)
    results = manager.list([0]*thr)
    threads = []
    for i in xrange(len(lists)):
        geneSymbolList = lists[i]
        calcList = []
        threads.append(mp.Process(target=calc_Correlations,args=(miGradeExpDic,miNormalExpDic,geneGradeExpDic, geneNormalExpDic, targetDic, corrDic,pValDic, geneSymbolList, results, i, miNormalMin, geneNormalMin)))
        threads[-1].deamon=True
    
    for i in xrange(0,len(threads)):
        threads[i].start()

    for i in xrange(0,len(threads)):
        threads[i].join()

    totalPval = sum(results)
            
    fdrDic = adjust_Pvalue2(pValDic, totalPval)
    
    write_DEG(miDEG_Dic, geneDEG_Dic, outputD)
    for n in xrange(len(grades)):
        iterN = n + 1
        for tempGrades in itertools.combinations(grades,iterN):
            write_temp_result3(miDEG_Dic, geneDEG_Dic, miDirectionDic, geneDirectionDic, miFoldDic, geneFoldDic, fdrDic, corrDic, tempGrades, totalsamples, infoDic, targetDic, allTargetDic, survivalDic, goDic, outputD, knownGenes, knownMis, grades, lineDic, matchDic)
            
            write_finalTargetInfo(targetDic, miExpDic, geneExpDic, totalsamples, tempGrades, outputD)

def cmp1r(a1,a2): return cmp(a2[1],a1[1])

def adjust_Pvalue2(pDic, totalPval):
    adjustedPvalueDic = dict(); estimatedFDRDic = dict()
    allitems = pDic.items()
    allitems.sort(cmp1r)

    n = totalPval
    minimum = 1
    for item in allitems:
        key = item[0]
        pValue = item[1]
        adjustedPvalue = pValue * totalPval / n
        if adjustedPvalue < minimum:
            estimatedFDRDic[key] = adjustedPvalue
            minimum = adjustedPvalue
        else:
            estimatedFDRDic[key] = minimum
        n -= 1
    return estimatedFDRDic

def write_CorrGeneExp(corrfdrDic, geneExpDic, totalsamples, outputD):
    writer = open(outputD + 'Exp_CorrelatedGenes.txt','w')
    writer.write('GeneSymbol' + '\t' + '\t'.join(totalsamples) + '\n')
    totalNewName = []
    for grade in corrfdrDic.keys():
        totalNewName += corrfdrDic[grade]
    geneSymbols = set(map(lambda x: x.split('\t')[0], totalNewName))
    for geneSymbol in geneSymbols:
        writer.write(geneSymbol + '\t' + '\t'.join(map(str,geneExpDic[geneSymbol])) + '\n')
    writer.close()

def write_DEG(miDEG_Dic, geneDEG_Dic, outputD):
    for grade in miDEG_Dic.keys():
        writer_DEGene = open(outputD+'/DEGene_' + grade + '.txt','w')
        writer_DEMi = open(outputD+'/DEMi_' + grade + '.txt','w')
        writer_DEGene.write('GeneSymbol' + '\t' + 'DEGene_Direction' + '\n')
        writer_DEMi.write('miRNA' + '\t' + 'DEMi_Direction' + '\n')
        for direction in miDEG_Dic[grade].keys():
            for id in miDEG_Dic[grade][direction].keys():
                writer_DEMi.write(id + '\t' + direction + '\n')
            for id in geneDEG_Dic[grade][direction].keys():
                writer_DEGene.write(id + '\t' + direction + '\n')
        writer_DEGene.close()
        writer_DEMi.close()

def write_temp_result3(miDEG_Dic, geneDEG_Dic, miDirectionDic, geneDirectionDic, miFoldDic, geneFoldDic, fdrDic, corrDic, targetGrades, totalsamples, infoDic, targetDic, allTargetDic, survivalDic, goDic, outputD, knownGenes, knownMis, totalGrades, lineDic, matchDic):

    survivalList = ['Cat_RFS_Uni','Cat_RFS_Multi','Cat_OS_Uni','Cat_OS_Multi','TCGA_RFS_Uni','TCGA_RFS_Multi','TCGA_OS_Uni','TCGA_OS_Multi']

    writer_Interaction = open(outputD + '/TargetInteractions.txt','w')
    writer_correlation = open(outputD + '/Correlation_' + '_'.join(targetGrades) + '.txt','w')
    writer_DEGene = open(outputD + '/IntersectDEGene_' + '_'.join(targetGrades) + '.txt','w')
    writer_DEMi = open(outputD + '/IntersectDEMi_' + '_'.join(targetGrades) + '.txt','w')
    writer_DE = open(outputD + '/TargetDE_' + '_'.join(targetGrades) + '.txt','w')
    writer_surv = open(outputD + '/Survival_' + '_'.join(targetGrades) + '.txt','w')

    writer_go = open(outputD + '/GO_' + '_'.join(targetGrades) + '_CS-0.1.txt','w')
    writer_go_isomir = open(outputD + '/GO_' + '_'.join(targetGrades) + '_CS-0.1_IsomiR.txt','w')

    writer_correlation.write('GeneSymbol' + '\t' + 'txnID' + '\t' + 'DEGene_Direction' + '\t' + '\t'.join(map(lambda x: 'DEGene_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'GeneFoldchange_' + x, totalGrades)) + '\t' + 'miRNA' + '\t' + 'DEMi_Direction' + '\t' + '\t'.join(map(lambda x: 'DEMi_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'MiFoldchange_' + x, totalGrades)) + '\t' + 'Correlation' + '\t' + 'FDR' + '\t' + 'ContextScore' + '\t' + 'Survival' + '\t' + '\t'.join(map(lambda x: 'Gene_' + x, survivalList)) + '\t' + '\t'.join(map(lambda x: 'miRNA_' + x, survivalList)) + '\n')
    writer_DEGene.write('GeneSymbol' + '\t' + 'txnID' + '\t' + 'DEGene_Direction' + '\t' + '\t'.join(map(lambda x: 'DEG_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'Foldchange_' + x, totalGrades)) + '\n')
    writer_DEMi.write('miRNA' + '\t' + 'DEMi_Direction' + '\t' + '\t'.join(map(lambda x: 'DEMi_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'Foldchange_' + x, totalGrades)) + '\n')
    writer_DE.write('GeneSymbol' + '\t' + 'txnID' + '\t' + 'DEGene_Direction' + '\t' + '\t'.join(map(lambda x: 'DEGene_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'GeneFoldchange_' + x, totalGrades)) + '\t' + 'miRNA' + '\t' + 'DEMi_Direction' + '\t' + '\t'.join(map(lambda x: 'DEMi_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'MiFoldchange_' + x, totalGrades)) + '\t' + 'ContextScore' + '\n')
    writer_surv.write('GeneSymbol' + '\t' + 'txnID' + '\t' + 'DEGene_Direction' + '\t' + '\t'.join(map(lambda x: 'DEGene_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'GeneFoldchange_' + x, totalGrades)) + '\t' + 'miRNA' + '\t' + 'DEMi_Direction' + '\t' + '\t'.join(map(lambda x: 'DEMi_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'MiFoldchange_' + x, totalGrades)) + '\t' + 'ContextScore' + '\t' + 'Survival' + '\t' + '\t'.join(map(lambda x: 'Gene_' + x, survivalList)) + '\t' + '\t'.join(map(lambda x: 'miRNA_' + x, survivalList)) + '\t' + 'UTR_starts' + '\t' + 'UTR_ends' + '\t' + 'Targeted_by_Canonical' + '\n')
    writer_go.write('GeneSymbol' + '\t' + 'txnID' + '\t' + 'DEGene_Direction' + '\t' + '\t'.join(map(lambda x: 'DEGene_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'GeneFoldchange_' + x, totalGrades)) + '\t' + 'miRNA' + '\t' + 'DEMi_Direction' + '\t' + '\t'.join(map(lambda x: 'DEMi_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'MiFoldchange_' + x, totalGrades)) + '\t' + 'ContextScore' + '\t' + 'Survival' + '\t' + '\t'.join(map(lambda x: 'Gene_' + x, survivalList)) + '\t' + '\t'.join(map(lambda x: 'miRNA_' + x, survivalList))+ '\t' + 'GO' + '\t' + 'UTR_starts' + '\t' + 'UTR_ends' + '\n')
    writer_go_isomir.write('GeneSymbol' + '\t' + 'txnID' + '\t' + 'DEGene_Direction' + '\t' + '\t'.join(map(lambda x: 'DEGene_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'GeneFoldchange_' + x, totalGrades)) + '\t' + 'miRNA' + '\t' + 'DEMi_Direction' + '\t' + '\t'.join(map(lambda x: 'DEMi_pvalue_' + x, totalGrades)) + '\t' + '\t'.join(map(lambda x: 'MiFoldchange_' + x, totalGrades)) + '\t' + 'ContextScore' + '\t' + 'Survival' + '\t' + '\t'.join(map(lambda x: 'Gene_' + x, survivalList)) + '\t' + '\t'.join(map(lambda x: 'miRNA_' + x, survivalList))+ '\t' + 'GO' + '\t' + 'UTR_starts' + '\t' + 'UTR_ends' + '\n')

    for gene in targetDic.keys():
        for miRNA in targetDic[gene].keys():
            writer_Interaction.write(gene + '\t' + miRNA + '\t' + str(targetDic[gene][miRNA]) + '\n')
    gene_target_list = []
    for grade in targetGrades:
        tempList = []
        for direction in geneDEG_Dic[grade].keys():
            for id in geneDEG_Dic[grade][direction].keys():
                if geneDEG_Dic[grade][direction][id] > 0.05: continue
                tempList.append(id)
        if len(gene_target_list) == 0:
            gene_target_list += tempList
        else:
            gene_target_list = set(gene_target_list).intersection(set(tempList))
    for grade in totalGrades:
        if grade in targetGrades: continue
        tempList = []
        for direction in geneDEG_Dic[grade].keys():
            for id in geneDEG_Dic[grade][direction].keys():
                if geneDEG_Dic[grade][direction][id] > 0.05: continue
                tempList.append(id)
        gene_target_list = set(gene_target_list).difference(set(tempList))

    # write DEG result
    for geneSymbol in gene_target_list:
        direction = []; pvalue = []; foldchange = []
        for grade in totalGrades:
            foldchange.append(geneFoldDic[grade][geneSymbol])
            if geneDirectionDic[geneSymbol].has_key(grade):
                direction.append(geneDirectionDic[geneSymbol][grade])
                pvalue.append(geneDEG_Dic[grade][geneDirectionDic[geneSymbol][grade]][geneSymbol])
            else:
                direction.append('-')
                pvalue.append('NA')
        writer_DEGene.write(geneSymbol + '\t' + geneSymbol + '\t' + ','.join(direction) + '\t' + '\t'.join(map(str,pvalue)) + '\t' + '\t'.join(map(str,foldchange)) + '\n')
    print "Total DEGene :", len(gene_target_list)
    print "Known gene in DEGene :", set(gene_target_list) & set(knownGenes)

    target_list = []
    nontargets = []
    notingo = []
    for geneSymbol in gene_target_list:
        if not targetDic.has_key(geneSymbol):
            nontargets.append(geneSymbol)
            continue
        miRNAs = targetDic[geneSymbol]
        for miRNA in miRNAs:
            newName = geneSymbol + '\t' + miRNA
            if targetDic[geneSymbol][miRNA] > -0.1: continue

            test_DEG = []
            for grade in targetGrades:
                DEmi = 0
                for direction in miDEG_Dic[grade].keys():
                    if miDEG_Dic[grade][direction].has_key(miRNA):
                        DEmi = 1
                    else:
                        pass
                if DEmi == 0:
                    test_DEG.append(0)
                else:
                    if miDirectionDic[miRNA][grade] != geneDirectionDic[geneSymbol][grade]:
                        test_DEG.append(1)
                    else:
                        test_DEG.append(0)
            if test_DEG.count(0) > 0:
                continue
            else:
                target_list.append(newName)
    print "No targeted genes :",len(nontargets)
    mi_target_list = set(map(lambda x: x.split('\t')[1], target_list))

    # write DEMi
    for miRNA in mi_target_list:
        direction = []; pvalue = []; foldchange = []
        for grade in totalGrades:
            foldchange.append(miFoldDic[grade][miRNA])
            if miDirectionDic[miRNA].has_key(grade):
                direction.append(miDirectionDic[miRNA][grade])
                pvalue.append(miDEG_Dic[grade][miDirectionDic[miRNA][grade]][miRNA])
            else:
                direction.append('_')
                pvalue.append('NA')
        writer_DEMi.write(miRNA + '\t' + ','.join(direction) + '\t' + '\t'.join(map(str,pvalue)) + '\t' + '\t'.join(map(str,foldchange)) + '\n')
    print "Total DEMi :", len(mi_target_list)
    print "Known miRNA in DEMi :", set(mi_target_list) & set(knownMis)

    # write DEG + correlation result
    corrGene = []; corrMi = []
    survivalGene = []; survivalMi = []
    for newName in target_list:
        geneSymbol, miRNA = newName.split('\t')

        geneDirection = []; miDirection = []
        genePvalue = []; miPvalue = []
        geneFold = []; miFold = []
        for grade in totalGrades:
            geneFold.append(geneFoldDic[grade][geneSymbol])
            if geneDirectionDic[geneSymbol].has_key(grade):
                geneDirection.append(geneDirectionDic[geneSymbol][grade])
                genePvalue.append(geneDEG_Dic[grade][geneDirectionDic[geneSymbol][grade]][geneSymbol])
            else:
                geneDirection.append('-')
                genePvalue.append('NA')
            miFold.append(miFoldDic[grade][miRNA])
            if miDirectionDic[miRNA].has_key(grade):
                miDirection.append(miDirectionDic[miRNA][grade])
                miPvalue.append(miDEG_Dic[grade][miDirectionDic[miRNA][grade]][miRNA])
            else:
                miDirection.append('-')
                miPvalue.append('NA')
        writer_DE.write(geneSymbol + '\t' + geneSymbol + '\t' + ','.join(geneDirection) + '\t' + '\t'.join(map(str,genePvalue)) + '\t' + '\t'.join(map(str,geneFold)) + '\t' + miRNA + '\t' + ','.join(miDirection) + '\t' + '\t'.join(map(str,miPvalue)) + '\t' + '\t'.join(map(str,miFold)) + '\t' + str(targetDic[geneSymbol][miRNA]) + '\n')

        geneSurvPs = survivalDic[geneSymbol]['P']
        geneHazards = survivalDic[geneSymbol]['hazard']
        miSurvPs = survivalDic[miRNA]['P']
        miHazards = survivalDic[miRNA]['hazard']
        survWritten = 0
        if 'Cancer' in targetGrades:
            survPass = []
            miSurvPvalue = []; geneSurvPvalue = []
            for i in xrange(len(geneSurvPs)):
                geneSurvP = geneSurvPs[i]
                miSurvP = miSurvPs[i]
                geneHazard = geneHazards[i]
                miHazard = miHazards[i]
                if geneHazard == "" or miHazard == "": continue
                if geneSurvP <= 0.05 and miSurvP <= 0.05 and geneHazard != miHazard:
                    survPass.append(survivalList[i])
                if geneSurvP <= 0.05:
                    geneSurvPvalue.append(geneSurvP)
                else:
                    geneSurvPvalue.append('NA')
                if miSurvP <= 0.05:
                    miSurvPvalue.append(miSurvP)
                else:
                    miSurvPvalue.append('NA')
            if len(survPass) == 0: continue
            survivalGene.append(geneSymbol)
            survivalMi.append(miRNA)
            if miRNA.count('_') == 0:
                targetedByCano = 0
            else:
                canoMi = miRNA.split('_')[0]
                if allTargetDic[geneSymbol].has_key(canoMi):
                    targetedByCano = allTargetDic[geneSymbol][canoMi]
                else:
                    targetedByCano = 0
            writer_surv.write(geneSymbol + '\t' + geneSymbol + '\t' + ','.join(geneDirection) + '\t' + '\t'.join(map(str,genePvalue)) + '\t' + '\t'.join(map(str,geneFold)) + '\t' + miRNA + '\t' + ','.join(miDirection) + '\t' + '\t'.join(map(str,miPvalue)) + '\t' + '\t'.join(map(str,miFold)) + '\t' + str(targetDic[geneSymbol][miRNA]) + '\t' + ','.join(survPass) + '\t' + '\t'.join(map(str,geneSurvPvalue)) + '\t' + '\t'.join(map(str,miSurvPvalue)) + '\t' + ','.join(map(lambda x: x[0], lineDic[geneSymbol][miRNA])) + ',\t' + ','.join(map(lambda x: x[1], lineDic[geneSymbol][miRNA])) + ',' + '\t' + str(targetedByCano) + '\n')
            survWritten = 1 
            if fdrDic.has_key(newName):
                if fdrDic[newName] <= 0.05 and corrDic[newName] < 0:
                    corrGene.append(geneSymbol)
                    corrMi.append(miRNA)
                    writer_correlation.write(geneSymbol + '\t' + geneSymbol + '\t' + ','.join(geneDirection) + '\t' + '\t'.join(map(str,genePvalue)) + '\t' + '\t'.join(map(str,geneFold)) + '\t' + miRNA + '\t' + ','.join(miDirection) + '\t' + '\t'.join(map(str,miPvalue)) + '\t' + '\t'.join(map(str,miFold)) + '\t' + str(round(corrDic[newName],4)) + '\t' + str(round(fdrDic[newName],4)) + '\t' + str(targetDic[geneSymbol][miRNA]) + '\t' + ','.join(survPass) + '\t' + '\t'.join(map(str,geneSurvPvalue)) + '\t' + '\t'.join(map(str,miSurvPvalue)) + '\n')

            if not goDic.has_key(geneSymbol): continue
            if targetDic[geneSymbol][miRNA] > -0.1: continue
            writer_go.write(geneSymbol + '\t' + geneSymbol + '\t' + ','.join(geneDirection) + '\t' + '\t'.join(map(str,genePvalue)) + '\t' + '\t'.join(map(str,geneFold)) + '\t' + miRNA + '\t' + ','.join(miDirection) + '\t' + '\t'.join(map(str,miPvalue)) + '\t' + '\t'.join(map(str,miFold)) + '\t' + str(targetDic[geneSymbol][miRNA]) + '\t' + ','.join(survPass) + '\t' + '\t'.join(map(str,geneSurvPvalue)) + '\t' + '\t'.join(map(str,miSurvPvalue)) + '\t' + ','.join(goDic[geneSymbol]) + '\t' + ','.join(map(lambda x: x[0], lineDic[geneSymbol][miRNA])) + ',\t' + ','.join(map(lambda x: x[1], lineDic[geneSymbol][miRNA])) + ',\n')
            if miRNA.count('_') == 2:
                writer_go_isomir.write(geneSymbol + '\t' + geneSymbol + '\t' + ','.join(geneDirection) + '\t' + '\t'.join(map(str,genePvalue)) + '\t' + '\t'.join(map(str,geneFold)) + '\t' + miRNA + '\t' + ','.join(miDirection) + '\t' + '\t'.join(map(str,miPvalue)) + '\t' + '\t'.join(map(str,miFold)) + '\t' + str(targetDic[geneSymbol][miRNA]) + '\t' + ','.join(survPass) + '\t' + '\t'.join(map(str,geneSurvPvalue)) + '\t' + '\t'.join(map(str,miSurvPvalue)) + '\t' + ','.join(goDic[geneSymbol]) + '\t' + ','.join(map(lambda x: x[0], lineDic[geneSymbol][miRNA])) + ',\t' + ','.join(map(lambda x: x[1], lineDic[geneSymbol][miRNA])) + ',\n')

        else: pass

        if not goDic.has_key(geneSymbol) or survWritten: 
            continue

    print "Known genes in correlated Genes :", set(corrGene)&set(knownGenes)
    print "Known miRNAs in correlated miRNAs :", set(corrMi)&set(knownMis)
    print "Known genes in survival isoforms :", set(survivalGene)&set(knownGenes)
    print "Known miRNAs in survival miRNAs :", set(survivalMi)&set(knownMis)

    writer_Interaction.close()
    writer_correlation.close()
    writer_DEGene.close()
    writer_DEMi.close()
    writer_DE.close()
    writer_surv.close()
    writer_go.close()
    writer_go_isomir.close()
    
def write_finalTargetInfo(targetDic, miExpDic, geneExpDic, totalsamples, targetGrades, outputD):
    refF = outputD + '/Survival_' + '_'.join(targetGrades) + '.txt'
    tempDic = dict(); allIso = []; allCano = []; matchedCano = []
    f = open(refF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        geneSymbol = line[1]
        mi = line[5]
        if mi.count('_') > 0:
            allIso.append(mi)
            matchedCano.append(mi.split('_')[0])
        else:
            allCano.append(mi)
        if not tempDic.has_key(geneSymbol): tempDic[geneSymbol] = []
        tempDic[geneSymbol].append(mi)
    matchedIso = filter(lambda x: x.split('_')[0] in allCano, allIso)
    print matchedCano
    print matchedIso
    writer_info = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_info.txt','w')
    matchedCanoOnly = 0; matchedIsoOnly = 0; matchedBoth = 0
    canoOnly = 0; isoOnly = 0; both = 0
    for geneSymbol in tempDic.keys():
        cano = filter(lambda x: not '_' in x, tempDic[geneSymbol])
        iso = filter(lambda x: '_' in x, tempDic[geneSymbol])

        if len(cano) == 0:
            isoOnly += 1
        elif len(iso) == 0:
            canoOnly += 1
        else:
            both += 1
        
        mCano = set(tempDic[geneSymbol]) & set(matchedCano)
        mIso = set(tempDic[geneSymbol]) & set(matchedIso)
        if len(mCano) == 0 and len(mIso) == 0: pass
        elif len(mCano) == 0:
            matchedIsoOnly += 1
        elif len(mIso) == 0:
            matchedCanoOnly += 1
        else:
            matchedBoth += 1
    writer_info.write('Only targeted by canonical miRNA :' + '\t' + str(canoOnly) + '\n')
    writer_info.write('\t'+'matched canonical miRNA:' + '\t' + str(matchedCanoOnly) + '\n')
    writer_info.write('Only targeted by isomiR :' + '\t' + str(isoOnly) + '\n')
    writer_info.write('\t'+'matched isomiR:'+'\t'+str(matchedIsoOnly) + '\n')
    writer_info.write('Both :' + '\t' + str(both) + '\n')
    writer_info.write('\t'+'matched both:'+'\t'+str(matchedBoth) +'\n')
    writer_info.write('\n')
    writer_info.write('Enrichment of both ([[both,isomiR only],[cano only, others]],alternative="greater"): '+str(stats.fisher_exact([[both,isoOnly],[canoOnly,13143-both-isoOnly-canoOnly]],alternative="greater")[1])+'\n')
    writer_info.write('Enrichment of matched both: '+str(stats.fisher_exact([[matchedBoth,matchedIsoOnly],[matchedCanoOnly,13143-matchedBoth-matchedIsoOnly-matchedCanoOnly]],alternative="greater")[1])+'\n')
    writer_info.close()

    writer_gene = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_GeneExp.txt','w')
    writer_gene_cano = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_GeneExp_byCanonical.txt','w')
    writer_gene_matchedCano = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_GeneExp_byMatchedCanonical.txt','w')
    writer_gene_iso = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_GeneExp_byIsomiRs.txt','w')
    writer_gene_matchedIso = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_GeneExp_byMatchedIsomiRs.txt','w')
    writer_cwe = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_ContextscoreWeightedExp.txt','w')
    writer_mi = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_SummedMi.txt','w')
    writer_eachMi = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_Mi.txt','w')
    writer_corr = open(outputD + '/Survival_' + '_'.join(targetGrades) + '_Correlation.txt','w')

    header1 = 'Name' + '\t' + '\t'.join(totalsamples)
    header2 = 'Name' + '\t' + 'Correlation_SummedMi' + '\t' + 'Correlation_CWE'
    writer_gene.write(header1 + '\n')
    writer_gene_cano.write(header1 + '\n')
    writer_gene_matchedIso.write(header1 + '\n')
    writer_gene_matchedCano.write(header1 + '\n')
    writer_gene_iso.write(header1 + '\n')
    writer_cwe.write(header1 + '\n')
    writer_mi.write(header1 + '\n')
    writer_eachMi.write(header1 + '\n')
    writer_corr.write(header2 + '\n')
    hitMis = []
    for geneSymbol in tempDic.keys():
        writer_gene.write(geneSymbol + '\t' + '\t'.join(map(str, geneExpDic[geneSymbol])) + '\n')
        if len(filter(lambda x: x.count('_') == 0, tempDic[geneSymbol])) > 0:
            writer_gene_cano.write(geneSymbol + '\t' + '\t'.join(map(str, geneExpDic[geneSymbol])) + '\n')
        if len(filter(lambda x: x.count('_') > 0, tempDic[geneSymbol])) > 0:
            writer_gene_iso.write(geneSymbol + '\t' + '\t'.join(map(str,geneExpDic[geneSymbol])) + '\n')
        if len(set(tempDic[geneSymbol]) & set(matchedCano)) > 0:
            writer_gene_matchedCano.write(geneSymbol + '\t' + '\t'.join(map(str,geneExpDic[geneSymbol])) + '\n')
        if len(set(tempDic[geneSymbol]) & set(matchedIso)) > 0:
            writer_gene_matchedIso.write(geneSymbol + '\t' + '\t'.join(map(str,geneExpDic[geneSymbol])) + '\n')
        CWE = []
        SummedMi = []
        for miRNA in tempDic[geneSymbol]:
            if len(SummedMi) == 0:
                SummedMi += miExpDic[miRNA]
                CWE += map(lambda x: targetDic[geneSymbol][miRNA] * x, miExpDic[miRNA])
            else:
                for i in xrange(len(totalsamples)):
                    #print miRNA, i
                    SummedMi[i] += miExpDic[miRNA][i]
                    CWE[i] += targetDic[geneSymbol][miRNA] * miExpDic[miRNA][i]
            hitMis.append(miRNA)
        writer_cwe.write(geneSymbol + '\t' + '\t'.join(map(str,CWE)) + '\n')
        writer_mi.write(geneSymbol + '\t' + '\t'.join(map(str,SummedMi)) + '\n')
        writer_corr.write(geneSymbol + '\t' + str(stats.pearsonr(geneExpDic[geneSymbol], SummedMi)[0]) + '\t' + str(stats.pearsonr(geneExpDic[geneSymbol], CWE)[0]) + '\n')

    for miRNA in set(hitMis):
        writer_eachMi.write(miRNA + '\t' + '\t'.join(map(str,miExpDic[miRNA])) + '\n')
    writer_gene.close()
    writer_gene_cano.close()
    writer_gene_matchedCano.close()
    writer_gene_iso.close()
    writer_gene_matchedIso.close()
    writer_cwe.close()
    writer_mi.close()
    writer_eachMi.close()
    writer_corr.close()

if __name__=='__main__':
    from rpy2.robjects.vectors import FloatVector
    from rpy2.robjects.packages import importr
    import multiprocessing as mp
    from scipy import stats
    import numpy as np
    import itertools
    import sys
    import os

    thr = 50
    catMiOriF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
    catMiF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Family_7mer_Expression_neoplasm_Qnorm_RPMFiltered_Intersect.txt'
    catIsoformF = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Isoform_QuantileNormalized.txt' # Quatile normalized and for getting major isoform
    catGeneF = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Qnorm_RPKMFiltered_Intersect.txt'
    catInfoF = '/home/seokju/Project/LiverCancer/reference/Catholic_info_v3.txt'
    catOutputD = '/home/seokju/Project/LiverCancer/Catholic/3.TargetAnalysis/9.TargetAnalysis_noPVTT/TotalCorrelation_AllTarget_GeneLevel_include6mer_AllSurv/'
    tcgaMiOriF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
    tcgaMiF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Family_7mer_Expression_neoplasm_Qnorm_RPMFiltered_Intersect.txt'
    tcgaGeneF = '/home/seokju/Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_neoplasm_Qnorm_RPKMFiltered_Intersect.txt'
    tcgaInfoF = '/home/seokju/Project/LiverCancer/reference/TCGA_LIHC_Clinical_UCSC_v3_for_profiling.txt'
    tcgaOutputD = '/home/seokju/Project/LiverCancer/TCGA_rev1/3.TargetAnalysis/9.TargetAnalysis_noPVTT/TotalCorrelation_AllTarget_GeneLevel_include6mer_AllSurv/'

    tsinghwaMiOriF = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt'
    tsinghwaMiF = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Family_7mer_Expression_neoplasm_Qnorm_RPMFiltered_Intersect.txt'
    tsinghwaGeneF = '/home/seokju/Project/LiverCancer/Tsinghua/2.mRNA/3.ExpressionProfiling_noPVTT/Bitseq_ResultTotal_3PseqUpdatedFa_GeneLevelExp_All_Qnorm_RPKMFiltered_Intersect.txt'
    tsinghwaInfoF = '/home/seokju/Project/LiverCancer/Tsinghua/info/Match_GSM_SRR_Sample.txt'
    tsinghwaOutputD = '/home/seokju/Project/LiverCancer/Tsinghua/3.TargetAnalysis/9.TargetAnalysis_noPVTT/TotalCorrelation_AllTarget_GeneLevel_include6mer_AllSurv/'

    if not os.path.exists(catOutputD): os.makedirs(catOutputD)
    if not os.path.exists(tcgaOutputD): os.makedirs(tcgaOutputD)
    if not os.path.exists(tsinghwaOutputD): os.makedirs(tsinghwaOutputD)
    targetF = '/home/seokju/Project/LiverCancer/Catholic/3.TargetAnalysis/7.TargetScan_noPVTT/result/context_scores.txt'
    knownF = '/home/seokju/Project/LiverCancer/reference/miRNA_oncogenes_targeting_pairs_HCC.txt'
    repreF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/miRDeep2_Family_7mer_Infomation.txt'
    gff3 = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hsa.hg19.gff3'

    catMi_RFS_Uni_F = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS_uni/'
    catMi_RFS_Multi_F = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS/'
    catMi_OS_Uni_F = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS_uni/'
    catMi_OS_Multi_F = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS/'
    tcgaMi_RFS_Uni_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS_uni/'
    tcgaMi_RFS_Multi_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS/'
    tcgaMi_OS_Uni_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS_uni'
    tcgaMi_OS_Multi_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS/'

    catGene_RFS_Uni_F = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS_uni/'    
    catGene_RFS_Multi_F = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS/'    
    catGene_OS_Uni_F = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS_uni/'
    catGene_OS_Multi_F = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS/'
    tcgaGene_RFS_Uni_F ='/home/seokju/Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS_uni/'    
    tcgaGene_RFS_Multi_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/RFS/'
    tcgaGene_OS_Uni_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS_uni/'    
    tcgaGene_OS_Multi_F = '/home/seokju/Project/LiverCancer/TCGA_rev1/2.mRNA/3.ExpressionProfiling_noPVTT/SurvivalAnalysis/OS/'


    refFlatF = '/home/seokju/Data/gencode.v19.annotation.3PseqUpdated_rev1.refFlat'
    refFlatF_ori = '/home/seokju/Data/gencode.v19.annotation.refFlat'
    goF = '/home/seokju/Data/gencode.v19.annotation.3PseqUpdated_rev1_UniqSymbols_GO_SignalingPathway.txt'

    catMiOriExpDic, catTotalsamples, catMiNormalMin = get_expTable(catMiOriF, 'miRNA', 'catholic')
    tcgaMiOriExpDic, tcgaTotalsamples, tcgaMiNormalMin = get_expTable(tcgaMiOriF, 'miRNA', 'TCGA')
    tsinghwaMiOriExpDic, tsinghwaTotalsamples, tsinghwaMiNormalMin = get_expTable(tsinghwaMiOriF, 'miRNA', 'Tsinghwa')

    catMiExpDic, catTotalsamples, catMiNormalMin = get_expTable(catMiF, 'miRNA','catholic')
    tcgaMiExpDic, tcgaTotalsamples, tcgaMiNormalMin = get_expTable(tcgaMiF, 'miRNA','TCGA')
    tsinghwaMiExpDic, tsinghwaTotalsamples, tsinghwaMiNormalMin = get_expTable(tsinghwaMiF, 'miRNA','Tsinghwa')

    catGeneExpDic, catTotalsamples, catGeneNormalMin = get_expTable(catGeneF, 'mRNA','catholic')
    catIsoformExpDic, catsamples, catIsoformNormalMin = get_expTable(catIsoformF, 'mRNA','catholic')
    tcgaGeneExpDic, tcgaTotalsamples, tcgaGeneNormalMin = get_expTable(tcgaGeneF, 'mRNA','TCGA')
    tsinghwaGeneExpDic, tsinghwaTotalsamples, tsinghwaGeneNormalMin = get_expTable(tsinghwaGeneF, 'mRNA','Tsinghwa')

    catInfoDic = get_infoDic(catInfoF)
    tcgaInfoDic = get_infoDic(tcgaInfoF)
    tsinghwaInfoDic = get_infoDic(tsinghwaInfoF)

    matchDic, matchDic2 = get_matchDic(refFlatF)
    # isoform level
    isoformDic = get_SameUTRIsoforms(refFlatF_ori, matchDic)
    # gene level
    majorDic = get_MajorIsoform(matchDic2, catIsoformExpDic)
    targetDic, allTargetDic, lineDic = get_TargetRelationship(targetF, matchDic, isoformDic, majorDic)
    goDic = get_GO(goF)
    annoDic = get_Anno(gff3)
    repreDic = get_representative(repreF)
    knownGenes, knownMis = get_Known(knownF, annoDic, repreDic, matchDic2)

    catMiOriDEG_Dic, catMiOriDirectionDic, catMiOriFoldDic = get_DEG(catMiOriExpDic, catTotalsamples, catInfoDic, catMiNormalMin)
    catMiDEG_Dic, catMiDirectionDic, catMiFoldDic = get_DEG(catMiExpDic, catTotalsamples, catInfoDic, catMiNormalMin)
    catMiDEG_Dic, catMiDirectionDic, catMiFoldDic = filter_further(catMiDEG_Dic, catMiDirectionDic, catMiFoldDic, catMiOriDEG_Dic, catMiOriDirectionDic)
    catGeneDEG_Dic, catGeneDirectionDic, catGeneFoldDic = get_DEG(catGeneExpDic, catTotalsamples, catInfoDic, catGeneNormalMin)

    tcgaMiOriDEG_Dic, tcgaMiOriDirectionDic, tcgaMiOriFoldDic = get_DEG(tcgaMiOriExpDic, tcgaTotalsamples, tcgaInfoDic, tcgaMiNormalMin)
    tcgaMiDEG_Dic, tcgaMiDirectionDic, tcgaMiFoldDic = get_DEG(tcgaMiExpDic, tcgaTotalsamples, tcgaInfoDic, tcgaMiNormalMin)
    tcgaMiDEG_Dic, tcgaMiDirectionDic, tcgaMiFoldDic = filter_further(tcgaMiDEG_Dic, tcgaMiDirectionDic, tcgaMiFoldDic, tcgaMiOriDEG_Dic, tcgaMiOriDirectionDic)
    tcgaGeneDEG_Dic, tcgaGeneDirectionDic, tcgaGeneFoldDic = get_DEG(tcgaGeneExpDic, tcgaTotalsamples, tcgaInfoDic, tcgaGeneNormalMin)

    tsinghwaMiOriDEG_Dic, tsinghwaMiOriDirectionDic, tsinghwaMiOriFoldDic = get_DEG(tsinghwaMiOriExpDic, tsinghwaTotalsamples, tsinghwaInfoDic, tsinghwaMiNormalMin)
    tsinghwaMiDEG_Dic, tsinghwaMiDirectionDic, tsinghwaMiFoldDic = get_DEG(tsinghwaMiExpDic, tsinghwaTotalsamples, tsinghwaInfoDic, tsinghwaMiNormalMin)
    tsinghwaMiDEG_Dic, tsinghwaMiDirectionDic, tsinghwaMiFoldDic = filter_further(tsinghwaMiDEG_Dic, tsinghwaMiDirectionDic, tsinghwaMiFoldDic, tsinghwaMiOriDEG_Dic, tsinghwaMiOriDirectionDic)
    tsinghwaGeneDEG_Dic, tsinghwaGeneDirectionDic, tsinghwaGeneFoldDic = get_DEG(tsinghwaGeneExpDic, tsinghwaTotalsamples, tsinghwaInfoDic, tsinghwaGeneNormalMin)

    # Intersect DEGs
    catMiDEG_Dic, tcgaMiDEG_Dic, tsinghwaMiDEG_Dic = intersect_DEGs(catMiDEG_Dic, tcgaMiDEG_Dic, tsinghwaMiDEG_Dic)
    catGeneDEG_Dic, tcgaGeneDEG_Dic, tsinghwaGeneDEG_Dic = intersect_DEGs(catGeneDEG_Dic, tcgaGeneDEG_Dic, tsinghwaGeneDEG_Dic)
    exit()
    survivalDic = dict()
    updateSurvival(catMi_RFS_Uni_F, survivalDic, 'univariate',0)
    updateSurvival(catMi_RFS_Multi_F, survivalDic, 'multivariate',1)
    updateSurvival(catMi_OS_Uni_F, survivalDic, 'univariate',2)
    updateSurvival(catMi_OS_Multi_F, survivalDic, 'multivariate',3)
    updateSurvival(tcgaMi_RFS_Uni_F, survivalDic, 'univariate',4)
    updateSurvival(tcgaMi_RFS_Multi_F, survivalDic, 'multivariate',5)
    updateSurvival(tcgaMi_OS_Uni_F, survivalDic, 'univariate',6)
    updateSurvival(tcgaMi_OS_Multi_F, survivalDic, 'multivariate',7)

    updateSurvival(catGene_RFS_Uni_F, survivalDic, 'univariate',0)
    updateSurvival(catGene_RFS_Multi_F, survivalDic, 'multivariate',1)
    updateSurvival(catGene_OS_Uni_F, survivalDic, 'univariate',2)
    updateSurvival(catGene_OS_Multi_F, survivalDic, 'multivariate',3)
    updateSurvival(tcgaGene_RFS_Uni_F, survivalDic, 'univariate',4)
    updateSurvival(tcgaGene_RFS_Multi_F, survivalDic, 'multivariate',5)
    updateSurvival(tcgaGene_OS_Uni_F, survivalDic, 'univariate',6)
    updateSurvival(tcgaGene_OS_Multi_F, survivalDic, 'multivariate',7)
    get_Correlations(catMiExpDic, catGeneExpDic, catMiDirectionDic, catGeneDirectionDic, catMiFoldDic, catGeneFoldDic, catTotalsamples, catInfoDic, catMiDEG_Dic, catGeneDEG_Dic, targetDic, allTargetDic, thr, catMiNormalMin, catGeneNormalMin, survivalDic, goDic, catOutputD, knownGenes, knownMis, lineDic, matchDic)
    get_Correlations(tcgaMiExpDic, tcgaGeneExpDic, tcgaMiDirectionDic, tcgaGeneDirectionDic, tcgaMiFoldDic, tcgaGeneFoldDic, tcgaTotalsamples, tcgaInfoDic, tcgaMiDEG_Dic, tcgaGeneDEG_Dic, targetDic, allTargetDic, thr, tcgaMiNormalMin, tcgaGeneNormalMin, survivalDic, goDic, tcgaOutputD, knownGenes, knownMis, lineDic, matchDic)
    get_Correlations(tsinghwaMiExpDic, tsinghwaGeneExpDic, tsinghwaMiDirectionDic, tsinghwaGeneDirectionDic, tsinghwaMiFoldDic, tsinghwaGeneFoldDic, tsinghwaTotalsamples, tsinghwaInfoDic, tsinghwaMiDEG_Dic, tsinghwaGeneDEG_Dic, targetDic, allTargetDic, thr, tsinghwaMiNormalMin, tsinghwaGeneNormalMin, survivalDic, goDic, tsinghwaOutputD, knownGenes, knownMis, lineDic, matchDic)

