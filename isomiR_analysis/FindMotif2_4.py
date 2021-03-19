def cmp1(a1,a2): return cmp(a1[1],a2[1])

def get_Fold(foldF, heterogeneityDic):
    seqDic = dict(); foldDic = dict()
    f = open(foldF); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        lines = chunk.split('\n')
        id = lines[0].strip()
        if not heterogeneityDic.has_key(id): continue
        seq = lines[1].strip()
        fold = lines[2].strip().split(' ')[0]
        seqDic[id] = seq
        foldDic[id] = fold
    return seqDic, foldDic

def get_mGHGscore(ghgF):
    ghgDic = dict()
    f = open(ghgF); lines = f.readlines(); f.close()
    for line in lines[1:]:
        line = line.strip().split('\t')
        fp, tp, score = line 
        ghgDic[fp+tp] = float(score)
    return ghgDic

def get_Heterogeneity(heterogeneityF, dataType, offsetDirection):
    heterogeneityDic = dict(); frequencyDic = dict(); changedDic = dict()
    Items = []
    f = open(heterogeneityF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        id = line[0]
        if id in ['hsa-mir-3605','hsa-mir-215']: continue
        arm = line[1]
        if offsetDirection == "both":
            catHeterogeneity = line[2]
            tcgaHeterogeneity = line[3]
            tsinghwaHeterogeneity = line[4]
        elif offsetDirection == "minus":
            catHeterogeneity = line[9]
            tcgaHeterogeneity = line[11]
            tsinghwaHeterogeneity = line[13]
        elif offsetDirection == "plus":
            catHeterogeneity = line[10]
            tcgaHeterogeneity = line[12]
            tsinghwaHeterogeneity = line[14]

        frequencyTrend = line[6]
        cat_corr = line[14]
        tcga_corr = line[15]
        tsinghwa_corr = line[16]
        if arm == 'FP':
            if dataType == 'Catholic':
                if catHeterogeneity == 'NA': continue
                heterogeneityDic[id] = [float(catHeterogeneity)]
                Items.append([id,float(catHeterogeneity)])
            elif dataType == 'TCGA':
                if tcgaHeterogeneity == 'NA': continue
                heterogeneityDic[id] = [float(tcgaHeterogeneity)]
                Items.append([id,float(tcgaHeterogeneity)])
            elif dataType == 'Tsinghwa':
                if tsinghwaHeterogeneity == 'NA': continue
                heterogeneityDic[id] = [float(tsinghwaHeterogeneity)]
                Items.append([id,float(tsinghwaHeterogeneity)])
            elif dataType == 'All':
                #if tcgaHeterogeneity == 'NA' or tcgaHeterogeneity == 'NA' or tsinghwaHeterogeneity == 'NA': continue
                #heterogeneityDic[id] = [np.mean(map(float,[catHeterogeneity, tcgaHeterogeneity, tsinghwaHeterogeneity]))]
                #Items.append([id,np.mean(map(float,[catHeterogeneity, tcgaHeterogeneity, tsinghwaHeterogeneity]))])
                if catHeterogeneity == 'NA' or tsinghwaHeterogeneity == 'NA': continue
                #if cat_corr == 'F' or tsinghwa_corr == 'F': continue

                #if frequencyTrend == 'other': continue
                heterogeneityDic[id] = [np.mean(map(float,[catHeterogeneity, tsinghwaHeterogeneity]))]
                Items.append([id,np.mean(map(float,[catHeterogeneity, tsinghwaHeterogeneity]))])
                #frequencyDic[id] = frequencyTrend
            #catItems.append([id,catHeterogeneity])
        changed = int(line[7])

        #changedDic[id + '_' +arm] = changed
        if not changedDic.has_key(id): changedDic[id] = dict()
        changedDic[id][arm] = changed
        
    #print heterogeneityDic['hsa-mir-4454']
    #exit()
    #catItems.sort(cmp1)
    Items = heterogeneityDic.items()
    medianFrequency = np.median(map(lambda x: x[1], Items))
    #if dataType != "All":
        #medianFrequency = np.median(map(lambda x: x[1], Items))
    for item in Items:
        id, heterogeneity = item 
        if heterogeneity < medianFrequency:
            frequencyDic[id] = 'Low'
        else:
            frequencyDic[id] = 'High'
    #else:
    #    pass
    #for i in xrange(len(catItems)):
    #    id, heterogeneity = catItems[i]
    #    if i < len(catItems) / 2:
    #        frequencyDic[id] = 'Low'
    #    else:
    #        frequencyDic[id] = 'High'

    return heterogeneityDic, frequencyDic, medianFrequency, changedDic

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
    #pairDic['hsa-mir-1234'] = ['hsa-miR-1234-3p']
    pairDic['hsa-mir-939'] = ['hsa-miR-939-5p','hsa-miR-939-3p']
    return pairDic, coordDic

def get_Sequence(seqF):
    seqDic = dict()
    f = open(seqF); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        lines = chunk.split('\n')
        id = lines[0].split(' ')[0]
        seq = ''.join(lines[1:])
        seqDic[id] = seq
    return seqDic

def anno_CleavageSite(preSeqDic, matureSeqDic, structureDic, pairDic, coordDic, changedDic, cleavageBy):
    # annotation index of cleavage site on 30nt lengthend hairpin structure (mirbase v21)!
    indexDic = dict(); allIndexDic = dict(); cleavageIndexDic = dict()
    # allIndexDic : Drosha cleavage site in 5' and 3'
    # cleavageIndexDic : Drosha cleavage site in 5' and its complemently pairing site
    for pre in pairDic.keys():
        if not preSeqDic.has_key(pre): continue
        #if not changedDic.has_key(pre): continue
        preSeq = preSeqDic[pre]
        if pre == 'hsa-mir-1237': continue
        #spreSeq = spreSeqDic[pre]
        preStructure = structureDic[pre]
        pairs = pairDic[pre]
        chr, start, end, strand = coordDic[pre]
        if len(pairs) == 2:
            if len(filter(lambda x: x.endswith('-5p'), pairs)) + len(filter(lambda x: x.endswith('-3p'), pairs)) == 2:
                fpMi = filter(lambda x: x.endswith('-5p'), pairs)[0]
                tpMi = filter(lambda x: x.endswith('-3p'), pairs)[0]

                fpSeq = matureSeqDic[fpMi]
                if changedDic.has_key(pre): 
                    if changedDic[pre].has_key('FP'):
                        changedFP = changedDic[pre]['FP']
                    else:
                        changedFP = 0
                    if changedDic[pre].has_key('TP'):
                        changedTP = changedDic[pre]['TP']
                    else:
                        changedTP = 0
                else:
                    changedFP = 0
                    changedTP = 0
                #tempCompile = re.compile(fpSeq)
                fpIndexIter = re.finditer(r'(?=(%s))' % re.escape(fpSeq), preSeq)#tempCompile.finditer(preSeq)
                fpIndex = map(lambda x: x.span(), fpIndexIter)
                if len(fpIndex) > 1:
                    print pre
                    exit()
                fpStartIndex = preSeq.index(fpSeq) + changedFP
                FpUnpairNum = 0; FpPairNum = 0
                TpUnpairNum = 0; TpPairNum = 0
                hitTpIndex = 0
                for i in xrange(len(preSeq)):
                    if i < fpStartIndex: continue
                    else:
                        if preStructure[i] == '.' and FpPairNum == 0 and TpPairNum == 0:
                            FpUnpairNum += 1 
                        elif preStructure[i] == '(' and TpPairNum == 0:
                            FpPairNum += 1 
                        elif preStructure[i] == ')' and FpPairNum != 0:
                            TpPairNum += 1 
                        if TpPairNum == FpPairNum and FpPairNum != 0 and TpPairNum != 0:
                            hitTpIndex = i + 2 + FpUnpairNum
                            break
                        else:
                            pass

                if cleavageBy == 'Drosha':
                    #print pre
                    #print preSeq
                    #print matureSeqDic[fpMi], matureSeqDic[tpMi]
                    indexDic[pre] = [preSeq.index(matureSeqDic[fpMi])+changedFP, preSeq.index(matureSeqDic[tpMi])+changedTP + len(matureSeqDic[tpMi]) - 1]
                    #allIndexDic[pre] = [preSeq.index(matureSeqDic[fpMi])+changedFP, preSeq.index(matureSeqDic[tpMi])+changedTP + len(matureSeqDic[tpMi]) - 1]
                    allIndexDic[pre] = [fpStartIndex, hitTpIndex-2]
                    cleavageIndexDic[pre] = [preSeq.index(matureSeqDic[fpMi])+changedFP, hitTpIndex]

                elif cleavageBy == 'Dicer':
                    indexDic[pre] = [preSeq.index(matureSeqDic[fpMi]) + len(matureSeqDic[fpMi]) - 1, preSeq.index(matureSeqDic[tpMi])]
                else:
                    print "Drosha or Dicer !!"
                    exit()
            else:
                print pre,pairs
                exit()
        else:
            solo = pairs[0]
            f = 0; t = 0
            soloSeq = matureSeqDic[solo]
            print pre
            print preSeq
            print soloSeq
            if solo.endswith('-5p'): f = 1
            elif solo.endswith('-3p'): t = 1
            else:
                if preSeq.index(soloSeq) > len(preSeq) - (preSeq.index(soloSeq) + len(soloSeq)): t = 1
                else: f = 1
            if pre in ['hsa-mir-4487','hsa-mir-941-5','hsa-mir-941-2','hsa-mir-941-3','hsa-mir-941-4','hsa-mir-3142']:
                f = 0; t = 1
            if f == 1 and t == 0:
                if changedDic.has_key(pre):
                    if changedDic[pre].has_key('FP'):
                        changedFP = changedDic[pre]['FP']
                    else:
                        changedFP = 0
                    if changedDic[pre].has_key('TP'):
                        changedTP = changedDic[pre]['TP']
                    else:
                        changedTP = 0
                else:
                    changedFP = 0
                    changedTP = 0

                #tempCompile = re.compile(soloSeq)
                fpIndexIter = re.finditer(r'(?=(%s))' % re.escape(soloSeq), preSeq)#tempCompile.finditer(preSeq)
                fpIndex = map(lambda x: x.span(), fpIndexIter)
                if len(fpIndex) > 1 and not pre in ['hsa-mir-8061']: 
                    print pre
                    exit()
                if pre == 'hsa-mir-8061':
                    fpStartIndex = preSeq.rindex(soloSeq) + changedFP
                    #print fpStartIndex
                    #exit()
                else:
                    fpStartIndex = preSeq.index(soloSeq) + changedFP

                # should I predict non-annotated miRNA ?
                #fpStartIndex = preSeq.index(soloSeq)
                FpUnpairNum = 0; FpPairNum = 0
                TpUnpairNum = 0; TpPairNum = 0
                IndexFromDroshaCleavage = 0; FpUnpairInMatureNum = 0
                hitTpIndex = 0
                for i in xrange(len(preSeq)):
                    if i < fpStartIndex: continue
                    else:
                        if preStructure[i] == '.' and FpPairNum == 0 and TpPairNum == 0:
                            FpUnpairNum += 1
                        elif preStructure[i] == '(' and TpPairNum == 0:
                            FpPairNum += 1
                        elif preStructure[i] == ')' and FpPairNum != 0:
                            TpPairNum += 1
                        elif preStructure[i] == '(' and IndexFromDroshaCleavage > fpStartIndex + len(soloSeq) - 1 and TpPairNum == 0:
                            TpGapNum += 1
                        if TpPairNum == FpPairNum and FpPairNum != 0 and TpPairNum != 0:
                            hitTpIndex = i + 2 + FpUnpairNum
                            break
                        else:
                            pass
                        IndexFromDroshaCleavage += 1

                if cleavageBy == 'Drosha':
                    cleavageIndexDic[pre] = [fpStartIndex, hitTpIndex - 2]
                    allIndexDic[pre] = [fpStartIndex, hitTpIndex-2]
                    indexDic[pre] = [fpStartIndex, 'NA']
                elif cleavageBy == "Dicer":
                    allIndexDic[pre] = [fpStartIndex+len(soloSeq)-1, hitTpIndex-len(soloSeq)+1]
                    indexDic[pre] = [fpStartIndex+len(soloSeq)-1,'NA']
                else:
                    print "Drosha or Dicer !!!"
                    exit()
            elif f == 0 and t == 1:
                if changedDic.has_key(pre):
                    if changedDic[pre].has_key('FP'):
                        changedFP = changedDic[pre]['FP']
                    else:
                        changedFP = 0
                    if changedDic[pre].has_key('TP'):
                        changedTP = changedDic[pre]['TP']
                    else:
                        changedTP = 0
                else:
                    changedFP = 0
                    changedTP = 0

                #tempCompile = re.compile(soloSeq)
                fpIndexIter = re.finditer(r'(?=(%s))' % re.escape(soloSeq), preSeq)#tempCompile.finditer(preSeq)
                fpIndex = map(lambda x: x.span(), fpIndexIter)
                if len(fpIndex) > 1 and not pre in ['hsa-mir-8075','hsa-mir-4487','hsa-mir-941-5','hsa-mir-941-2','hsa-mir-941-3','hsa-mir-941-4','hsa-mir-3142']: 
                    print pre
                    exit()
                if pre == 'hsa-mir-8075':
                    tpEndIndex = preSeq.index(soloSeq)+changedTP + len(soloSeq) - 1
                elif pre == 'hsa-mir-3142':
                    tpEndIndex = fpIndex[1][0]+changedTP + len(soloSeq) - 1
                else:
                    tpEndIndex = preSeq.rindex(soloSeq)+changedTP + len(soloSeq) - 1

                TpUnpairNum = 0; TpPairNum = 0
                FpUnpairNum = 0; FpPairNum = 0
                hitFpIndex = 0
                for i in range(len(preSeq))[::-1]:
                    if tpEndIndex < i: continue
                    else:
                        if preStructure[i] == '.' and TpPairNum == 0 and FpPairNum == 0:
                            TpUnpairNum += 1
                        elif preStructure[i] == ')' and FpPairNum == 0:
                            TpPairNum += 1
                        elif preStructure[i] == '(' and TpPairNum != 0:
                            FpPairNum += 1
                        if FpPairNum == TpPairNum and FpPairNum != 0 and TpPairNum != 0:
                            hitFpIndex = i + 2 - TpUnpairNum
                            break
                        else:
                            pass

                if cleavageBy == 'Drosha':
                    cleavageIndexDic[pre] = [hitFpIndex, tpEndIndex]
                    allIndexDic[pre] = [hitFpIndex, tpEndIndex-2]
                    indexDic[pre] = ['NA', tpEndIndex]
                elif cleavageBy == 'Dicer':
                    allIndexDic[pre] = [hitFpIndex+len(soloSeq)-1,tpEndIndex-len(soloSeq)+1]
                    indexDic[pre] = ['NA', tpEndIndex-len(soloSeq)+1]
            else:
                print pre, solo
                exit()
        # find cleavage index ( 5' Dorhsa cleavage site and their pairing site)
        #f = ''; t = ''
        #for solo in pairs:
        #    soloSeq = matureSeqDic[solo]
        #    if solo.endswith('-5p'): f = solo
        #    elif solo.endswith('-3p'): t = solo
        #    else:
        #        if preSeq.index(soloSeq) > len(preSeq) - (preSeq.index(soloSeq) + len(soloSeq)): t = solo
        #        else: f = solo
        #if f 
    #print len(cleavageIndexDic.keys())
    #print cleavageIndexDic
    #exit()
    for pri in allIndexDic.keys():
        priSeq = preSeqDic[pri]
        preSeq = priSeq[allIndexDic[pri][0]:allIndexDic[pri][1]+1]
        pairs = pairDic[pri]
        for pair in pairs:
            matureSeq = matureSeqDic[pair]
            if preSeq.startswith(matureSeq) or preSeq.endswith(matureSeq):
                continue
            else:
                #print pri, pair
                #print priSeq
                #print matureSeq
                pass
    #print allIndexDic['hsa-mir-8075']
    #print allIndexDic['hsa-mir-8061']
    #exit()
    return indexDic, allIndexDic, cleavageIndexDic

def seperate_Sequence(priSeqDic, structureDic, allIndexDic, rev, outputD, heterogeneityDic, frequencyDic, analysisType):

    matchDic_ori = dict(); mismatchDic_ori = dict()
    leftBulgeDic_ori = dict(); rightBulgeDic_ori = dict()
    leftInternalDic_ori = dict(); rightInternalDic_ori = dict()
    leftInternalDic_Symm_ori = dict(); rightInternalDic_Symm_ori = dict()
    leftInternalDic_Assy_ori = dict(); rightInternalDic_Assy_ori = dict()
    myStructureDic = dict(); myStructureDic_37 = dict()
    earlyPassed = []

    matchDic = dict(); mismatchDic = dict()
    leftBulgeDic = dict(); rightBulgeDic = dict()
    leftInternalDic = dict(); rightInternalDic = dict()
    leftInternalDic_Symm = dict(); rightInternalDic_Symm = dict()
    leftInternalDic_Assy = dict(); rightInternalDic_Assy = dict()

    newPriSeqDic = dict(); newStructureDic = dict()
    matchDic_noB = dict()
    leftBulgeDic_noB = dict(); rightBulgeDic_noB = dict()
    leftInternalDic_noB = dict(); rightInternalDic_noB = dict()
    leftInternalDic_Symm_noB = dict(); rightInternalDic_Symm_noB = dict()
    leftInternalDic_Assy_noB = dict(); rightInternalDic_Assy_noB = dict()
    loopIndexDic_noB = dict(); newIndexDic_noB = dict()

    loopIndexDic = dict(); newIndexDic = dict()

    finalSequenceDic_37 = dict() # It will only have -13 ~ 23 index of no bulged sequence (no asymmetric internal loop)
                                 # , seperated left -13 ~ 23 from Drosha cleavage site and complement sequence
    finalStructureDic_37 = dict()
    finalSequenceDic_All = dict()
    finalStructureDic_All = dict()
    finalSequenceDic_Loop = dict() # -13 ~ loop
    overlap = []

    bulgeInCleavageSite = dict()
    assyInternalLoopInCleavageSite = []
    assyInternalLoopInCleavageSiteLargerInFP = []
    passed = []
    #trans = maketrans('()',')(')
    #loopCompile = re.compile(r"\(\.+\)")
    #mismatchCompile = re.compile(r"\.+")
    for id in allIndexDic.keys():
        if rev == False:
            preStart, preEnd = allIndexDic[id]
            priSeq = priSeqDic[id]
            priStructure = structureDic[id]
        else:
            priSeq = priSeqDic[id][::-1]
            priStructure = structureDic[id][::-1].translate(trans)
            tempPreStart, tempPreEnd = allIndexDic[id]
            preStart = len(priSeq) - 1 - tempPreEnd
            preEnd = len(priSeq) - 1 - tempPreStart
        newIndexDic[id] = [preStart,preEnd]
        loopIndexIter = loopCompile.finditer(priStructure)
        loopIndex = map(lambda x: x.span(), loopIndexIter)[0]
        #print loopIndex
        loopIndexDic[id] = loopIndex

        mismatchIndexIter = mismatchCompile.finditer(priStructure)
        mismatchIndex = map(lambda x: x.span(), mismatchIndexIter)
        mismatchDic[id] = map(lambda x: [x[0], x[1]-1], mismatchIndex)
        # Find all mismatches
        # and then, define internal loop or bulge
        # first, measure against terminal loop
        match, leftBulge, rightBulge, leftInternal, rightInternal, leftInternal_Symm, rightInternal_Symm, leftInternal_Assy, rightInternal_Assy = Find_mismatch(loopIndex, priStructure)
        #if id == 'hsa-mir-126':
        #    exit()
        '''
        if rev == True:
            hit = 0
            for tempIndexes in match + rightInternal:
                if preEnd in tempIndexes:
                    hit += 1
            if hit != 1:
                earlyPassed.append(id)
                continue
            #if not preEnd in 
        '''
        matchDic_ori[id] = match
        leftBulgeDic_ori[id] = leftBulge; rightBulgeDic_ori[id] = rightBulge
        leftInternalDic_ori[id] = leftInternal; rightInternalDic_ori[id] = rightInternal
        leftInternalDic_Symm_ori[id] = leftInternal_Symm; rightInternalDic_Symm_ori[id] = rightInternal_Symm
        leftInternalDic_Assy_ori[id] = leftInternal_Assy; rightInternalDic_Assy_ori[id] = rightInternal_Assy

        # Assign bulge, internal loop, mismatch to structure
        myStructure = list(priStructure)
        for tempMatch in match:
            for tempIndex in tempMatch:
                myStructure[tempIndex] = 'M'
        for tempBulge in leftBulge + rightBulge:
            myStructure[tempBulge] = 'B'
        for tempInternal in leftInternal + rightInternal:
            for tempIndex in tempInternal:
                myStructure[tempIndex] = 'I'
        for tempIndex in xrange(len(myStructure)):
            if myStructure[tempIndex] == '.':
                myStructure[tempIndex] = 'm'
        print id        
        print priSeq
        print priStructure
        print ''.join(myStructure)
        print match
        print leftBulge, rightBulge
        print leftInternal, rightInternal
        print leftInternal_Symm, rightInternal_Symm
        print leftInternal_Assy, rightInternal_Assy
        #exit()
        #if id == 'hsa-mir-3591': exit()
        if analysisType == 'Original':
            #oriSequence_FP = priSeq[preStart-13-10:preStart+25+1]
            oriSequence_FP = priSeq[preStart-30:preStart+30]
            #oriSequence_TP = priSeq[loopIndex[1] + loopIndex[0] - (preStart+25+1):loopIndex[1] + loopIndex[0] - (preStart+23+1)+37+10]
            #oriSequence_TP = priSeq[preEnd-25:preEnd+24]
            oriSequence_TP = priSeq[preEnd-29:preEnd+31]
            #oriStructure_FP = priStructure[preStart-13-10:preStart+25+1]
            oriStructure_FP = priStructure[preStart-30:preStart+30]
            #oriStructure_TP = priStructure[loopIndex[1] + loopIndex[0] - (preStart+25+1):loopIndex[1] + loopIndex[0] - (preStart+23+1)+37+10]
            #oriStructure_TP = priSeq[preEnd-25:preEnd+24]
            oriStructure_TP = priSeq[preEnd-29:preEnd+31]
            #print oriSequence_FP
            #print oriSequence_TP[::-1]
            #print oriStructure_FP
            #print oriStructure_TP[::-1]
            #exit()
            finalSequenceDic_37[id] = [oriSequence_FP, oriSequence_TP]
            finalStructureDic_37[id] = [oriStructure_FP, oriStructure_TP]

            finalSequenceDic_All[id] = priSeq
            finalStructureDic_All[id] = priStructure
            
            myStructureDic[id] = ''.join(myStructure)
            #myStructureDic_37[id] = [''.join(myStructure)[preStart-13-10:preStart+25+1],''.join(myStructure)[loopIndex[1] + loopIndex[0] - (preStart+25+1):loopIndex[1] + loopIndex[0] - (preStart+23+1)+37+10]]
            myStructureDic_37[id] = [''.join(myStructure)[preStart-13-10:preStart+25+1],''.join(myStructure)[preEnd-25:preEnd+24]]
            #print myStructureDic
            #exit()
            oriSequence_FP_Loop = priSeq[preStart-13:loopIndex[0]+1]
            oriSequence_TP_Loop = priSeq[loopIndex[1]-1:loopIndex[1]+len(oriSequence_FP_Loop)-1]
            finalSequenceDic_Loop[id] = [oriSequence_FP_Loop,oriSequence_TP_Loop]

            loopIndexDic_noB[id] = loopIndex
            newIndexDic_noB[id] = [preStart, preEnd]
            continue
        #tend1 = 0; tend2 = 0
        for i in xrange(len(leftInternal_Assy)):
            tempLeftAssyLoop = leftInternal_Assy[i]
            tempRightAssyLoop = rightInternal_Assy[i]
            print '\t left  :', map(lambda x: priSeq[x], tempLeftAssyLoop)[::-1], tempLeftAssyLoop[::-1]
            print '\t right :', map(lambda x: priSeq[x], tempRightAssyLoop)[::-1], tempRightAssyLoop[::-1]
            if len(tempLeftAssyLoop) > len(tempRightAssyLoop):
                leftBulge += tempLeftAssyLoop[::-1][-1*(len(tempLeftAssyLoop)-len(tempRightAssyLoop)):]
                tempIndexLeft = leftInternal.index(tempLeftAssyLoop)
                leftInternal.pop(tempIndexLeft)
                leftInternal.append(tempLeftAssyLoop[::-1][:-1*(len(tempLeftAssyLoop)-len(tempRightAssyLoop))][::-1])
                leftInternal_Symm.append(tempLeftAssyLoop[::-1][:-1*(len(tempLeftAssyLoop)-len(tempRightAssyLoop))][::-1])
                tempIndexRight = rightInternal.index(tempRightAssyLoop)
                rightInternal.pop(tempIndexRight)
                rightInternal.append(tempRightAssyLoop)
                rightInternal_Symm.append(tempRightAssyLoop)
            else:
                rightBulge += tempRightAssyLoop[::-1][-1*(len(tempRightAssyLoop)-len(tempLeftAssyLoop)):]
                tempIndexRight = rightInternal.index(tempRightAssyLoop)
                rightInternal.pop(tempIndexRight)
                rightInternal.append(tempRightAssyLoop[::-1][:-1*(len(tempRightAssyLoop)-len(tempLeftAssyLoop))][::-1])
                rightInternal_Symm.append(tempRightAssyLoop[::-1][:-1*(len(tempRightAssyLoop)-len(tempLeftAssyLoop))][::-1])
                tempIndexLeft = leftInternal.index(tempLeftAssyLoop)
                leftInternal.pop(tempIndexLeft)
                leftInternal.append(tempLeftAssyLoop)
                leftInternal_Symm.append(tempLeftAssyLoop)
            #continue
            # I decided to consider remained nucleotide whici is in assymetric internal loop and nearby terminal loop as mismatch.
            # So, assymetric internal loop is divided into mismatch and symmetric internal loop
            ######################################
            #if preStart in tempLeftAssyLoop:
            #    #tend1 += 1
            #    assyInternalLoopInCleavageSite.append([id,tempLeftAssyLoop,tempRightAssyLoop])
            #    if len(tempLeftAssyLoop) > len(tempRightAssyLoop):
            #        #tend2 += 1
            #        assyInternalLoopInCleavageSiteLargerInFP.append([id,tempLeftAssyLoop,tempRightAssyLoop])
            #        passed.append(id)
            #        continue
        #if tend1 != 0:
        #    assyInternalLoopInCleavageSite.append(id)
        #if tend2 != 0:
        #    assyInternalLoopInCleavageSiteLargerInFP.append(id)
        print "After remove assy internal loop"
        print leftBulge, rightBulge
        print leftInternal, rightInternal
        print leftInternal_Assy, rightInternal_Assy
        print leftInternal_Symm, rightInternal_Symm

        matchDic[id] = match
        leftBulgeDic[id] = leftBulge; rightBulgeDic[id] = rightBulge
        leftInternalDic[id] = leftInternal; rightInternalDic[id] = rightInternal
        leftInternalDic_Symm[id] = leftInternal_Symm; rightInternalDic_Symm[id] = rightInternal_Symm
        leftInternalDic_Assy[id] = leftInternal_Assy; rightInternalDic_Assy[id] = rightInternal_Assy
        preStart_noB = preStart; preEnd_noB = preEnd
        newSequence_noB = list(priSeq)
        newStructure_noB = list(priStructure)

        if preStart in leftBulge:
            #bulgeInCleavageSite.append(id)
            bulgeInCleavageSite[id] = [preStart, priSeq[preStart-1:preStart+2], leftBulge]
            passed.append(id)
            #print id
            #exit()
            continue

        for i in xrange(len(rightBulgeDic[id])):
            tempIndex = rightBulgeDic[id][::-1][i]
            if newStructure_noB[tempIndex] != '.':
                print id, tempIndex
                print priStructure
                exit()
            else:
                newStructure_noB[tempIndex] = ''
                newSequence_noB[tempIndex] = ''
                if tempIndex < preStart:
                    preStart_noB -= 1
                if tempIndex < preEnd_noB:
                    preEnd_noB -= 1
        for i in xrange(len(leftBulgeDic[id])):
            tempIndex = leftBulgeDic[id][i]
            if newStructure_noB[tempIndex] != '.':
                print id, tempIndex
                print priStructure
                exit()
            else:
                newStructure_noB[tempIndex] = ''
                newSequence_noB[tempIndex] = ''
                if tempIndex < preStart:
                    preStart_noB -= 1
                if tempIndex < preEnd:
                    preEnd_noB -= 1
        newSequence_noB = ''.join(newSequence_noB)
        newStructure_noB = ''.join(newStructure_noB)

        loopIndex_noB = map(lambda x: x.span(), loopCompile.finditer(newStructure_noB))[0]
        match_noB, leftBulge_noB, rightBulge_noB, leftInternal_noB, rightInternal_noB, leftInternal_Symm_noB, rightInternal_Symm_noB, leftInternal_Assy_noB, rightInternal_Assy_noB = Find_mismatch(loopIndex_noB, newStructure_noB)

        if len(leftInternal_Assy_noB) > 0 or len(rightInternal_Assy_noB) > 0:
            print "After remove bulges"
            print leftInternal_Assy_noB, rightInternal_Assy_noB
            exit()
        if leftInternal_noB != leftInternal_Symm_noB or rightInternal_noB != rightInternal_Symm_noB:
            print "No match internal and internal_Symm"
            exit()

        newPriSeqDic[id] = newSequence_noB
        newStructureDic[id] = newStructure_noB
        matchDic_noB[id] = match_noB
        leftBulgeDic_noB[id] = []; rightBulgeDic_noB[id] = []
        leftInternalDic_noB[id] = leftInternal_noB; rightInternalDic_noB[id] = rightInternal_noB
        leftInternalDic_Symm_noB[id] = leftInternal_Symm_noB; rightInternalDic_Symm_noB[id] = rightInternal_Symm_noB
        leftInternalDic_Assy_noB[id] = leftInternal_Assy_noB; rightInternalDic_Assy_noB[id] = rightInternal_Assy_noB
        #loopIndexDic_noB[id] = loopIndex_noB
        #newIndexDic_noB[id] = [preStart_noB, preEnd_noB]
        
        # After remove bulge
        print "new miRNA !!"
        print "Cleavage :", preStart_noB, preEnd_noB, loopIndex_noB
        print newSequence_noB
        print newStructure_noB
        print newSequence_noB[preStart_noB-13:preStart_noB+23+1]
        print newSequence_noB[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1):loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37][::-1]
        print newStructure_noB[preStart_noB-13:preStart_noB+23+1], newStructure_noB[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1):loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37]
        if loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1)*2 < 0:
            overlap.append(id)
        #if id == 'hsa-mir-411':
        #    exit()
        # Now, we can use just leftInternalDic and rightInternalDic, and matchDic

        newSequence_noB_FP = newSequence_noB[preStart_noB-13-10:preStart_noB+23+1]
        newSequence_noB_TP = newSequence_noB[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1):loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37+10]
        newStructure_noB_FP = newStructure_noB[preStart_noB-13-10:preStart_noB+23+1]
        newStructure_noB_TP = newStructure_noB[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1):loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37+10]

        newSequence_noB_FP_Loop = newSequence_noB[preStart_noB-13:loopIndex_noB[0] + 1]
        newSequence_noB_TP_Loop = newSequence_noB[loopIndex_noB[1]-1:loopIndex_noB[1] + len(newSequence_noB_FP_Loop)-1]
        if analysisType == 'noB':
            finalSequenceDic_37[id] = [newSequence_noB_FP,newSequence_noB_TP]
            finalStructureDic_37[id] = [newStructure_noB_FP, newStructure_noB_TP]

            finalSequenceDic_All[id] = newSequence_noB
            finalStructureDic_All[id] = newStructure_noB

            finalSequenceDic_Loop[id] = [newSequence_noB_FP_Loop,newSequence_noB_TP_Loop]
            loopIndexDic_noB[id] = loopIndex_noB
            newIndexDic_noB[id] = [preStart_noB, preEnd_noB]
        print "To Loop Position !!"
        print newSequence_noB_FP_Loop
        print newSequence_noB_TP_Loop
        
    print overlap
    print "Bulge in Drosha cleavage site:", bulgeInCleavageSite
    for id in bulgeInCleavageSite.keys():
        if heterogeneityDic.has_key(id):
            print id, heterogeneityDic[id], frequencyDic[id], bulgeInCleavageSite[id]
    print "Assymetric internal loop in Drosha cleavage site:", assyInternalLoopInCleavageSite
    print "Larger 5' arm assymetric internal loop than 3' arm in Drosha cleavage site:", assyInternalLoopInCleavageSiteLargerInFP
    #print len(set(heterogeneityDic.keys()) - set(passed))
    #exit()
    '''
    # Find basal junction
    if rev == True:
        pass
    else:
        junctionDic = dict()
        notAnnotate = []
        # Find from match dictionary ?
        # it would be better to find junction because mismatch dictionary is not consistent in indexing
        for id in matchDic.keys():
            if id in passed: continue
            #starts = map(lambda x: x[0], matchDic[id])
            
            # find match after mismatch
            #print len(matchDic[id])
            matchedSiteFromSS = matchDic[id][::-1]
            #print matchedSiteFromSS
            matchedSiteAfterMM = []
            wasStart = -1; wasEnd = -1; mismatch = 0
            for i in xrange(len(matchedSiteFromSS)):
                nowStart = matchedSiteFromSS[i][0]
                nowEnd = matchedSiteFromSS[i][1]
                if i == 0: # First matched site
                    matchedSiteAfterMM.append([nowStart,nowEnd])
                    wasStart = nowStart
                    wasEnd = nowEnd
                    continue
                if nowStart - wasStart == 1 and nowEnd - wasEnd == -1: # continuously matched site is passed
                    wasStart = nowStart
                    wasEnd = nowEnd
                    continue
                else: # others are appended to the list
                    wasStart = nowStart
                    wasEnd = nowEnd
                    matchedSiteAfterMM.append([nowStart,nowEnd]) 
    
            #nearJunction = filter(lambda x: 13 <= x[0] <= 21, matchDic[id])
            #tempMid = 17 # -13 index from Drosha cleavage site is same with 17 index from first nt
            
            #tend14 = 0; tend15 = 0; tend16 = 0; tend17 = 0; tend18 = 0; tend19 = 0; tend20 = 0
            #print '>'+id
            #print priSeqDic[id]
            #print structureDic[id]
            #if 17 in matchedSiteAfterMM:
                # correction by UG motif
            if id == 'hsa-mir-215':
                junctionDic[id] = filter(lambda x: x[0] == 19, matchedSiteAfterMM)[0]
            elif id == 'hsa-mir-99a':
                junctionDic[id] = filter(lambda x: x[0] == 14, matchedSiteAfterMM)[0]
            elif id == 'hsa-mir-204':
                junctionDic[id] = filter(lambda x: x[0] == 18, matchedSiteAfterMM)[0]
            elif id == 'hsa-mir-3613':
                junctionDic[id] = filter(lambda x: x[0] == 21, matchedSiteAfterMM)[0]
            elif id == 'hsa-mir-3591':
                junctionDic[id] = filter(lambda x: x[0] == 18, matchedSiteAfterMM)[0]
            elif id == 'hsa-mir-653':
                junctionDic[id] = filter(lambda x: x[0] == 19, matchedSiteAfterMM)[0]
            elif id == 'hsa-mir-618':
                junctionDic[id] = filter(lambda x:x [0] == 20, matchedSiteAfterMM)[0]
            else:
                if len(filter(lambda x: x[0] == 17, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 17, matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 16, matchedSiteAfterMM)) > 0 and len(filter(lambda x: x[0] == 18, matchedSiteAfterMM)) > 0:
                    print id, 16, "and", 18
                    print matchedSiteAfterMM
                elif len(filter(lambda x: x[0] == 16, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 16, matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 18, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 18, matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 15, matchedSiteAfterMM)) > 0 and len(filter(lambda x: x[0] == 19, matchedSiteAfterMM)) > 0:
                    print id, 15, "and", 19
                    print matchedSiteAfterMM
                elif len(filter(lambda x: x[0] == 15, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 15, matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 19, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 19 ,matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 14, matchedSiteAfterMM)) > 0 and len(filter(lambda x: x[0] == 20, matchedSiteAfterMM)) > 0:
                    print id, 14, "and", 20
                    print matchedSiteAfterMM
                elif len(filter(lambda x: x[0] == 14, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 14 ,matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 20, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 20, matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 13, matchedSiteAfterMM)) > 0 and len(filter(lambda x: x[0] == 21, matchedSiteAfterMM)) > 0:
                    print id, 13, "and", 21
                    print matchedSiteAfterMM
                elif len(filter(lambda x: x[0] == 13, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 13 ,matchedSiteAfterMM)[0]
                elif len(filter(lambda x: x[0] == 21, matchedSiteAfterMM)) > 0:
                    junctionDic[id] = filter(lambda x: x[0] == 21, matchedSiteAfterMM)[0]
                else:
                    print id, "No in 13~21"
                    print matchedSiteAfterMM
                    if id == 'hsa-mir-4662a':
                        junctionDic[id] = filter(lambda x: x[0] == 12, matchedSiteAfterMM)[0]
                        print "Junction of", id, ": 12"
                    else:
                        notAnnotate.append(id)
        # write stem length
        # stem length will be sum of the matched, bulged(both arm), symmetric internal loop, asymmetric internal loop(big one)
        #print matchDic
        #print junctionDic
        #exit()
        for id in allIndexDic.keys():
            if id in passed: continue
            if id in notAnnotate: continue
            #matchedPosition = filter(lambda x: junctionDic[id][0] <= x[0] and x[1] <= junctionDic[id][1], matchDic[id])
            #bulgeSite = filter(lambda x: junctionDic[id][0] <= x <= junctionDic[id][1], leftBulgeDic[id]) + filter(lambda x: junctionDic[id][0] <= x <= junctionDic[id][1], rightBulgeDic[id])
            #internalLoop_asym = filter(lambda x: 
            #print matchDic[id]
            #print junctionDic[id]
            #print matchedPosition
        #print leftInternalDic_Assy
        #exit()
        # write junction distribution
        #if rev == True:
        #    pass
        #else:
        writer = open(outputD + '/BasalJunctionDistribution.txt','w')
        writer_junctionUG = open(outputD + '/BasalJunctionUG_Distribution.txt','w')
        header = 'Index' + '\t' + 'Hit' + '\n'
        writer.write(header)
        writer_junctionUG.write('Index' + '\t' + 'UG1' + '\t' + 'UG2' + '\n')
        tempDic1 = dict(); tempDic2 = dict(); tempDic3 = dict()#temp = 0
        for id in junctionDic.keys():
            if id in passed: continue
            if id in notAnnotate: continue
            #junction = allIndexDic[id][0] - junctionDic[id][0]
            junction = junctionDic[id][0] - allIndexDic[id][0]
            if junction == -13:
                print id#, leftBulgeDic[id], rightBulgeDic[id]
                print priSeqDic[id]
                print structureDic[id][:allIndexDic[id][0]], "5p arm"
                print structureDic[id][:allIndexDic[id][1]+1], "3p arm"
                for tempIndex in leftBulgeDic[id]:
                    print structureDic[id][:tempIndex+1], "5p bulge"
                for tempIndex in rightBulgeDic[id]:
                    print structureDic[id][:tempIndex+1], "3p bulge"
                print structureDic[id][:junctionDic[id][0]+1],"5p junction"
                print structureDic[id][:junctionDic[id][1]+1],"3p junction"
                print structureDic[id]
            #priStructure = structureDic[id]
            priSeq = priSeqDic[id]
            # find UG at annotated junction
            if priSeq[junctionDic[id][0]-1:junctionDic[id][0]+1] == 'UG':
                junctionUG = 1
            else:
                junctionUG = 0
            # find all UG at each position
            #for tempIndex in map(lambda x: x[0], junctionDic.values()):
            for tempInd in range(12,22):
                tempIndex = tempInd - allIndexDic[id][0]
                if not tempDic3.has_key(tempIndex): tempDic3[tempIndex] = 0
                if priSeq[tempInd-1:tempInd+1] == 'UG':
                    tempDic3[tempIndex] += 1
                else:
                    continue
            #if priSeq[16:18] == 'UG':
            #    thirteenUG = 1
            #else:
            #    thirteenUG = 0
            if not tempDic1.has_key(junction): tempDic1[junction] = 0
            tempDic1[junction] += 1
            if not tempDic2.has_key(junction): tempDic2[junction] = 0
            tempDic2[junction] += junctionUG
            #if not tempDic3.has_key(junction): tempDic3[junction] = 0
            #temp += thirteenUG
        #for junction in tempDic1.keys():
        for junction in range(min(tempDic1.keys()), max(tempDic1.keys())+1):
            if not tempDic1.has_key(junction):
                temp1 = 0
            else:
                temp1 = tempDic1[junction]
            if not tempDic2.has_key(junction):
                temp2 = 0
            else:
                temp2 = tempDic2[junction]
            writer.write(str(junction) + '\t' + str(temp1) + '\n')
            #writer_junctionUG.write(str(junction) + '\t' + str(tempDic2[junction]) + '\t')
            #if junction == -13:
            #    writer_junctionUG.write(str(temp) + '\n')
            #else:
            #    writer_junctionUG.write('0' + '\n')
            writer_junctionUG.write(str(junction) + '\t' + str(temp2) + '\t' + str(tempDic3[junction]) + '\n')
        writer.close()
        writer_junctionUG.close()
    '''
    
    # write mid stem length distribution ( Drosha cleavage site ~ loop )
    if rev == False:
        writer = open(outputD + '/LoopPositionDistribution_FP.txt','w')
    else:
        writer = open(outputD + '/LoopPositionDistribution_TP.txt','w')
    writer.write('LoopPosition' + '\t' + 'Number' + '\n')
    loopPositionDic = dict()
    for id in allIndexDic.keys():
        loopStart, loopEnd = loopIndexDic[id]
        preStart, preEnd = allIndexDic[id]
        if rev == False:
            preStart, preEnd = allIndexDic[id]
            priSeq = priSeqDic[id]
            priStructure = structureDic[id]
        else:
            priSeq = priSeqDic[id][::-1]
            priStructure = structureDic[id][::-1].translate(trans)
            tempPreStart, tempPreEnd = allIndexDic[id]
            preStart = len(priSeq) - 1 - tempPreEnd
            preEnd = len(priSeq) - 1 - tempPreStart
        loopIndexIter = loopCompile.finditer(priStructure)
        loopIndex = map(lambda x: x.span(), loopIndexIter)[0]
        
        loopPosition = loopIndex[0] - preStart
        if not loopPositionDic.has_key(loopPosition): loopPositionDic[loopPosition] = 0
        loopPositionDic[loopPosition] += 1

    for loopPosition in sorted(loopPositionDic.keys()):
        writer.write(str(loopPosition) + '\t' + str(loopPositionDic[loopPosition]) + '\n')
    writer.close()
    
    # write bulge distribution
    minimum = []; maximum = []
    minimum = -1 * max(map(lambda x: x[0], newIndexDic.values()))
    maximum = max(map(lambda x: len(priSeqDic[x]) - 1 - newIndexDic[x][0], newIndexDic.keys()))
    if analysisType == 'Original':
        pass
    else:
        minimum_noB = -1 * max(map(lambda x: x[0], newIndexDic_noB.values()))
        maximum_noB = max(map(lambda x: len(newPriSeqDic[x]) - 1 - newIndexDic_noB[x][0], newIndexDic_noB.keys()))

    # write loop sequences 
    # from -5 nt of loop index to +5 nt of loop

    if rev == False:
        writer_High = open(outputD + '/LoopSequence_FP_high.fa','w')
        writer_Low = open(outputD + '/LoopSequence_FP_low.fa','w')
    else:
        writer_High = open(outputD + '/LoopSequence_TP_high.fa','w')
        writer_Low = open(outputD + '/LoopSequence_TP_low.fa','w')


    #if rev == 'False':
    #    writer = open(outputD + '/LoopSequence_FP.fa','w')
    #elif rev == 'True':
    #    writer = open(outputD + '/LoopSequence_TP.fa','w')

    for id in frequencyDic.keys():
        if frequencyDic[id] == 'High':
            writer_High.write('>' + id + '\n' + priSeqDic[id][loopIndexDic[id][0] - 10: loopIndexDic[id][1] + 10] + '\n')
        else:
            writer_Low.write('>' + id + '\n' + priSeqDic[id][loopIndexDic[id][0] - 10: loopIndexDic[id][1]+10] + '\n')
    writer_High.close()
    writer_Low.close()
    
    # original miRNAs
    write_mismatchInfo(rev, 'Ori', outputD, frequencyDic, minimum, maximum, newIndexDic, loopIndexDic, priSeqDic, structureDic, leftBulgeDic_ori, rightBulgeDic_ori, leftInternalDic_ori, rightInternalDic_ori, leftInternalDic_Assy_ori, rightInternalDic_Assy_ori, leftInternalDic_Symm_ori, rightInternalDic_Symm_ori, earlyPassed)
    print "Original", minimum, maximum
    #print len(leftBulgeDic.keys())
    #exit()
    # sepearate asymmetric internal loop into symmetric internal loop and bulge
    if analysisType == 'noB':
        write_mismatchInfo(rev, 'AssyToBulge', outputD, frequencyDic, minimum, maximum, newIndexDic, loopIndexDic, priSeqDic, structureDic, leftBulgeDic, rightBulgeDic, leftInternalDic, rightInternalDic, leftInternalDic_Assy, rightInternalDic_Assy, leftInternalDic_Symm, rightInternalDic_Symm, [])
        print "Original", minimum, maximum

        # remove bulges
        write_mismatchInfo(rev, 'RemoveBulge', outputD, frequencyDic, minimum_noB, maximum_noB, newIndexDic_noB, loopIndexDic_noB, newPriSeqDic, newStructureDic, leftBulgeDic_noB, rightBulgeDic_noB, leftInternalDic_noB, rightInternalDic_noB, leftInternalDic_Assy_noB, rightInternalDic_Assy_noB, leftInternalDic_Symm_noB, rightInternalDic_Symm_noB, passed)
        print "Removal of bulge", minimum_noB, maximum_noB

    for id in newIndexDic:
        print id, newIndexDic[id]
    '''
    if rev == True:
        prefix = '1.Mismatch_TP'
    else:
        prefix = '1.Mismatch_FP'
    #writer = open('./result/1.Mismatch.txt','w')
    #writer_Bulge = open('./result/1.Mismatch_bulge.txt','w')
    #writer_Internal = open('./result/1.Mismatch_internal_loop.txt','w')
    #writer_Internal_symm = open('./result/1.Mismatch_internal_loop_symmetric.txt','w')
    #writer_Internal_assy = open('./result/1.Mismatch_internal_loop_assymmetric.txt','w')
    writer = open(outputD + prefix + '.txt','w')
    writer_Bulge = open(outputD + prefix + '_bulge.txt','w')
    writer_Internal = open(outputD + prefix + '_internal_loop.txt','w')
    writer_Internal_symm = open(outputD + prefix + '_internal_loop_symmetric.txt','w')
    writer_Internal_assy = open(outputD + prefix + '_internal_loop_assymetric.txt','w')
    #header = 'Index' + '\t' + '\t'.join(map(str,range(0,maximum+1))) + '\n'
    header = 'Index' + '\t' + '\t'.join(map(str,range(minimum, maximum+1))) + '\n'
    writer.write(header)
    writer_Bulge.write(header)
    writer_Internal.write(header)
    writer_Internal_symm.write(header)
    writer_Internal_assy.write(header)
    for id in allIndexDic.keys():
        if id in passed: continue
        if rev == False:
            preStart, preEnd = allIndexDic[id]
            priSeq = priSeqDic[id]
            priStructure = structureDic[id]
        else:
            priSeq = priSeqDic[id][::-1]
            priStructure = structureDic[id][::-1].translate(trans)
            tempPreStart, tempPreEnd = allIndexDic[id]
            preStart = len(priSeq) - 1 - tempPreEnd
            preEnd = len(priSeq) - 1 - tempPreStart
        
        #loopIndex = loopIndexDic[id]
        loopIndexIter = loopCompile.finditer(priStructure)
        loopIndex = map(lambda x: x.span(), loopIndexIter)[0]
        writer.write(id)
        writer_Bulge.write(id)
        writer_Internal.write(id)
        writer_Internal_symm.write(id)
        writer_Internal_assy.write(id)
        for n in range(minimum, maximum+1):
        #for i in range(0,maximum+1):
            i = n + preStart
            if i < 0 or i > len(priSeqDic[id]) - 1:
                writer.write('\t' + '1')
            else:
                if priStructure[i] == '.':
                    writer.write('\t' + '1')
                else:
                    writer.write('\t' + '0')
            if i < 0 or i > len(priSeqDic[id]) - 1:
                #writer.write('\t' + '1')
                writer_Bulge.write('\t' + '0')
                writer_Internal.write('\t' + '0')
                writer_Internal_symm.write('\t' + '0')
                writer_Internal_assy.write('\t' + '0')
            elif loopIndex[0] <= i < loopIndex[1]:
                #writer.write('\t' + '1')
                writer_Bulge.write('\t' + '0')
                writer_Internal.write('\t' + '0')
                writer_Internal_symm.write('\t' + '0')
                writer_Internal_assy.write('\t' + '0')
            else:
                #if i in leftBulgeDic[id] or i in rightBulgeDic[id] or len(filter(lambda x: i in x, leftInternalDic[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic[id])) > 0:
                #    #writer.write('\t' + '1')
                #else:
                #    #writer.write('\t' + '0')
                if i in leftBulgeDic[id] or i in rightBulgeDic[id]:
                    writer_Bulge.write('\t' + '1')
                    #writer.write('\t' + '1')
                else:
                    writer_Bulge.write('\t' + '0')
                    #writer.write('\t' + '0')
                if len(filter(lambda x: i in x, leftInternalDic[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic[id])) > 0:
                    writer_Internal.write('\t' + '1')
                else:
                    writer_Internal.write('\t' + '0')
                if len(filter(lambda x: i in x, leftInternalDic_Assy[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic_Assy[id])) > 0:
                    #writer_Internal.write('\t' + '1')
                    writer_Internal_assy.write('\t' + '1')
                else:
                    #writer_Internal.write('\t' + '0')
                    writer_Internal_assy.write('\t' + '0')
                if len(filter(lambda x: i in x, leftInternalDic_Symm[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic_Symm[id])) > 0:
                    #writer_Internal.write('\t' + '1')
                    writer_Internal_symm.write('\t' + '1')
                else:
                    #writer_Internal.write('\t' + '0')
                    writer_Internal_symm.write('\t' + '0')
        writer.write('\n')
        writer_Bulge.write('\n')
        writer_Internal.write('\n')
        writer_Internal_assy.write('\n')
        writer_Internal_symm.write('\n')
    writer_Bulge.close()
    writer_Internal.close()
    writer_Internal_symm.close()
    writer_Internal_assy.close()
    '''
    print [loopIndexDic[x][1] - loopIndexDic[x][0] for x in frequencyDic.keys() if frequencyDic[x] == 'High']
    print [loopIndexDic[x][1] - loopIndexDic[x][0] for x in frequencyDic.keys() if frequencyDic[x] == 'Low']

    print stats.mannwhitneyu([loopIndexDic[x][1] - loopIndexDic[x][0] for x in frequencyDic.keys() if frequencyDic[x] == 'High'],[loopIndexDic[x][1] - loopIndexDic[x][0] for x in frequencyDic.keys() if frequencyDic[x] == 'Low'])
    return finalSequenceDic_37, finalStructureDic_37, finalSequenceDic_Loop, finalSequenceDic_All, finalStructureDic_All, loopIndexDic_noB, newIndexDic_noB, myStructureDic, myStructureDic_37

def Find_mismatch(loopIndex, priStructure):
    match = []
    leftBulge = []; rightBulge = []
    leftInternal = []; rightInternal = []
    leftInternal_Symm = []; rightInternal_Symm = []
    leftInternal_Assy = []; rightInternal_Assy = []

    tempLeftIndex = loopIndex[0]; tempRightIndex = loopIndex[1] - 1
    tempLeftBulge = []; tempRightBulge = []
    hit = 0
    #for i in xrange(loopIndex[0]+5):
    for i in xrange(len(priStructure)):
        print i
        print '\tMatch :', match
        print '\tBulge :', leftBulge, rightBulge
        print '\tInternal Loop :', leftInternal, rightInternal
        if hit == priStructure.count('('): break
        if priStructure[tempLeftIndex] == '(' and priStructure[tempRightIndex] == ')':
            match.append([tempLeftIndex, tempRightIndex])
            if len(tempLeftBulge) == 0 and len(tempRightBulge) == 0:
                pass
            elif len(tempLeftBulge) > 0 and len(tempRightBulge) == 0:
                leftBulge += tempLeftBulge
            elif len(tempLeftBulge) == 0 and len(tempRightBulge) > 0:
                rightBulge += tempRightBulge
            else:
                if len(tempLeftBulge) == len(tempRightBulge):
                    leftInternal_Symm.append(tempLeftBulge)
                    leftInternal.append(tempLeftBulge)
                    rightInternal_Symm.append(tempRightBulge)
                    rightInternal.append(tempRightBulge)
                else:
                    leftInternal_Assy.append(tempLeftBulge)
                    leftInternal.append(tempLeftBulge)
                    rightInternal_Assy.append(tempRightBulge)
                    rightInternal.append(tempRightBulge)
            tempRightBulge = []
            tempLeftBulge = []
            tempRightIndex += 1
            tempLeftIndex -= 1
            hit += 1
        elif priStructure[tempLeftIndex] == '(' and priStructure[tempRightIndex] == '.':
            if len(tempLeftBulge) > 0 and len(tempRightBulge) == 0:
                leftBulge += tempLeftBulge
                tempRightBulge.append(tempRightIndex)
                tempLeftBulge = []
            else:
                tempRightBulge.append(tempRightIndex)
            tempRightIndex += 1
        elif priStructure[tempLeftIndex] == '.' and priStructure[tempRightIndex] == ')':
            if len(tempLeftBulge) == 0 and len(tempRightBulge) > 0:
                rightBulge += tempRightBulge
                tempLeftBulge.append(tempLeftIndex)
                tempRightBulge = []
            else:
                tempLeftBulge.append(tempLeftIndex)
            tempLeftIndex -= 1
        else:
            tempLeftBulge.append(tempLeftIndex)
            tempRightBulge.append(tempRightIndex)
            tempRightIndex += 1
            tempLeftIndex -= 1
    return match, leftBulge, rightBulge, leftInternal, rightInternal, leftInternal_Symm, rightInternal_Symm, leftInternal_Assy, rightInternal_Assy

def write_mismatchInfo(rev, structureType, outputD, frequencyDic, minimum, maximum, indexDic, loopIndexDic, priSeqDic, structureDic, leftBulgeDic, rightBulgeDic, leftInternalDic, rightInternalDic, leftInternalDic_Assy, rightInternalDic_Assy, leftInternalDic_Symm, rightInternalDic_Symm, passed):
    if rev == True:
        prefix = 'MismatchInformation_TP_'+structureType# + '_'
    else:
        prefix = 'MismatchInformation_FP_'+structureType# + '_'
    
    writer_All = open(outputD + prefix + '.txt','w')
    writer_All_High = open(outputD + prefix + '_high.txt','w')
    writer_All_Low = open(outputD + prefix + '_low.txt','w')
    #writer_Junction = open(outputD + prefix + '.txt','w')
    #writer_Junction_High = open(outputD + prefix + '_high.txt','w')
    #writer_Junction_Low = open(outputD + prefix + '_Low.txt','w')
    writer_Bulge = open(outputD + prefix + '_bulge.txt','w')
    writer_Bulge_High = open(outputD + prefix + '_bulge_high.txt','w')
    writer_Bulge_Low = open(outputD + prefix + '_bulge_low.txt','w')
    writer_without_Bulge = open(outputD + prefix + '_withoutBulge.txt','w')
    writer_without_Bulge_High = open(outputD + prefix + '_withoutBulge_high.txt','w')
    writer_without_Bulge_Low = open(outputD + prefix + '_withoutBulge_low.txt','w')
    writer_Internal = open(outputD + prefix + '_internal_loop.txt','w')
    writer_Internal_High = open(outputD + prefix + '_internal_loop_high.txt','w')
    writer_Internal_Low = open(outputD + prefix + '_internal_loop_low.txt','w')
    writer_Internal_symm = open(outputD + prefix + '_internal_loop_symmetric.txt','w')
    writer_Internal_symm_High = open(outputD + prefix + '_internal_loop_symmetric_high.txt','w')
    writer_Internal_symm_Low = open(outputD + prefix + '_internal_loop_symmetric_low.txt','w')
    writer_Internal_assy = open(outputD + prefix + '_internal_loop_asymmetric.txt','w')
    writer_Internal_assy_High = open(outputD + prefix + '_internal_loop_asymmetric_high.txt','w')
    writer_Internal_assy_Low = open(outputD + prefix + '_internal_loop_asymmetric_low.txt','w')
    writer_without_Bulge_Internal = open(outputD + prefix + '_withoutBulgeInternalLoop.txt','w')
    writer_without_Bulge_Internal_High = open(outputD + prefix + '_withoutBulgeInternalLoop_high.txt','w')
    writer_without_Bulge_Internal_Low = open(outputD + prefix + '_withoutBulgeInternalLoop_low.txt','w')
    header = 'Index' + '\t' + '\t'.join(map(str,range(minimum,maximum+1))) + '\n'

    writer_All.write(header)
    writer_All_High.write(header)
    writer_All_Low.write(header)
    writer_Bulge.write(header)
    writer_Bulge_High.write(header)
    writer_Bulge_Low.write(header)
    writer_without_Bulge.write(header)
    writer_without_Bulge_High.write(header)
    writer_without_Bulge_Low.write(header)
    writer_Internal.write(header)
    writer_Internal_High.write(header)
    writer_Internal_Low.write(header)
    writer_Internal_symm.write(header)
    writer_Internal_symm_High.write(header)
    writer_Internal_symm_Low.write(header)
    writer_Internal_assy.write(header)
    writer_Internal_assy_High.write(header)
    writer_Internal_assy_Low.write(header)
    writer_without_Bulge_Internal.write(header)
    writer_without_Bulge_Internal_High.write(header)
    writer_without_Bulge_Internal_Low.write(header)
    for id in indexDic.keys():
        if id in passed: continue
        if not frequencyDic.has_key(id): continue
        if rev == False:
            preStart, preEnd = indexDic[id]
            priSeq = priSeqDic[id]
            priStructure = structureDic[id]
        else:
            priSeq = priSeqDic[id][::-1]
            priStructure = structureDic[id][::-1].translate(trans)
            #tempPreStart, tempPreEnd = indexDic[id]
            preStart, preEnd = indexDic[id]
            #preStart = len(priSeq) - 1 - tempPreEnd
            #preEnd = len(priSeq) - 1 - tempPreStart
        #if id == 'hsa-mir-3591':
        #    print indexDic[id], preStart, preEnd
        #    exit()
        loopIndexIter = loopCompile.finditer(priStructure)
        loopIndex = map(lambda x: x.span(), loopIndexIter)[0]
        writer_All.write(id)
        writer_Bulge.write(id)
        writer_without_Bulge.write(id)
        writer_Internal.write(id)
        writer_Internal_assy.write(id)
        writer_Internal_symm.write(id)
        writer_without_Bulge_Internal.write(id)
        if frequencyDic[id] == 'High':
            writer_All_High.write(id)
            writer_Bulge_High.write(id)
            writer_without_Bulge_High.write(id)
            writer_Internal_High.write(id)
            writer_Internal_assy_High.write(id)
            writer_Internal_symm_High.write(id)
            writer_without_Bulge_Internal_High.write(id)
        elif frequencyDic[id] == 'Low':
            writer_All_Low.write(id)
            writer_Bulge_Low.write(id)
            writer_without_Bulge_Low.write(id)
            writer_Internal_Low.write(id)
            writer_Internal_assy_Low.write(id)
            writer_Internal_symm_Low.write(id)
            writer_without_Bulge_Internal_Low.write(id)
        else:
            print id, frequencyDic[id]
            exit()
        for n in range(minimum, maximum+1):
            i = n+preStart
            hit = '\t' + '1'
            no = '\t' + '0'
            if i < 0 or i > len(priSeqDic[id]) - 1:
                writer_All.write(hit)
                writer_without_Bulge.write(hit)
                writer_without_Bulge_Internal.write(hit)
                if frequencyDic[id] == 'High':
                    writer_All_High.write(hit)
                    writer_without_Bulge_High.write(hit)
                    writer_without_Bulge_Internal_High.write(hit)
                elif frequencyDic[id] == 'Low':
                    writer_All_Low.write(hit)
                    writer_without_Bulge_Low.write(hit)
                    writer_without_Bulge_Internal_Low.write(hit)
            elif loopIndex[0] + 1 <= i < loopIndex[1] -1:
                writer_All.write(hit)
                writer_without_Bulge.write(hit)
                writer_without_Bulge_Internal.write(hit)
                if frequencyDic[id] == 'High':
                    writer_All_High.write(hit)
                    writer_without_Bulge_High.write(hit)
                    writer_without_Bulge_Internal_High.write(hit)
                elif frequencyDic[id] == 'Low':
                    writer_All_Low.write(hit)
                    writer_without_Bulge_Low.write(hit)
                    writer_without_Bulge_Internal_Low.write(hit)
            else:
                if priStructure[i] == '.':
                    if i in leftBulgeDic[id] or i in rightBulgeDic[id]:
                        writer_without_Bulge.write(no)
                        writer_without_Bulge_Internal.write(no)
                        if frequencyDic[id] == 'High':
                            writer_without_Bulge_High.write(no)
                            writer_without_Bulge_Internal_High.write(no)
                        else:
                            writer_without_Bulge_Low.write(no)
                            writer_without_Bulge_Internal_Low.write(no)
                    elif len(filter(lambda x: i in x, leftInternalDic[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic[id])) > 0:
                        writer_without_Bulge.write(hit)
                        writer_without_Bulge_Internal.write(no)
                        if frequencyDic[id] == 'High':
                            writer_without_Bulge_High.write(hit)
                            writer_without_Bulge_Internal_High.write(no)
                        else:
                            writer_without_Bulge_Low.write(hit)
                            writer_without_Bulge_Internal_Low.write(no)
                    else:
                        writer_without_Bulge.write(hit)
                        writer_without_Bulge_Internal.write(hit)
                        if frequencyDic[id] == 'High':
                            writer_without_Bulge_High.write(hit)
                            writer_without_Bulge_Internal_High.write(hit)
                        else:
                            writer_without_Bulge_Low.write(hit)
                            writer_without_Bulge_Internal_Low.write(hit)

                    writer_All.write(hit)
                    if frequencyDic[id] == 'High':
                        writer_All_High.write(hit)
                    else:
                        writer_All_Low.write(hit)
                else:
                    writer_All.write(no)
                    writer_without_Bulge.write(no)
                    writer_without_Bulge_Internal.write(no)
                    if frequencyDic[id] == 'High':
                        writer_All_High.write(no)
                        writer_without_Bulge_High.write(no)
                        writer_without_Bulge_Internal_High.write(no)
                    else:
                        writer_All_Low.write(no)
                        writer_without_Bulge_Low.write(no)
                        writer_without_Bulge_Internal_Low.write(no)

            if i < 0 or i > len(priSeqDic[id]) - 1:
                writer_Bulge.write(no)
                writer_Internal.write(no)
                writer_Internal_symm.write(no)
                writer_Internal_assy.write(no)
                if frequencyDic[id] == 'High':
                    writer_Bulge_High.write(no)
                    writer_Internal_High.write(no)
                    writer_Internal_symm_High.write(no)
                    writer_Internal_assy_High.write(no)
                else:
                    writer_Bulge_Low.write(no)
                    writer_Internal_Low.write(no)
                    writer_Internal_symm_Low.write(no)
                    writer_Internal_assy_Low.write(no)
            elif loopIndex[0] + 1 <= i < loopIndex[1] - 1:
                writer_Bulge.write(no)
                writer_Internal.write(no)
                writer_Internal_assy.write(no)
                writer_Internal_symm.write(no)
                if frequencyDic[id] == 'High':
                    writer_Bulge_High.write(no)
                    writer_Internal_High.write(no)
                    writer_Internal_assy_High.write(no)
                    writer_Internal_symm_High.write(no)
                else:
                    writer_Bulge_Low.write(no)
                    writer_Internal_Low.write(no)
                    writer_Internal_assy_Low.write(no)
                    writer_Internal_symm_Low.write(no)
            else:
                if i in leftBulgeDic[id] or i in rightBulgeDic[id]:
                    writer_Bulge.write(hit)
                    if frequencyDic[id] == 'High':
                        writer_Bulge_High.write(hit)
                    else:
                        writer_Bulge_Low.write(hit)
                else:
                    writer_Bulge.write(no)
                    if frequencyDic[id] == 'High':
                        writer_Bulge_High.write(no)
                    else:
                        writer_Bulge_Low.write(no)
                if len(filter(lambda x: i in x, leftInternalDic[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic[id])) > 0:
                    writer_Internal.write(hit)
                    if frequencyDic[id] == 'High':
                        writer_Internal_High.write(hit)
                    else:
                        writer_Internal_Low.write(hit)
                else:
                    writer_Internal.write(no)
                    if frequencyDic[id] == 'High':
                        writer_Internal_High.write(no)
                    else:
                        writer_Internal_Low.write(no)
                if len(filter(lambda x: i in x, leftInternalDic_Assy[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic_Assy[id])) > 0:
                    writer_Internal_assy.write(hit)
                    if frequencyDic[id] == 'High':
                        writer_Internal_assy_High.write(hit)
                    else:
                        writer_Internal_assy_Low.write(hit)
                else:
                    writer_Internal_assy.write(no)
                    if frequencyDic[id] == 'High':
                        writer_Internal_assy_High.write(no)
                    else:
                        writer_Internal_assy_Low.write(no)
                if len(filter(lambda x: i in x, leftInternalDic_Symm[id])) > 0 or len(filter(lambda x: i in x, rightInternalDic_Symm[id])) > 0:
                    writer_Internal_symm.write(hit)
                    if frequencyDic[id] == 'High':
                        writer_Internal_symm_High.write(hit)
                    else:
                        writer_Internal_symm_Low.write(hit)
                else:
                    writer_Internal_symm.write(no)
                    if frequencyDic[id] == 'High':
                        writer_Internal_symm_High.write(no)
                    else:
                        writer_Internal_symm_Low.write(no)
        writer_All.write('\n')
        writer_Bulge.write('\n')
        writer_without_Bulge.write('\n')
        writer_Internal.write('\n')
        writer_Internal_assy.write('\n')
        writer_Internal_symm.write('\n')
        writer_without_Bulge_Internal.write('\n')
        if frequencyDic[id] == 'High':
            writer_All_High.write('\n')
            writer_Bulge_High.write('\n')
            writer_without_Bulge_High.write('\n')
            writer_Internal_High.write('\n')
            writer_Internal_assy_High.write('\n')
            writer_Internal_symm_High.write('\n')
            writer_without_Bulge_Internal_High.write('\n')
        else:
            writer_All_Low.write('\n')
            writer_Bulge_Low.write('\n')
            writer_without_Bulge_Low.write('\n')
            writer_Internal_Low.write('\n')
            writer_Internal_assy_Low.write('\n')
            writer_Internal_symm_Low.write('\n')
            writer_without_Bulge_Internal_Low.write('\n')


    writer_All.close()
    writer_All_High.close()
    writer_All_Low.close()
    writer_Bulge.close()
    writer_Bulge_High.close()
    writer_Bulge_Low.close()
    writer_without_Bulge.close()
    writer_without_Bulge_High.close()
    writer_without_Bulge_Low.close()
    writer_Internal.close()
    writer_Internal_High.close()
    writer_Internal_Low.close()
    writer_Internal_symm.close()
    writer_Internal_symm_High.close()
    writer_Internal_symm_Low.close()
    writer_Internal_assy.close()
    writer_Internal_assy_High.close()
    writer_Internal_assy_Low.close()
    writer_without_Bulge_Internal.close()
    writer_without_Bulge_Internal_High.close()
    writer_without_Bulge_Internal_Low.close()

    return

def write_Mismatch(newIndexDic, loopIndexDic, maximum, leftBulgeDic, rightBulgeDic, leftInternalDic, rightInternalDic, leftInternalDic_Symm, rightInternalDic_Symm, leftInternalDic_Assy, rightInternalDic_Assy, outputD, prefix, suffix):
    writer = open(outputD + '/' + prefix + '_' + suffix + '.txt','w')
    writer_Bulge = open(outputD + '/' + prefix + '_bulge_' + suffix + '.txt','w')
    writer_Internal = open(outputD + '/' + prefix + '_internal_loop_' + suffix + '.txt','w')
    writer_Internal_symm = open(outputD + '/' + prefix + '_internal_loop_symmetric_' + suffix + '.txt','w')
    writer_Internal_assy = open(outputD + '/' + prefix + '_internal_loop_asymmetric_' + suffix + '.txt','w')
    
    header = 'Index' + '\t' + '\t'.join(map(str,range(0,maximum+1)))
    writer.write(header + '\n')
    writer_Bulge.write(header + '\n')
    writer_Internal.write(header + '\n')
    writer_Internal_symm.write(header + '\n')
    writer_Internal_assy.write(header + '\n')

    for id in newIndexDic.keys():
        preStart, preEnd = newIndexDic[id]

        loopIndex = loopIndexDic[id]
        writer.write(id)
        writer_Bulge.write(id)
        writer_Internal.write(id)
        writer_Internal_symm.write(id)
        writer_Internal_assy.write(id)

        #for i in range(0,maximum+1):
        #    if i < 0 or i > len(priSeqDic[id]) - 1:

    writer.close()
    writer_Bulge.close()
    writer_Internal.close()
    writer_Internal_symm.close()
    writer_Internal_assy.close()
    
def write_result(priSeqDic, finalSequenceDic_37, finalStructureDic_37, finalSequenceDic_Loop, heterogeneityDic, frequencyDic, medianFrequency, outputD, myStructureDic, myStructureDic_37):

    trans2 = maketrans('AGCU','UCGA')
    # make table for finding total motifs
    maxLen = max(map(lambda x: len(finalSequenceDic_Loop[x][0]), finalSequenceDic_Loop.keys()))
    totalValues = []
    writer_dist = open(outputD + '/ValueDist.txt','w')
    writer_dist.write('Value' + '\n')
    for i in xrange(3):
        moduleN = i + 1
        for n in xrange(0, maxLen - i):
            combinations = map(lambda x: ''.join(x), sorted(itertools.product('AGCU',repeat=moduleN)))
            antisenses = map(lambda x: x.translate(trans2), combinations)
            tempDic = dict()
            for combination in combinations:
                for antisense in antisenses:
                    tempDic[combination + '_' + antisense] = []
            for miRNA in finalSequenceDic_Loop.keys():
                if not heterogeneityDic.has_key(miRNA): continue
                seqFP, seqTP = finalSequenceDic_Loop[miRNA]
                seqTP_rv = seqTP[::-1]

                #if len(seqFP) >= n+moduleN: continue
                tempSeqFP = seqFP[n:n+moduleN]
                tempSeqTP = seqTP_rv[n:n+moduleN]
                #print miRNA, n, n+moduleN,tempSeqFP, tempSeqTP
                #print seqFP, seqTP
                if len(tempSeqFP) != moduleN: continue
                heterogeneity = heterogeneityDic[miRNA][0] # Catholic frequency
                #totalValues.append(math.log(heterogeneity,2) - math.log(medianFrequency,2))
                #writer_dist.write(str(math.log(heterogeneity,2) - math.log(medianFrequency,2)) + '\n')
                #tempDic[tempSeqFP + '_' + tempSeqTP].append(math.log(heterogeneity,2) - math.log(medianFrequency,2))
                #print miRNA, heterogeneity
                totalValues.append(heterogeneity)
                writer_dist.write(str(heterogeneity) + '\n')
                tempDic[tempSeqFP + '_' + tempSeqTP].append(heterogeneity)
                #totalValues.append(heterogeneity - medianFrequency)
                #writer_dist.write(str(heterogeneity - medianFrequency) + '\n')
                #tempDic[tempSeqFP + '_' + tempSeqTP].append(heterogeneity - medianFrequency)
                #if heterogeneity - medianFrequency > 5:
                #    print miRNA
            if not os.path.exists(outputD + '/OffsetTotal/'): os.makedirs(outputD + '/OffsetTotal/')
            writer = open(outputD + '/OffsetTotal/All_N' + str(moduleN) + '-' + str(n+1) + '_FPheterogeneity.txt','w')
            writer_mean = open(outputD + '/OffsetTotal/miRNA_N' + str(moduleN) + '-' + str(n+1) + '_FPheterogeneity.txt','w')
            writer.write('3p_5p' + '\t' + '\t'.join(map(lambda x: x.replace('U','T'), combinations)) + '\n')
            writer_mean.write('3p_5p' + '\t' + '\t'.join(map(lambda x: x.replace('U','T'), combinations)) + '\n')
            for antisense in antisenses:
                temp = []; tempMean = []
                for combination in combinations:
                    if len(tempDic[combination + '_' + antisense]) > 0:
                        temp.append(','.join(map(str,tempDic[combination + '_' + antisense])))
                        tempMean.append(np.mean(tempDic[combination + '_' + antisense]))
                    else:
                        temp.append('NA')
                        tempMean.append('NA')
                writer.write(antisense.replace('U','T') + '\t' + '\t'.join(temp) + '\n')
                writer_mean.write(antisense.replace('U','T') + '\t' + '\t'.join(map(str,tempMean)) + '\n')
            writer.close()
            writer_mean.close()

    #print stats.mannwhitneyu([[loopIndexDic[x][1] - loopIndexDic[x][0] for x in frequencyDic.keys() if frequencyDic[x] == 'High'],[loopIndexDic[x][1] - loopIndexDic[x][0] for x in frequencyDic.keys() if frequencyDic[x] == 'Low']])
    print "Maximum :", max(totalValues)
    print "Minimum :", min(totalValues)
    writer_dist.close()

    bpDistDic = dict()

    for fp in ['A','G','C','U']:
        for tp in ['A','G','C','U']:
            bpDistDic[fp + '-' + tp] = {'High':[],'Low':[]}
            for i in xrange(60):
                bpDistDic[fp+'-'+tp]['High'].append(0)
                bpDistDic[fp+'-'+tp]['Low'].append(0)
    for id in finalSequenceDic_37.keys():
        if not frequencyDic.has_key(id): continue
        fpSeq, antiTpSeq = finalSequenceDic_37[id]
        #print fpSeq
        #print antiTpSeq
        tpSeq = antiTpSeq[::-1]
        #fpSeq = fpSeq[4:]
        #print fpSeq, tpSeq
        for i in xrange(len(fpSeq)):
            if frequencyDic[id] == 'High':
                bpDistDic[fpSeq[i] + '-' + tpSeq[i]]['High'][i] += 1
            elif frequencyDic[id] == 'Low':
                bpDistDic[fpSeq[i] + '-' + tpSeq[i]]['Low'][i] += 1
    
    writer_high_All = open(outputD + '/Cat_All_HighFrequency.fa','w')
    writer_low_All = open(outputD + '/Cat_All_LowFrequency.fa','w')
    writer_high_37 = open(outputD + '/Cat_37_HighFrequency.fa','w')
    writer_low_37 = open(outputD + '/Cat_37_LowFrequency.fa','w')
    writer_high_37_fp = open(outputD + '/Cat_37_HighFrequency_FP.fa','w')
    writer_high_37_tp = open(outputD + '/Cat_37_HighFrequency_TP.fa','w')
    writer_high_37_all = open(outputD + '/Cat_37_HighFrequency_All.fa','w')
    writer_high_37_bp = open(outputD + '/Cat_37_HighFrequency_BP.txt','w')
    writer_low_37_fp  = open(outputD + '/Cat_37_LowFrequency_FP.fa','w')
    writer_low_37_tp  = open(outputD + '/Cat_37_LowFrequency_TP.fa','w')
    writer_low_37_all = open(outputD + '/Cat_37_LowFrequency_All.fa','w')
    writer_low_37_bp = open(outputD + '/Cat_37_LowFrequency_BP.txt','w')
    writer_top_25_fp = open(outputD + '/Cat_37_Top25_FP.fa','w')
    writer_top_25_tp = open(outputD + '/Cat_37_Top25_TP.fa','w')
    writer_bot_75_fp = open(outputD + '/Cat_37_Bot75_FP.fa','w')
    writer_bot_75_tp = open(outputD + '/Cat_37_Bot75_TP.fa','w')
    writer_bot_25_fp = open(outputD + '/Cat_37_Bot25_FP.fa','w')
    writer_bot_25_tp = open(outputD + '/Cat_37_Bot25_TP.fa','w')

    writer_back_all = open(outputD + '/Cat_37_Background_All.fa','w')
    writer_back_fp = open(outputD + '/Cat_37_Background_FP.fa','w')
    writer_back_tp = open(outputD + '/Cat_37_Background_TP.fa','w')

    writer_high_37_bp.write('Position' + '\t' + 'Pair' + '\t' + 'Number' + '\n')
    writer_low_37_bp.write('Position' + '\t' + 'Pair' + '\t' + 'Number' + '\n')
    
    for id in finalSequenceDic_37.keys():
        writer_back_all.write('>' + id + '\n' + finalSequenceDic_37[id][0] + finalSequenceDic_37[id][1][::-1] + '\n')
        writer_back_fp.write('>' + id + '\n' + finalSequenceDic_37[id][0] + '\n')
        writer_back_tp.write('>' + id + '\n' + finalSequenceDic_37[id][1][::-1] + '\n')
        if len(finalSequenceDic_37[id][0]) != 60:
            print "FP is wrong !"
            print id, finalSequenceDic_37[id][0], finalSequenceDic_37[id][1][::-1]
        if len(finalSequenceDic_37[id][1][::-1]) != 60:
            print "TP is wrong !"
            print id, finalSequenceDic_37[id][0], finalSequenceDic_37[id][1][::-1]

    for id in frequencyDic.keys():
        if frequencyDic[id] == 'High':
            writer_high_All.write('>' + id + '\n' + priSeqDic[id] + '\n' + myStructureDic[id] + '\n')
            if finalSequenceDic_37.has_key(id):
                writer_high_37.write('>' + id + '\n' + ''.join(finalSequenceDic_37[id]) + '\n')
                writer_high_37_fp.write('>' + id + '\n' + finalSequenceDic_37[id][0] + '\n')
                writer_high_37_tp.write('>' + id + '\n' + finalSequenceDic_37[id][1][::-1] + '\n')
                #writer_high_37_all.write('>' + id + ' ' + finalStructureDic_37[id][0] + finalStructureDic_37[id][1] + '\n' + finalSequenceDic_37[id][0] + finalSequenceDic_37[id][1][::-1] + '\n')
                writer_high_37_all.write('>' + id + ' ' + myStructureDic_37[id][0] + myStructureDic_37[id][1][::-1] + '\n' + finalSequenceDic_37[id][0] + finalSequenceDic_37[id][1][::-1] + '\n')
        elif frequencyDic[id] == 'Low':
            writer_low_All.write('>' + id + '\n' + priSeqDic[id] + '\n' + myStructureDic[id] + '\n')
            if finalSequenceDic_37.has_key(id):
                writer_low_37.write('>' + id + '\n' + ''.join(finalSequenceDic_37[id]) + '\n')
                writer_low_37_fp.write('>' + id + '\n' + finalSequenceDic_37[id][0] + '\n')
                writer_low_37_tp.write('>' + id + '\n' + finalSequenceDic_37[id][1][::-1] + '\n')
                #writer_low_37_all.write('>' + id + ' ' + finalStructureDic_37[id][0] + finalStructureDic_37[id][1] + '\n' + finalSequenceDic_37[id][0] + finalSequenceDic_37[id][1][::-1] + '\n')
                writer_low_37_all.write('>' + id + ' ' + myStructureDic_37[id][0] + myStructureDic_37[id][1][::-1] + '\n' + finalSequenceDic_37[id][0] + finalSequenceDic_37[id][1][::-1] + '\n')

    for bp in bpDistDic.keys():
        #for i in xrange(len(bpDistDic[bp])):
        #    position = i + 1
        #    writer_37_bp.write(str(position) + '\t' + bp + '\t' + str(bpDistDic[bp][i]) + '\t' + frequencyDic[id
        for i in xrange(len(bpDistDic[bp]['High'])):
            position = i + 1
            writer_high_37_bp.write(str(position) + '\t' + bp + '\t' + str(bpDistDic[bp]['High'][i]) + '\n')
        for i in xrange(len(bpDistDic[bp]['Low'])):
            position = i+1
            writer_low_37_bp.write(str(position) + '\t' + bp + '\t' + str(bpDistDic[bp]['Low'][i]) + '\n')

    totalItems = heterogeneityDic.items()
    totalItems.sort(cmp1)
    writeN = 0
    for item in totalItems:
        id = item[0]
        if not finalSequenceDic_37.has_key(id): continue
        if writeN < 26:
            writer_top_25_fp.write('>' + id + '\n' + finalSequenceDic_37[id][0] + '\n')
            writer_top_25_tp.write('>' + id + '\n' + finalSequenceDic_37[id][1][::-1] + '\n')
        else:
            writer_bot_75_fp.write('>' + id + '\n' + finalSequenceDic_37[id][0] + '\n')
            writer_bot_75_tp.write('>' + id + '\n' + finalSequenceDic_37[id][1][::-1] + '\n')
        writeN += 1
    writeN2 = 0
    for item in totalItems[::-1]:
        id = item[0]
        if not finalSequenceDic_37.has_key(id): continue
        if writeN2 < 26:
            writer_bot_25_fp.write('>' + id + '\n' + finalSequenceDic_37[id][0] + '\n')
            writer_bot_25_tp.write('>' + id + '\n' + finalSequenceDic_37[id][1][::-1] + '\n')
        writeN2 += 1

    writer_high_All.close()
    writer_low_All.close()
    writer_high_37.close()
    writer_low_37.close()
    writer_high_37_fp.close()
    writer_high_37_tp.close()
    writer_high_37_all.close()
    writer_high_37_bp.close()
    writer_low_37_fp.close()
    writer_low_37_tp.close()
    writer_low_37_all.close()
    writer_low_37_bp.close()
    writer_top_25_fp.close()
    writer_top_25_tp.close()
    writer_bot_75_fp.close()
    writer_bot_75_tp.close()
    writer_bot_25_fp.close()
    writer_bot_25_tp.close()

    writer_back_all.close()
    writer_back_fp.close()
    writer_back_tp.close()

def write_knownMotif(finalSequenceDic_All, loopIndexDic_noB, newIndexDic_noB, frequencyDic, ghgDic, outputD):
    writer_basalUG_high = open(outputD + '/KnownMotifDistribution_basalUG_High.txt','w')
    writer_basalUG_low = open(outputD + '/KnownMotifDistribution_basalUG_Low.txt','w')
    writer_GHG_high = open(outputD + '/KnownMotifDistribution_mismatchedGHG_High.txt','w')
    writer_GHG_low = open(outputD + '/KnownMotifDistribution_mismatchedGHG_Low.txt','w')
    writer_apicalUGU_high = open(outputD + '/KnownMotifDistribution_apicalUGU_High.txt','w')
    writer_apicalUGU_low = open(outputD + '/KnownMotifDistribution_apicalUGU_Low.txt','w')
    writer_CNNC_high = open(outputD + '/KnownMotifDistribution_CNNC_High.txt','w')
    writer_CNNC_low = open(outputD + '/KnownMotifDistribution_CNNC_Low.txt','w')

    basalUGdist = range(-30, 29) #range(-24, -3)
    GHGdist = range(-28, 31) #range(-15,2)
    apicalUGUdist = range(-30, 28) #range(12,33)
    CNNCdist = range(-30, 27) #range(9,31)
    writer_basalUG_high.write('miRNA' + '\t' + '\t'.join(map(str,basalUGdist)) + '\n')
    writer_basalUG_low.write('miRNA' + '\t' + '\t'.join(map(str,basalUGdist)) + '\n')
    writer_GHG_high.write('miRNA' + '\t' + '\t'.join(map(lambda x: str(x),GHGdist)) + '\n')
    writer_GHG_low.write('miRNA' + '\t' + '\t'.join(map(lambda x: str(x),GHGdist)) + '\n')
    writer_apicalUGU_high.write('miRNA' + '\t' + '\t'.join(map(str,apicalUGUdist)) + '\n')
    writer_apicalUGU_low.write('miRNA' + '\t' + '\t'.join(map(str,apicalUGUdist)) + '\n')
    writer_CNNC_high.write('miRNA' + '\t' + '\t'.join(map(lambda x: str(x), CNNCdist)) + '\n')
    writer_CNNC_low.write('miRNA' + '\t' + '\t'.join(map(lambda x: str(x), CNNCdist)) + '\n')

    for id in finalSequenceDic_All.keys():
        if not frequencyDic.has_key(id): continue
        tempUG = []
        tempUGU = []
        tempGHG = []
        tempCNNC = []
        tempSequence = finalSequenceDic_All[id]
        preStart_noB, preEnd_noB = newIndexDic_noB[id]
        loopIndex_noB = loopIndexDic_noB[id]
        newSequence_noB_FP = tempSequence[preStart_noB-13:preStart_noB+23+1]
        newSequence_noB_TP = tempSequence[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1):loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37]
        #print id, preStart_noB, loopIndex_noB, loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37, newSequence_noB_TP
        #print tempSequence
        for tempIndex in basalUGdist:
            #if id == 'hsa-mir-379':
            #    print id, preStart_noB, loopIndex_noB
            #    print tempSequence
            #    print tempIndex, tempSequence[preStart_noB-tempIndex-1:preStart_noB-tempIndex+1]
            if tempSequence[preStart_noB+tempIndex:preStart_noB+tempIndex+2] == 'UG':
                tempUG.append(1)
            else:
                tempUG.append(0)
        for tempIndex in GHGdist:
            tempSense = tempSequence[preStart_noB+tempIndex:preStart_noB+tempIndex+3]
            #tempAntisense = tempSequence[loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37-1+tempIndex:loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1)+37-1+tempIndex+3]
            tempAntisense = tempSequence[preEnd_noB-tempIndex-2:preEnd_noB-tempIndex+1][::-1]
            #tempSense = newSequence_noB_FP[13+tempIndex:13+tempIndex+3]
            #tempAntisense = newSequence_noB_TP[::-1][13+tempIndex:13+tempIndex+3]
            #if tempSense[0] == 'C' and tempSense[2] == 'C' and tempAntisense[0] == 'G' and tempAntisense[1] in ['A','C','U'] and tempAntisense[2] == 'G': # strict !
            #if tempSense[0] + tempAntisense[0] in ['CG','UG'] and tempSense[1] + tempAntisense[1] in ['CU','UC','GA','GC','AU','AC','UA'] and tempSense[2] + tempAntisense[2] in ['CG','GC','AU','UA']: # like bartel
            if ghgDic[tempSense + tempAntisense] >= 60:
            #if tempSense[0] + tempAntisense[0] in ['CG','UG'] and \
            #   tempSense[1] + tempAntisense[1] in ['UC','CU','GA','GC','AU','AC','UA'] and \
            #   tempSense[2] + tempAntisense[2] in ['CG','GC','AU','UA']:
                tempGHG.append(1)
            else:
                tempGHG.append(0)
            #if id == 'hsa-mir-101-1':
                #print tempSequence[preStart_noB-tempIndex-3:preStart_noB-tempIndex+1]
            #    print tempSequence
            #    print preStart_noB, preEnd_noB, tempIndex, loopIndex_noB
            #    print tempSense, tempAntisense
            #    print tempSense[0] + tempAntisense[0], tempSense[1] + tempAntisense[1], tempSense[2] + tempAntisense[2]
                #exit()
                                                
        for tempIndex in apicalUGUdist:
            #if id == 'hsa-mir-641':
            #    print tempIndex, tempSequence[preStart_noB+tempIndex-1:preStart_noB+tempIndex+2]
            #    print tempSequence
            #    print preStart_noB
            if tempSequence[preStart_noB+tempIndex-1:preStart_noB+tempIndex+2] in ['UGU']:
                tempUGU.append(1)
            else:
                tempUGU.append(0)
        for tempIndex in CNNCdist:
            #tempSeq = tempSequence[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1)+37-1+tempIndex:loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1)+37-1+tempIndex+4] 
            tempSeq = tempSequence[preEnd_noB+tempIndex:preEnd_noB+tempIndex+4]
            if id == 'hsa-mir-27a':
                print tempSequence
                print preStart_noB, preEnd_noB, tempIndex, loopIndex_noB
                print tempSeq, tempIndex
                #exit()
            if len(tempSeq) < 4:
                tempCNNC.append(0)
                print id, "is too short to calc CNNC"
                print tempSeq, tempSequence
            else:
                if tempSeq[0] == tempSeq[3] == 'C':
                    tempCNNC.append(1)
                else:
                    tempCNNC.append(0)
        if frequencyDic[id] == 'High':
            writer_basalUG_high.write(id + '\t' + '\t'.join(map(str,tempUG)) + '\n')
            writer_GHG_high.write(id + '\t' + '\t'.join(map(str,tempGHG)) + '\n')
            writer_apicalUGU_high.write(id + '\t' + '\t'.join(map(str,tempUGU)) + '\n')
            writer_CNNC_high.write(id + '\t' + '\t'.join(map(str,tempCNNC)) + '\n')
        elif frequencyDic[id] == 'Low':
            writer_basalUG_low.write(id + '\t' + '\t'.join(map(str,tempUG)) + '\n')
            writer_GHG_low.write(id + '\t' + '\t'.join(map(str,tempGHG)) + '\n')
            writer_apicalUGU_low.write(id + '\t' + '\t'.join(map(str,tempUGU)) + '\n')
            writer_CNNC_low.write(id + '\t' + '\t'.join(map(str,tempCNNC)) + '\n')

    writer_basalUG_high.close()
    writer_basalUG_low.close()
    writer_GHG_high.close()
    writer_GHG_low.close()
    writer_apicalUGU_high.close()
    writer_apicalUGU_low.close()
    writer_CNNC_high.close()
    writer_CNNC_low.close()
    exit()

def write_CandMotif(finalSequenceDic_All, loopIndexDic_noB, newIndexDic_noB, frequencyDic, outputD):
    writer_high = open(outputD + '/Candidate_Motifs_inHigh.txt','w')
    writer_low = open(outputD + '/Candidate_Motifs_inLow.txt','w')
    header = 'miRNA' + '\t' + '-12::C-G' + '\t' + '1::U-N' + '\t' + '5::R-U' + '\t' + '12::GHY-YNG' + '\t' + '20::G-C'
    writer_high.write(header + '\n')
    writer_low.write(header + '\n')
    for id in finalSequenceDic_All.keys():
        if not frequencyDic.has_key(id): continue
        tempSequence = finalSequenceDic_All[id]
        preStart_noB, preEnd_noB = newIndexDic_noB[id]
        loopIndex_noB = loopIndexDic_noB[id]
        newSequence_noB_FP = tempSequence[preStart_noB-13:preStart_noB+23+1]
        newSequence_noB_TP = tempSequence[loopIndex_noB[1] + loopIndex_noB[0] - (preStart_noB+23+1):loopIndex_noB[1] + loopIndex_noB[0]-(preStart_noB+23+1)+37]
        
        tempList = []
        # -12::C-G
        if newSequence_noB_FP[1] == 'C' and newSequence_noB_TP[::-1][1] == 'G':
            tempList.append(1)
        else:
            tempList.append(0)
        # 1::T-N
        if newSequence_noB_FP[13] == 'U':
            tempList.append(1)
        else:
            tempList.append(0)
        # 5::R-U
        if newSequence_noB_FP[17] in ['A','G'] and newSequence_noB_TP[::-1][17] == 'U':
            tempList.append(1)
        else:
            tempList.append(0)
        # 12::GHY-YNG
        if newSequence_noB_FP[24] == 'G' and newSequence_noB_FP[25] != 'G' and newSequence_noB_FP[26] in ['C','U'] and newSequence_noB_TP[::-1][24] in ['C','U'] and newSequence_noB_TP[::-1][26] == 'G':
            tempList.append(1)
        else:
            tempList.append(0)
        # 20::G-C
        if newSequence_noB_FP[32] == 'G' and newSequence_noB_TP[::-1][32] == 'C':
            tempList.append(1)
        else:
            tempList.append(0)
        
        if frequencyDic[id] == 'High':
            writer_high.write(id + '\t' + '\t'.join(map(str,tempList)) + '\n')
        elif frequencyDic[id] == 'Low':
            writer_low.write(id + '\t' + '\t'.join(map(str,tempList)) + '\n')

    writer_high.close()
    writer_low.close()

if __name__=='__main__':
    import os
    import re
    import sys
    import math
    import itertools
    import numpy as np
    from scipy import stats
    from string import maketrans
    #foldF = '/home/seokju/Project/LiverCancer/src_rev1/miRDeep2/Isomir/Stratification/HairpIndex_MFE_new/SecondGeneratedSeqs_result'
    dataType = sys.argv[1] # Catholic | TCGA | Tsinghwa | All
    gradeType = sys.argv[2] # Normal | Cancer | Total
    analysisType = sys.argv[3] # Original | noB
    offsetDirection = sys.argv[4] # plus minus both
    #dataType = 'All'
    #gradeType = 'Total'
    #if gradeType == 'Normal':
    heterogeneityF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/'+gradeType+'/Frequency/AlleleFrequencyDist_OffsetPlusMiNum1_Frequency_IntersectSet_FDR.txt'
    ghgF = '/home/seokju/Project/Structure/dataset/mGHG_score.txt'
    #print heterogeneityF
    #exit()
    if dataType == 'Catholic':
        outputD = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5/'+gradeType+'/Frequency/result/'
        #heterogeneityF =
    elif dataType == 'TCGA':
        outputD = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/5.IsomiR/MotifAnalysis5/'+gradeType+'/Frequency/result/'
        #heterogeneityF = 
    elif dataType == 'Tsinghwa':
        outputD = '/home/seokju/Project/LiverCancer/Tsinghua/1.miRNA/5.IsomiR/MotifAnalysis5/'+gradeType+'/Frequency/result/'
        #heterogeneityF = 
    elif dataType == 'All':
        outputD = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/MotifAnalysis5_cutoff10_noTCGA/'+gradeType+'/Frequency/All_correlation_' + analysisType +'_' + offsetDirection+ '/'
    else:
        print "You should choice among 'Catholic', 'TCGA', and 'Tsinghwa'"
        exit()
    if not os.path.exists(outputD): os.makedirs(outputD)
    #inputType = 2  # 1 : Cat only in intersected set 2: Cat only in Cat only set
    #if inputType == 1:
    #    outputD = '/home/seokju/Project/LiverCancer/src_rev1/miRDeep2/Isomir/Stratification/renew/result_1/'
    #    heterogeneityF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/TotalIsomirStats_cutoff10_new4/Drosha/Normal/Frequency/AlleleFrequencyDist_OffsetPlusMiNum1_Frequency_IntersectSet.txt'
    #else:
    #    outputD = '/home/seokju/Project/LiverCancer/src_rev1/miRDeep2/Isomir/Stratification/renew/result_2/'
    #    heterogeneityF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/5.IsomiR/TotalIsomirStats_cutoff10_new5/Drosha/Normal/Frequency/AlleleFrequencyDist_OffsetPlusMiNum1_Frequency_IntersectSet.txt'

    foldF = '/home/seokju/Project/LiverCancer/src_rev1/miRDeep2/Isomir/Stratification/renew/HairpIndex_MFE3/SecondGeneratedSeqs_result'
    pairF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/hsa.hg19.gff3'
    matureF = '/home/seokju/Project/LiverCancer/reference/mirbase_v21/miRNA.fa'

    trans = maketrans('()',')(')
    loopCompile = re.compile(r"\(\.+\)")
    mismatchCompile = re.compile(r"\.+")
    
    if not os.path.exists(outputD): os.makedirs(outputD)
    heterogeneityDic, frequencyDic, medianFrequency, changedDic = get_Heterogeneity(heterogeneityF, dataType, offsetDirection)
    ghgDic = get_mGHGscore(ghgF)
    #print len(frequencyDic.keys())
    #exit()
    priSeqDic, structureDic = get_Fold(foldF, heterogeneityDic)
    #print len(priSeqDic.keys())
    #exit()
    pairDic, coordDic = get_pairDic(pairF)
    matureSeqDic = get_Sequence(matureF)    
    #print len(changedDic.keys())
    indexDic, allIndexDic, cleavageIndexDic = anno_CleavageSite(priSeqDic, matureSeqDic, structureDic, pairDic, coordDic, changedDic, 'Drosha')
    #print allIndexDic['hsa-mir-17']
    #print len(allIndexDic.keys())
    #exit()
    #print allIndexDic['hsa-mir-15a']
    #exit()
    finalSequenceDic_37,finalStructureDic_37,finalSequenceDic_Loop, finalSequenceDic_All, finalStructureDic_All, loopIndexDic_noB, newIndexDic_noB, myStructureDic, myStructureDic_37 = seperate_Sequence(priSeqDic, structureDic, allIndexDic, False, outputD, heterogeneityDic, frequencyDic, analysisType)
    #finalSequenceDic_37,finalStructureDic_37,finalSequenceDic_Loop, finalSequenceDic_All, finalStructureDic_All, loopIndexDic_noB, newIndexDic_noB, myStructureDic, myStructureDic_37 = seperate_Sequence(priSeqDic, structureDic, allIndexDic, True, outputD, heterogeneityDic, frequencyDic, analysisType)
    #print myStructureDic
    #exit()
    #print len(set(frequencyDic.keys()))
    #print len(set(finalSequenceDic_37.keys()) & set(frequencyDic.keys()))
    #exit()

    # no overhang of 3P cleavage site
    # it is not good because 2nt overhang is better to be anntated in terms of basal junction.
    #outputD = './result_2/'
    #seperate_Sequence(priSeqDic, structureDic, cleavageIndexDic, True, outputD)
    #seperate_Sequence(priSeqDic, structureDic, cleavageIndexDic, False, outputD)

    write_result(priSeqDic, finalSequenceDic_37, finalStructureDic_37, finalSequenceDic_Loop, heterogeneityDic, frequencyDic, medianFrequency, outputD, myStructureDic, myStructureDic_37)
    write_knownMotif(finalSequenceDic_All, loopIndexDic_noB, newIndexDic_noB, frequencyDic, ghgDic, outputD)
    #write_CandMotif(finalSequenceDic_All, loopIndexDic_noB, newIndexDic_noB, frequencyDic, outputD)

    #print outputD
