#import sys

#inputF = sys.argv[1]
#outputD = sys.argv[2]
#airF = '/home/seokju/Project/LiverCancer/src_rev1/ContextScore/AIRs_Huh7.txt'
#contextscoreF = '/home/seokju/Project/LiverCancer/Catholic/3.TargetAnalysis/7.TargetScan_v5/result/context_scores.txt'
#refFlatF = '/home/seokju/Data/gencode.v19.annotation.3PseqUpdated_rev1.refFlat'
#utrF = '/home/seokju/Project/LiverCancer/src_rev1/ContextScore/UTR_Sequences_IsoSelected.txt'

def get_utrF(utrF):
    utrDic = dict()
    for line in open(utrF):
        line = line.strip().split('\t')
        txnID = line[0]
        species = line[1]
        if not species == '9606': continue
        sequence = line[2].replace('-','')
        utrDic[txnID] = len(sequence)
    return utrDic
#print utrDic
def get_airF(airF):
    airDic = dict()
    for line in open(airF):
        #tempList = []
        line = line.strip().split('\t')
        txnID = line[0]
        start = int(line[1])
        end = int(line[2])
        air = float(line[3]) / 100
        if not airDic.has_key(txnID): airDic[txnID] = []
        airDic[txnID].append([start, air])
        airDic[txnID].append([end+1, air])
    return airDic

#airDic = get_airF(airF)

def get_matchDic(refFlatF):
    matchDic = dict()
    for line in open(refFlatF):
        line = line.strip().split('\t')
        geneSymbol = line[0]
        txnID = line[1]
        matchDic[txnID] = geneSymbol
    return matchDic

#matchDic = get_matchDic(refFlatF)

def get_contextscore(contextscoreF, airDic, matchDic, utrDic):
    scoreDic = dict(); newAirDic = dict(); typeDic = dict(); newUtrDic = dict()
    f = open(contextscoreF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        txnID = line[0]
        geneSymbol = matchDic[txnID]
        miRNA = line[2]
        siteType = line[3]
        start = int(line[4])
        end = int(line[5])
        score = float(line[30])
        #if score == 0: continue
        if not scoreDic.has_key(geneSymbol): scoreDic[geneSymbol] = dict()
        if not scoreDic[geneSymbol].has_key(miRNA): scoreDic[geneSymbol][miRNA] = 0
        scoreDic[geneSymbol][miRNA] += score
        if not newAirDic.has_key(geneSymbol): newAirDic[geneSymbol] = airDic[txnID]
        if not typeDic.has_key(geneSymbol): typeDic[geneSymbol] = dict()
        if not typeDic[geneSymbol].has_key(miRNA): typeDic[geneSymbol][miRNA] = []
        typeDic[geneSymbol][miRNA].append([start, end, siteType])
        if not newUtrDic.has_key(geneSymbol): newUtrDic[geneSymbol] = utrDic[txnID]
    return scoreDic, newAirDic, typeDic, newUtrDic

#scoreDic, newAirDic, typeDic = get_contextscore(contextscoreF, airDic, matchDic)

def get_Candidates(inputF):
    candDic = dict()
    f = open(inputF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.strip().split('\t')
        geneSymbol = line[0]
        miRNA = line[5]
        if not candDic.has_key(geneSymbol): candDic[geneSymbol] = []
        candDic[geneSymbol].append(miRNA)
    return candDic

if __name__=='__main__':
    import sys

    inputF = sys.argv[1]
    outputD = sys.argv[2]
    airF = '/home/seokju/Project/LiverCancer/src_rev1/ContextScore/AIRs_Huh7.txt'
    contextscoreF = '/home/seokju/Project/LiverCancer/Catholic/3.TargetAnalysis/7.TargetScan_noPVTT/result/context_scores.txt'
    refFlatF = '/home/seokju/Data/gencode.v19.annotation.sorted.refFlat'
    utrF = '/home/seokju/Project/LiverCancer/src_rev1/ContextScore/UTR_Sequences.txt'


    utrDic = get_utrF(utrF)
    print utrDic['ENST00000002165.6']
    #print utrDici['ENST00000302219.6']
    airDic = get_airF(airF)
    matchDic = get_matchDic(refFlatF)
    scoreDic, newAirDic, typeDic, newUtrDic = get_contextscore(contextscoreF, airDic, matchDic, utrDic)
    candDic = get_Candidates(inputF)

    for geneSymbol in candDic.keys():
        for miRNA in candDic[geneSymbol]:
            writer_site = open(outputD + geneSymbol + '_' + miRNA + '_site.txt', 'w')
            writer_site.write('miRNA' + '\t' + 'start' + '\t' + 'end' + '\t' + 'type' + '\n')
            for site in typeDic[geneSymbol][miRNA]:
                writer_site.write(miRNA + '\t' + '\t'.join(map(str, site)) + '\n')
            writer_site.close()
            writer_air = open(outputD + geneSymbol + '_' + miRNA + '_air.txt', 'w')
            writer_air.write('x' + '\t' + 'y' + '\n')
            for air in newAirDic[geneSymbol]:
                writer_air.write('\t'.join(map(str,air)) + '\n')
            writer_air.write(str(air[0]) + '\t' + '0' + '\n')
            writer_air.write(str(newUtrDic[geneSymbol]) + '\t' + '0' + '\n')
            writer_air.close()
        
