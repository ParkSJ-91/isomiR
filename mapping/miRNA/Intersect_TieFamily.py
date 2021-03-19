import numpy as np

catTable = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
tcgaTable = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered.txt'
tsinghwaTable = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/miRDeep2_Expression_All_Qnorm_RPMFiltered.txt'
datatype = 'miRNA'

catIntersectTable = catTable.split('.txt')[0] + '_Intersect.txt'
tcgaIntersectTable = tcgaTable.split('.txt')[0] + '_Intersect.txt'
tsinghwaIntersectTable = tsinghwaTable.split('.txt')[0] + '_Intersect.txt'

catOutputD = '/'.join(catTable.split('/')[:-1]) + '/'
tcgaOutputD = '/'.join(tcgaTable.split('/')[:-1]) + '/'
tsinghwaOutputD = '/'.join(tsinghwaTable.split('/')[:-1]) + '/'

def get_ExpTable(expTable):
    expDic = dict(); medianDic = dict(); rankDic = dict()
    f = open(expTable); lines = f.readlines(); f.close()
    header = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        line = line.strip().split('\t')
        miRName = line[0]
        exp = map(float, line[1:])
        expDic[miRName] = exp
        medianDic[miRName] = np.median(exp)
    
    totalItems = medianDic.items()
    totalItems.sort(cmp1)
    totalItems = totalItems[::-1]
    length = len(totalItems)
    for i in xrange(len(totalItems)):
        item = totalItems.pop(0)
        rankDic[item[0]] = (i + 1) / length # weighted rank 
    return header, expDic, rankDic

def cmp1(a1,a2): return cmp(a1[1],a2[1])

catHeader, catDic, nothing = get_ExpTable(catTable)
tcgaHeader, tcgaDic, nothing = get_ExpTable(tcgaTable)
tsinghwaHeader, tsinghwaDic, nothing = get_ExpTable(tsinghwaTable)

catMi = catDic.keys()
tcgaMi = tcgaDic.keys()
tsinghwaMi = tsinghwaDic.keys()

print len(set(catMi) & set(tcgaMi) & set(tsinghwaMi))
print len(set(catMi) & set(tcgaMi))
writer_cat = open(catIntersectTable,'w')
writer_tcga = open(tcgaIntersectTable,'w')
writer_tsinghwa = open(tsinghwaIntersectTable,'w')

writer_cat.write(datatype + '\t' + '\t'.join(catHeader) + '\n')
writer_tcga.write(datatype + '\t' + '\t'.join(tcgaHeader) + '\n')
writer_tsinghwa.write(datatype + '\t' + '\t'.join(tsinghwaHeader) + '\n')

for name in set(catMi) & set(tcgaMi) & set(tsinghwaMi):
    writer_cat.write(name + '\t' + '\t'.join(map(str,catDic[name])) + '\n')
    writer_tcga.write(name + '\t' + '\t'.join(map(str,tcgaDic[name])) + '\n')
    writer_tsinghwa.write(name + '\t' + '\t'.join(map(str,tsinghwaDic[name])) + '\n')

writer_cat.close()
writer_tcga.close()
writer_tsinghwa.close()

if datatype == 'mRNA':
    exit()

catHeader,catExpDic, catRankDic = get_ExpTable(catIntersectTable)
tcgaHeader,tcgaExpDic, tcgaRankDic = get_ExpTable(tcgaIntersectTable)
tsinghwaHeader, tsinghwaExpDic, tsinghwaRankDic = get_ExpTable(tsinghwaIntersectTable)
type = catTable.split('/')[-1].split('_')[2]
catSeqF = '/'.join(catTable.split('/')[:-1]) + '/miRDeep2_Sequence.txt'

f = open(catSeqF); lines = f.readlines(); f.close()
seqDic = dict()
for line in lines:
    line = line.strip().split('\t')
    seqDic[line[0]] = line[1]


def get_Fam(expDic, seqDic):
    famDic = dict()
    for miRNA in expDic.keys():
        seed = seqDic[miRNA][1:8]
        if not famDic.has_key(seed): famDic[seed] = []
        famDic[seed].append(miRNA)
    return famDic
famDic = get_Fam(catExpDic, seqDic)

writer_catExp = open(catOutputD + 'miRDeep2_Family_7mer_Expression_'+type+'_Qnorm_RPMFiltered_Intersect.txt','w')
writer_catInfo = open(catOutputD + 'miRDeep2_Family_7mer_Infomation.txt','w')
writer_tcgaExp = open(tcgaOutputD + 'miRDeep2_Family_7mer_Expression_'+type+'_Qnorm_RPMFiltered_Intersect.txt','w')
writer_tcgaInfo = open(tcgaOutputD + 'miRDeep2_Family_7mer_Infomation.txt','w')
writer_tsinghwaExp = open(tsinghwaOutputD + 'miRDeep2_Family_7mer_Expression_'+type+'_Qnorm_RPMFiltered_Intersect.txt','w')
writer_tsinghwaInfo = open(tsinghwaOutputD + 'miRDeep2_Family_7mer_Infomation.txt','w')

writer_catExp.write('Representative_miRNA' + '\t' + '\t'.join(catHeader) + '\n')
writer_catInfo.write('Representative_miRNA' + '\t' + 'Seed' + '\t' + 'miRNAs' + '\n')
writer_tcgaExp.write('Representative_miRNA' + '\t' + '\t'.join(tcgaHeader) + '\n')
writer_tcgaInfo.write('Representative_miRNA' + '\t' + 'Seed' + '\t' + 'miRNAs' + '\n')
writer_tsinghwaExp.write('Representative_miRNA' + '\t' + '\t'.join(tsinghwaHeader) + '\n')
writer_tsinghwaInfo.write('Representative_miRNA' + '\t' + 'Seed' + '\t' + 'miRNAs' + '\n')

for seed in famDic.keys():
    catTotalExp = []
    tcgaTotalExp = []
    rank = []
    miRNAs = famDic[seed]
    for miRNA in miRNAs:
        rank.append([miRNA,catRankDic[miRNA]])
        rank.append([miRNA,tcgaRankDic[miRNA]])
        rank.append([miRNA,tsinghwaRankDic[miRNA]])
        if len(catTotalExp) == 0:
            catTotalExp = catExpDic[miRNA]
            tcgaTotalExp = tcgaExpDic[miRNA]
            tsinghwaTotalExp = tsinghwaExpDic[miRNA]
        else:
            for i in xrange(len(catExpDic[miRNA])):
                catTotalExp[i] += catExpDic[miRNA][i]

            for j in xrange(len(tcgaExpDic[miRNA])):
                tcgaTotalExp[j] += tcgaExpDic[miRNA][j]

            for k in xrange(len(tsinghwaExpDic[miRNA])):
                tsinghwaTotalExp[k] += tsinghwaExpDic[miRNA][k]

    rank.sort(cmp1)
    representativeMirna = rank[0][0]
    writer_catExp.write(representativeMirna + '\t' + '\t'.join(map(str,catTotalExp)) + '\n')
    writer_tcgaExp.write(representativeMirna + '\t' + '\t'.join(map(str,tcgaTotalExp)) + '\n')
    writer_tsinghwaExp.write(representativeMirna + '\t' + '\t'.join(map(str,tsinghwaTotalExp)) + '\n')
    writer_catInfo.write(representativeMirna + '\t' + seed + '\t' + '/'.join(miRNAs) + '\n')
    writer_tcgaInfo.write(representativeMirna + '\t' + seed + '\t' + '/'.join(miRNAs) + '\n')
    writer_tsinghwaInfo.write(representativeMirna + '\t' + seed + '\t' + '/'.join(miRNAs) + '\n')
writer_catExp.close()
writer_tcgaExp.close()
writer_tsinghwaExp.close()
writer_catInfo.close()
writer_tcgaInfo.close()
writer_tsinghwaInfo.close()
