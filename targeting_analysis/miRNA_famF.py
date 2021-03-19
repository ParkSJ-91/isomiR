######################################
###                                ###
###         made by seokju         ###
###           2016.09.10           ###
###                                ###
######################################

def get_sequenceDic(seqF, expF):
    expressedList = []; seqDic = dict()

    f = open(expF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        expressedMirna = line.split('\t')[0]
        expressedList.append(expressedMirna)

    f = open(seqF); lines = f.readlines(); f.close()
    for line in lines:
        line = line.strip('\n').split('\t')
        miRName = line[0]
        if not miRName in expressedList: continue
        sequence = line[1]
        seqDic[miRName] = sequence
    return seqDic

def get_miR_family(originalDic, heteroDic, mutationDic):
    famDic = dict()
    for miR_name in originalDic.keys():
        sequence = originalDic[miR_name]
        seed = sequence[1:8]
        if not famDic.has_key(seed): famDic[seed] = []
        famDic[seed].append(miR_name)
    for hetero_name in heteroDic.keys():
        sequence = heteroDic[hetero_name]
        seed = sequence[1:8]
        if not famDic.has_key(seed): famDic[seed] = []
        famDic[seed].append(hetero_name)
    for mutation_name in mutationDic.keys():
        sequence = mutationDic[mutation_name]
        seed = sequence[1:8]
        if not famDic.has_key(seed): famDic[seed] = []
        famDic[seed].append(mutation_name)
    return famDic

def get_miRFam(totalSeqDic):
    famDic= dict()
    for miRName in totalSeqDic.keys():
        seq = totalSeqDic[miRName]
        seed = seq[1:8]
        if not famDic.has_key(seed): famDic[seed] = []
        famDic[seed].append(miRName)
    return famDic

def get_conservedSeed(miR_FamInfoF, species):
    seedDic = dict()
    f = open(miR_FamInfoF); lines = f.readlines()[1:]; f.close()
    for line in lines:
        line = line.split('\t')
        seed = line[1]
        #species = line[2]
        if species == 'hsa':
            if line[2] == '9606': continue
        elif species == 'mmu':
            if line[2] == '10090': continue
        if not seedDic.has_key(seed): seedDic[seed] = []
        seedDic[seed].append(line[2])
    #for seed in seedDic.keys():
    print seedDic
    return seedDic

def write_result(outputD, totalSeqDic, seedDic, species):
    #total_miR_Dic = dict(originalDic.items() + heteroDic.items() + mutationDic.items())
    #famDic = get_miR_family(originalDic, heteroDic, mutationDic)
    #famDic = get_miR_family(totalSeqDic)
    famDic = get_miRFam(totalSeqDic)
    miR_fam_writer = open(outputD + 'miRNA_family_info.txt','w'); miR_seq_writer = open(outputD + 'miRNA_sequence_info.txt','w')
    if species == 'hsa':
        speciesID = '9606'
    elif species == 'mmu':
        speciesID = '10090'
    for seed in famDic.keys():
        miR_names = filter(lambda x: '-'.join(x.split('-')[1:]), famDic[seed])
        miR_names.sort()
        fam_name = '/'.join(miR_names)
        #species = sorted(list(set(seedDic[seed])))
        if seedDic.has_key(seed):
            species = sorted(list(set(seedDic[seed])))
            miR_fam_writer.write(fam_name + '\t' + seed + '\t' + speciesID + ';' + ';'.join(species) + '\n')
        else: miR_fam_writer.write(fam_name + '\t' + seed + '\t' + speciesID + '\n')
        for miR_name in famDic[seed]:
            miR_seq_writer.write(fam_name + '\t' + speciesID + '\t' + miR_name + '\t' + totalSeqDic[miR_name] + '\n')
    miR_fam_writer.close(); miR_seq_writer.close()

def main(outputD, seqF, expF, miR_FamInfoF, species):
    #originalDic = get_sequenceDic(originalF)
    #heteroDic = get_sequenceDic(heteroF)
    #mutationDic = get_sequenceDic(mutationF)
    totalSeqDic = get_sequenceDic(seqF,expF)
    seedDic = get_conservedSeed(miR_FamInfoF, species)
    write_result(outputD, totalSeqDic, seedDic, species)

if __name__=='__main__':
    import sys
    if len(sys.argv) != 5:
        #print 'Use example : miRNA_famF.py outputDirectory originalSequenceFile mutationSequenceFile heterogeneitySequenceFile'
        print "Usage : miRNA_famF.py output_directory type(catholic or TCGA) species(hsa or mmu)"
        print 'All file format : miRNA_name \t sequence'
        exit()
    expF = sys.argv[1]
    seqF = sys.argv[2]
    outputD = sys.argv[3]
    species = sys.argv[4]
    '''
    type = sys.argv[2]
    if type == 'catholic':
        seqF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling/miRDeep2_Sequence.txt'
        #expF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling/miRDeep2_Family_7mer_Expression_neoplasm_Qnorm_RPMFiltered_Intersect_RankSum_Intersect.txt'
        expF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling/miRDeep2_Expression_neoplasm_v2_Qnorm_RPMFiltered_Intersect.txt'
    elif type == 'TCGA':
        seqF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling/miRDeep2_Sequence.txt'
        expF = '/home/seokju/Project/LiverCancer/TCGA_rev1/1.miRNA/3.ExpressionProfiling/miRDeep2_Family_7mer_Expression_neoplasm_Qnorm_RPMFiltered_Intersect_RankSum_Intersect.txt'
    elif type == 'Intersection':
        seqF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling/miRDeep2_Sequence.txt'
        expF = '/home/seokju/Project/LiverCancer/Catholic/1.miRNA/3.ExpressionProfiling/miRDeep2_Expression_neoplasm_Qnorm_RPMFiltered_Intersect.txt'
    else:
        print "You should choice one between catholic or TCGA"
        exit()
    #originalF = sys.argv[2]
    #mutationF = sys.argv[3]
    #heteroF = sys.argv[4]
    '''
    miR_FamInfoF = '/home/seokju/Project/LiverCancer/reference/CompareTargetPredictionTools/ContextScore/miR_Family_Info.txt'
    main(outputD, seqF, expF, miR_FamInfoF, species)
    
