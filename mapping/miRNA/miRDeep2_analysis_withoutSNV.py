def get_sequenceDictionary(annotationFile):
    sequenceDictionary = dict()
    f = open(annotationFile); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        chunk = chunk.split('\n')
        miRName = chunk[0].split(' ')[0]
        sequence = chunk[1]
        sequenceDictionary[miRName] = sequence
    return sequenceDictionary

class mrd:
    def __init__(self, chunk, window):
        self._pairDic = dict(); self._sequenceDic = dict()
        self._lineDic = dict(); self._lineDic_monoU = dict()
        self._readDic = dict(); self._readDic_monoU = dict()
        self._doubleOffsetDic = dict()
        offsetNames = []; doubleOffsetNames = []
        lines = chunk.split('\n')
        names = map(lambda x: x.split(' ')[0], filter(lambda x: x.startswith('hsa-'),lines))
        self._preName = names[0]
        self._preSequence = filter(lambda x: len(x) > 0, filter(lambda x: x.startswith('pri_seq'),lines)[0].split(' '))[1]
        self._exp = filter(lambda x: len(x) > 0, filter(lambda x: x.startswith('exp'), lines)[0].split(' '))[1]
        self._matureNames = names[1:]

        self._pairDic[self._preName] = self._matureNames
        self._sequenceDic[self._preName] = self._preSequence

        reads = filter(lambda x: x.startswith('seq_'), lines)
        self._totalReads = sum(map(lambda x: int(x.split('_x')[1].split(' ')[0]), reads))
        reads_mmZero = filter(lambda x: int(x.split('\t')[1]) == 0, reads)
        reads_mmOne = filter(lambda x: int(x.split('\t')[1]) == 1, reads)
        reads_mmTwo = filter(lambda x: int(x.split('\t')[1]) == 2, reads)
        
        for matureName in self._matureNames:
            if matureName.endswith('-5p'):
                tempMarker = '5'
            elif matureName.endswith('-3p'):
                tempMarker = '3'
            else:
                tempMarker = 'M'
                
            for n in xrange(window*2+1):
                # get sequence !
                i = n - 5
                if self._exp.index(tempMarker) + i < 0: continue
                tempStartIndex = self._exp.index(tempMarker) + i
                tempEndIndex = self._exp.rindex(tempMarker) + i
                tempSequence = self._preSequence[tempStartIndex:tempEndIndex+1]
                # it is not seed ! just for compile heterogeneity miRNAs
                tempSeed = tempSequence[1:8] 
                if len(tempSequence) < 18: continue

                if i == 0:
                    new_miRName = matureName + '\t' + self._preName
                else:
                    new_miRName = matureName + '_' + tempSeed.upper() + '_' + str(i) + '\t' + self._preName
                    offsetNames.append(new_miRName)

                if not self._sequenceDic.has_key(new_miRName): self._sequenceDic[new_miRName] = tempSequence.upper()
                if not self._readDic.has_key(new_miRName): self._readDic[new_miRName] = 0
                if not self._readDic_monoU.has_key(new_miRName): self._readDic_monoU[new_miRName] = 0
                for tindex in xrange(window*2+1):
                    lastIndex = tindex - 5 
                    dOffsetName = matureName + '_' + str(i) + '_' + str(lastIndex) + '\t' + self._preName     
                    if not self._doubleOffsetDic.has_key(dOffsetName): self._doubleOffsetDic[dOffsetName] = 0
                if not self._lineDic.has_key(new_miRName): self._lineDic[new_miRName] = []
                if not self._lineDic_monoU.has_key(new_miRName): self._lineDic_monoU[new_miRName] = []
                # get read count !
                for read in reads_mmZero:
                    original = read
                    read = filter(lambda x: len(x)>0,read.split(' '))
                    read_name = read[0].split('_x')[0]
                    read_number = int(read[0].split('_x')[1])
                    read_sequence = read[1].split('\t')[0].strip('.')
                    read_sequence_start_in_read = read[1].index(read_sequence)
                    read_sequence_end_in_read = read[1].rindex(read_sequence) + len(read_sequence) - 1

                    if read_sequence_start_in_read == tempStartIndex:
                        self._readDic[new_miRName] += read_number
                        self._lineDic[new_miRName].append(original)

                        # get read count of double (3'UTR, 5'UTR) offset
                    for tindex in xrange(window*2+1):
                        lastIndex = tindex - 5
                        temptempEndIndex = self._exp.rindex(tempMarker) + lastIndex
                        dOffsetName = matureName + '_' + str(i) + '_' + str(lastIndex) + '\t' + self._preName
                        doubleOffsetNames.append(dOffsetName)

                        if temptempEndIndex == read_sequence_end_in_read and read_sequence_start_in_read == tempStartIndex:
                            self._doubleOffsetDic[dOffsetName] += read_number
                # get mono-uridylated read !
                for read in reads_mmOne:
                    original = read
                    read = filter(lambda x: len(x)>0,read.split(' '))
                    read_name = read[0].split('_x')[0]
                    read_number = int(read[0].split('_x')[1])
                    read_sequence = read[1].split('\t')[0].strip('.')
                    if read_sequence[-1] != 'U': continue
                    read_sequence_start_in_read = read[1].index(read_sequence)

                    if read_sequence_start_in_read == tempStartIndex:
                        self._readDic_monoU[new_miRName] += read_number
                        self._lineDic_monoU[new_miRName].append(original)

        tempNames = filter(lambda x: x.count('_') > 0, self._readDic.keys())
        alloffsetNames = filter(lambda x: len(x.split('_')[1]) != 2, tempNames)
        self._chunk = ''
        self._chunk += '>' + self._preName + '\n'
        for matureName in self._matureNames:
            if not matureName in self._readDic.keys(): continue
            offsetNames = filter(lambda x: x.split('_')[0] == matureName, alloffsetNames)
            self._chunk += matureName + '\t' + 'onset(perfect,monoU): ' + str(self._readDic[matureName]) + ',' + str(self._readDic_monoU[matureName]) + '\t' + 'offset(perfect,monoU): ' + str(sum(map(lambda x: self._readDic[x], offsetNames))) + ',' + str(sum(map(lambda x: self._readDic_monoU[x], offsetNames))) + '\n'
        self._chunk += self._preSequence + '\n'
        self._chunk += self._exp + '\n'
        for matureName in self._matureNames:
            if not matureName in self._readDic.keys(): continue
            for line in self._lineDic[matureName]:
                line = filter(lambda x: len(x) > 0, line.split('\t')[0].split(' '))
                read_number = line[0].split('_x')[1]
                query = line[1]
                self._chunk += query + '\t' + 'readCount: ' + read_number + '\t' + 'mismatch: 0' + '\t' + 'onset' + '\n'
            for line in self._lineDic_monoU[matureName]:
                line = filter(lambda x: len(x) > 0, line.split('\t')[0].split(' '))
                read_number = line[0].split('_x')[1]
                query = line[1]
                self._chunk += query + '\t' + 'readCount: ' + read_number + '\t' + 'mismatch: 1' + '\t' + 'onsetU' + '\n'
            for n in xrange(window*2+1):
                i = n-5
                if i == 0: continue
                new_miRName = filter(lambda x: x.split('_')[0] == matureName and x.count('_') > 0 and x.split('\t')[0].endswith('_'+str(i)), self._lineDic.keys())
                new_miRName = filter(lambda x: x.split('_')[1] != 2, new_miRName)
                if len(new_miRName) == 0: continue
                for line in self._lineDic[new_miRName[0]]:
                    line = filter(lambda x: len(x) > 0, line.split('\t')[0].split(' '))
                    read_number = line[0].split('_x')[1]
                    query = line[1]
                    self._chunk += query + '\t' + 'readCount: ' + read_number + '\t' + 'mismatch: 0' + '\t' + 'offset(' + str(i) + ')' + '\n'
                for line in self._lineDic_monoU[new_miRName[0]]:
                    line = filter(lambda x: len(x) > 0 , line.split('\t')[0].split(' '))
                    read_number = line[0].split('_x')[1]
                    query = line[1]
                    self._chunk += query + '\t' + 'readCount: ' + read_number + '\t' + 'mismatch: 1' + '\t' + 'offsetU(' + str(i) + ')' + '\n' 
    def preName(self): return self._preName
    def matureNames(self): return self._matureNames
    def preSequence(self): return self._preSequence
    def exp(self): return self._exp
    def totalReads(self): return self._totalReads
    def pairDic(self): return self._pairDic
    def sequenceDic(self): return self._sequenceDic
    def readDic(self): return self._readDic
    def readUDic(self): return self._readDic_monoU
    def doubleOffsetDic(self): return self._doubleOffsetDic
    def lineDic(self): return self._lineDic
    def lineUDic(self): return self._lineDic_monoU
    def chunk(self): return self._chunk

def get_miRDeep2(mrdF, window):
    totalSequenceDic = dict(); totalReadDic = dict(); totalMappedReads = 0#; totalLineDic = dict()
    new_chunks = []; totalReadMonoUDic = dict(); totalMappedReadsDic = dict()
    totalReadDoubleOffsetDic = dict()
    f = open(mrdF); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        # onset perfect expression profiling
        # heterogeneity perfect expression profiling
        # visualization above all
        parsed_chunk = mrd(chunk, window)
        totalMappedReadsOnPreMiRNA = parsed_chunk.totalReads()
        matureNames = parsed_chunk.matureNames()
        preName = parsed_chunk.preName()
        totalMappedReads += totalMappedReadsOnPreMiRNA
        for matureName in matureNames:
            name = matureName + '\t' + preName
            if not totalMappedReadsDic.has_key(name): totalMappedReadsDic[name] = 0
            totalMappedReadsDic[name] +=  totalMappedReadsOnPreMiRNA
        totalSequenceDic.update(parsed_chunk.sequenceDic())
        new_chunks.append(parsed_chunk.chunk())
        for miRName in parsed_chunk.readDic().keys():
            if not totalReadDic.has_key(miRName): totalReadDic[miRName] = 0
            totalReadDic[miRName] += parsed_chunk.readDic()[miRName]
            #if miRName == 'hsa-miR-3168':
            #    print totalReadDic['hsa-miR-3168']
        for miRName in parsed_chunk.readUDic().keys():
            if not totalReadMonoUDic.has_key(miRName): totalReadMonoUDic[miRName] = 0
            #print chunk
            #print parsed_chunk.readUDic()
            totalReadMonoUDic[miRName] += parsed_chunk.readUDic()[miRName]
        for miRName in parsed_chunk.doubleOffsetDic().keys():
            if not totalReadDoubleOffsetDic.has_key(miRName): totalReadDoubleOffsetDic[miRName] = 0
            totalReadDoubleOffsetDic[miRName] += parsed_chunk.doubleOffsetDic()[miRName]
    return totalSequenceDic, totalReadDic, totalMappedReads, new_chunks, totalReadMonoUDic, totalReadDoubleOffsetDic, totalMappedReadsDic

def write_miRDeep2(sample, outputD, totalSequenceDic, totalReadDic, totalMappedReads, new_chunks, totalReadMonoUDic, totalReadDoubleOffsetDic, totalMappedReadsDic):
    write_OnsetPerfect = open(outputD + sample + '_OnsetPerfect.txt','w')
    write_OnsetSNV = open(outputD + sample + '_OnsetSNV.txt','w')
    write_OffsetPerfect = open(outputD + sample + '_OffsetPerfect.txt','w')
    write_DoubleOffsetPerfect = open(outputD + sample + '_DoubleOffset.txt','w')
    write_OffsetSNV = open(outputD + sample + '_OffsetSNV.txt','w')
    write_Expression = open(outputD + sample + '_Expression.txt','w')
    writer_monoU = open(outputD + sample + '_monoU.txt','w')

    onsetperfectNames = sorted(filter(lambda x: x.count('_')==0, totalReadDic.keys()))
    offsetperfectNames = sorted(filter(lambda x: x.count('_')==2 and len(x.split('_')[1]) != 2, totalReadDic.keys()))
    for miRName in onsetperfectNames:
        write_OnsetPerfect.write(miRName + '\t' + totalSequenceDic[miRName] + '\t' + str(totalReadDic[miRName]) + '\t' + str(totalMappedReadsDic[miRName]) + '\t' + str(totalMappedReads) + '\n')

    for miRName in offsetperfectNames:
        miRNA, preName = miRName.split('\t')
        matureName = miRNA.split('_')[0] + '\t' + preName
        if matureName == 'hsa-miR-3168': continue
        heteroPos = miRName.split('_')[-1]
        matureSequence = totalSequenceDic[miRName]
        write_OffsetPerfect.write(matureName + '\t' + preName + '\t' + heteroPos + '\t' + matureSequence + '\t' + str(totalReadDic[miRName]) + '\t' + str(totalReadDic[matureName]) + '\t' + str(totalMappedReads) + '\n')

    for chunk in new_chunks:
        write_Expression.write(chunk)

    for miRName in totalReadMonoUDic.keys():
        miRNA, preName = miRName.split('\t')
        matureName = miRNA.split('_')[0] + '\t' + preName
        if matureName == 'hsa-miR-3168': continue
        writer_monoU.write(miRName + '\t' + str(totalReadMonoUDic[miRName]) + '\t' + str(totalReadDic[matureName]) + '\t' + str(totalMappedReads) + '\n')

    dOffsetNames = totalReadDoubleOffsetDic.keys()
    for miRName in onsetperfectNames:
        matureName, preName = miRName.split('\t')
        if matureName == 'hsa-miR-3168': continue
        for fwIndex in range(-5,6):
            readcounts = []
            for rvIndex in range(-5,6):
                if not totalReadDoubleOffsetDic.has_key(matureName + '_' + str(fwIndex) + '_' + str(rvIndex) + '\t' + preName):
                    continue
                readcounts.append(totalReadDoubleOffsetDic[matureName + '_' + str(fwIndex) + '_' + str(rvIndex) + '\t' + preName])
            if len(readcounts) == 0:
                continue
            write_DoubleOffsetPerfect.write(miRName + '_' + str(fwIndex) + '\t' + preName + '\t' + '\t'.join(map(str,readcounts)) + '\n')


    write_OnsetPerfect.close()
    write_OffsetPerfect.close()
    write_DoubleOffsetPerfect.close()
    write_OffsetSNV.close()
    write_OnsetSNV.close()
    write_Expression.close()
    writer_monoU.close()

def main(inputD, outputD, window, matureFa, type):
    matureSeqDic = get_sequenceDictionary(matureFa)
    if type == 'Catholic':
        samples = sorted(filter(lambda x: x.startswith('Sample_'), os.listdir(inputD)))
    elif type == 'TCGA':
        samples = sorted(filter(lambda x: x.startswith('TCGA'), os.listdir(inputD)))
    elif type == 'Tsinghua':
        samples = sorted(filter(lambda x: x.startswith('SRR'), os.listdir(inputD)))
    elif type == 'cellLine':
        samples = sorted(os.listdir(inputD))
    for sample in samples:
        print(sample)
        mrdF = inputD + sample + '/expression_analyses/' + sorted(os.listdir(inputD + sample + '/expression_analyses/'))[-1] + '/miRBase_new.mrd'
        totalSequenceDic, totalReadDic, totalMappedReads, new_chunks, totalReadMonoUDic, totalReadDoubleOffsetDic, totalMappedReadsDic = get_miRDeep2(mrdF, window)
        write_miRDeep2(sample, outputD, totalSequenceDic, totalReadDic, totalMappedReads, new_chunks, totalReadMonoUDic, totalReadDoubleOffsetDic, totalMappedReadsDic)
   
if __name__=='__main__':
    import sys
    import os
    if len(sys.argv) != 5:
        print("Usage : miRDeep2_analysis.py input_directory output_directory offset_window data_type")
        exit()
    inputD = sys.argv[1]
    outputD = sys.argv[2]
    window = int(sys.argv[3])
    type = sys.argv[4]
    matureFa = '/home/seokju/New/LiverCancer/data/mature.fa'
    if not os.path.exists(outputD): os.makedirs(outputD)
    main(inputD, outputD, window, matureFa,type)
