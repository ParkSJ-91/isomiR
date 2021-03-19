def parsing(mrdFile):
    mrdDic = dict()
    file = open(mrdFile); all = file.read(); file.close()
    chunks = all.split('>')[1:]
    for chunk in chunks:
        lines = chunk.split('\n')
        pre_miRNA = lines[0]
        total_read_count = lines[1].split(' ')[-1]
        seq_lines = filter(lambda x: 'seq_' in x, lines)
        if len(seq_lines) == 0: continue
        for seq_line in seq_lines:
            seq = seq_line.split(' ')[0]
            if not mrdDic.has_key(seq): mrdDic[seq] = []
            was_line = [pre_miRNA, total_read_count] 
            mrdDic[seq].append(was_line)
    return mrdDic

def comparing(mrdDic):
    newDic = dict() 
    for seqID in mrdDic.keys():
        if len(mrdDic[seqID]) == 1: continue
        total = sum(map(lambda x: int(x[1]), mrdDic[seqID]))
        for was_line in mrdDic[seqID]:
            if not newDic.has_key(was_line[0]): newDic[was_line[0]] = []
            seq_count = int(seqID.split('x')[1])
            new_line = [seqID, seqID.split('x')[0] + 'x' + str(int(round(seq_count*int(was_line[1])/float(total))))]
            newDic[was_line[0]].append(new_line)
    return newDic

def revising(mrdFile, newDic):
    writer = open(mrdFile.strip('.mrd') + '_new.mrd','w')
    f = open(mrdFile); all = f.read(); f.close()
    chunks = all.split('>')[1:]
    for pre_miRNA in newDic.keys():
        beRevisedChunk = filter(lambda x: pre_miRNA in x, chunks)[0]
        num = chunks.index(beRevisedChunk)
        for new_line in newDic[pre_miRNA]:
            wasChunk = chunks[num]
            chunks[num] = chunks[num].replace(new_line[0], new_line[1])
    for revised_chunk in chunks:
        writer.write('>' + revised_chunk)
    writer.close()

def main(mrdFile):
    mrdDic = parsing(mrdFile)
    newDic = comparing(mrdDic)
    revising(mrdFile, newDic)

if __name__=='__main__':
    import sys
    import os
    mrdFile = sys.argv[1]
    main(mrdFile)
