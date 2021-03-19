def unzip(inputF, name):
    unzip_cmd = 'gzip -d ' + inputF 
    print unzip_cmd
    p = Popen(unzip_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def make_dir(outputD, name):
    if not os.path.exists(outputD + '/1.preprocess/fastqc/'+name): os.makedirs(outputD + '/1.preprocess/fastqc/'+name)
    if not os.path.exists(outputD + '/1.preprocess/fastqc2/'+name): os.makedirs(outputD + '/1.preprocess/fastqc2/'+name)
    if not os.path.exists(outputD + '/1.preprocess/bowtie/'+name + '/log/'): os.makedirs(outputD + '/1.preprocess/bowtie/'+name+'/log/')
    if not os.path.exists(outputD + '/1.preprocess/mismatch/'+name): os.makedirs(outputD + '/1.preprocess/mismatch/'+name)
    if not os.path.exists(outputD + '/1.preprocess/seqtk/'+name+'/log/'): os.makedirs(outputD + '/1.preprocess/seqtk/'+name+'/log/')
    if not os.path.exists(outputD + '/1.preprocess/sickle/'+name+'/log/'): os.makedirs(outputD + '/1.preprocess/sickle/'+name+'/log/')

def preprocess_fastqc(inputFs, outputD, name, threads):
    fastqcD = outputD + '/1.preprocess/fastqc/' + name +'/'
    if len(inputFs) == 1: readtype = 'single'
    else: readtype = 'paired'
    fastqc_cmd = '/home/seokju/New/build/FastQC/fastqc -o ' + fastqcD + ' -f fastq -t ' + threads + ' -q '
    if readtype == 'paired': fastqc_cmd += inputFs[0] + ' ' + inputFs[1]
    else: fastqc_cmd += inputFs[0]
    print fastqc_cmd

    p = Popen(fastqc_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def preprocess_bowtie(inputFs, outputD, build_index, name, threads):
    bowtieD = outputD + '/1.preprocess/bowtie/' + name + '/'
    if len(inputFs) == 1: readtype = 'single'
    else: readtype = 'paired'
    bowtie_cmd = '/home/seokju/New/build/bowtie-1.2.2-linux-x86_64/bowtie -q --chunkmbs 1024 -p ' + threads + ' ' + build_index + ' '
    if readtype == 'paired': bowtie_cmd += '-1 ' + inputFs[0] + ' -2 ' + inputFs[1] + ' '
    else: bowtie_cmd += inputFs[0] + ' '
    bowtie_cmd += bowtieD + '/' + name + '.bw'
    print bowtie_cmd

    p=Popen(bowtie_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def preprocess_mismatchRate(inputFs, outputD, name):
    bowtieD = outputD + '/1.preprocess/bowtie/' + name +'/'
    mismatchD = outputD + '/1.preprocess/mismatch/' + name +'/'
    hitNum = 0; mismD = dict()
    f = open(bowtieD + '/' + name + '.bw'); lines = f.readlines(); f.close()
    for line in lines:
        line = line.split('\n')[0].split('\t')
        chr = line[2]; sense = line[1]; start = int(line[3])+1
        mism = line[-1]; mism = mism.split(','); mism=filter(lambda x: len(x)>1, mism)
        if len(mism)>0:
            for m in mism:
                pos, m = m.split(':'); postN,preN = m.split('>'); pos = int(pos)
                if not (mismD.has_key(pos)): mismD[pos]=[]
                mismD[pos].append((postN,preN))
        hitNum +=1
    mmrF = open(mismatchD + '/' + name + '.mmr','w')
    pos = mismD.keys(); pos.sort()
    for p in pos: mmrF.write(str(p) + '\t' + str(round(100*len(mismD[p])/(float)(hitNum),2)) + '\n')
    mmrF.close()

def preprocess_seqtk(inputFs, outputD, mmr_threshold, name):
    mmrF = outputD + '/1.preprocess/mismatch/' +name + '/' + name + '.mmr'
    seqtkD = outputD + '/1.preprocess/seqtk/' + name + '/'
    mmr_threshold = float(mmr_threshold)
    f = open(mmrF); lines = f.readlines(); f.close()
    readLength = len(lines)
    trimLeft = 0; trimRight = 0
    for l in xrange(readLength/2):
        line = lines[l].strip().split('\t')
        left = int(line[0]); lmmr = float(line[1])
        if lmmr >= mmr_threshold: trimLeft = left+1
    for r in reversed(xrange(readLength/2,readLength)):
        line = lines[r].strip().split('\t')
        right = int(line[0]); rmmr = float(line[1])
        if rmmr >= mmr_threshold: trimRight = readLength - right
    for fastqF in inputFs:
        if trimLeft == 0 and trimRight == 0:
            seqtk_cmd = 'ln -s ' + fastqF + ' ' + seqtkD + '/' + fastqF.split('/')[-1].rstrip('.fastq') + '.trimmed.fastq'
        else:
            seqtk_cmd = '/home/seokju/New/build/seqtk-1.0/seqtk trimfq'
            if trimLeft > 0: seqtk_cmd += ' -b ' + str(trimLeft)
            if trimRight > 0: seqtk_cmd += ' -e ' + str(trimRight)
            seqtk_cmd += ' ' + fastqF + ' > ' + seqtkD + '/' + fastqF.split('/')[-1].rstrip('.fastq') + '.trimmed.fastq'
        print seqtk_cmd
        p=Popen(seqtk_cmd, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.wait()

def preprocess_sickle(inputFs, outputD, quality_threshold, length_threshold, name):
    fastqcD = outputD + '/1.preprocess/fastqc/' + name + '/'
    seqtkD = outputD + '/1.preprocess/seqtk/' + name + '/'
    sickleD = outputD + '/1.preprocess/sickle/' + name + '/'
    if len(inputFs) == 2: readtype = 'paired'
    else: readtype = 'single'
    quality_type = 'sanger'
    if inputFs[0].endswith('fastq'):
        #if readtype == 'single':
        tempName1 = inputFs[0].split('/')[-1].split('.fastq')[0]
        if readtype == 'paired':
            tempName2 = inputFs[1].split('/')[-1].split('.fastq')[0]
    elif inputFs[0].endswith ('.fq'):
        tempName1 = inputFs[0].split('/')[-1].split('.fq')[0]
        if readtype == 'paired':
            tempName2 = inputFs[1].split('/')[-1].split('.fq')[0]

    sickle_cmd = '/share/apps/programs/sickle/1.2/sickle '
    if readtype == 'paired': sickle_cmd += 'pe -f ' + seqtkD + '/' + tempName1 +'.trimmed.fastq -r ' + seqtkD + '/' + tempName2 +'.trimmed.fastq '
    else: sickle_cmd += 'se -f ' + seqtkD + '/' + name + '.trimmed.fastq '
    sickle_cmd += '-t ' + quality_type + ' -q ' + str(quality_threshold) + ' -l ' + str(length_threshold)
    if readtype == 'paired': sickle_cmd += ' -o ' + sickleD + '/' + tempName1 +'.trimmed.filtered.fastq -p ' + sickleD + '/' + tempName2 + '.trimmed.filtered.fastq -s ' + sickleD + '/' + name + '_s.trimmed.filtered.fastq'
    else: sickle_cmd += ' -o ' + sickleD + '/' + name + '.trimmed.filtered.fastq'
    print sickle_cmd

    p=Popen(sickle_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

    if readtype == 'paired':
        outputFs = [sickleD + '/' + tempName1 + '.trimmed.filtered.fastq', sickleD + '/' + tempName2 + '.trimmed.filtered.fastq']
    else:
        outputFs = [sickleD + '/' + tempName1 + '.trimmed.filtered.fastq']
    return outputFs

def preprocess_fastqc2(outputFs, outputD, name, threads):
    fastqcD2 = outputD + '/1.preprocess/fastqc2/' + name + '/'
    if len(outputFs) == 1: readtype = 'single'
    else: readtype = 'paired'
    fastqc2_cmd =  '/home/seokju/New/build/FastQC/fastqc -o ' + fastqcD2 + ' -f fastq -t ' + threads + ' -q '
    if readtype == 'paired': fastqc2_cmd += outputFs[0] + ' '+ outputFs[1]
    else: fastqc2_cmd += outputFs[0]
    print fastqc2_cmd
    
    p=Popen(fastqc2_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def main(inputD, samples, outputD, threads, fastaF):

    for sample in samples:
        inputFs = sorted(filter(lambda x: 'fastq' in x or 'fq' in x, os.listdir(inputD + sample)))
        inputFs = map(lambda x: inputD + sample + '/' + x, inputFs)
        name = sample
        for inputF in inputFs:
            if inputF.endswith('gz'):
                unzip(inputF, name)
            else:
                pass
        inputFs = map(lambda x: inputD + sample + '/' + x, sorted(filter(lambda y: y.endswith('fastq') or y.endswith('fq'), os.listdir(inputD + sample))))
        make_dir(outputD, name)
        print inputFs
        preprocess_fastqc(inputFs, outputD, name, threads)
        preprocess_bowtie(inputFs, outputD, build_index, name, threads)
        preprocess_mismatchRate(inputFs, outputD, name)
        preprocess_seqtk(inputFs, outputD, mmr_threshold, name)
        outputFs = preprocess_sickle(inputFs, outputD, quality_threshold, length_threshold,  name)
        preprocess_fastqc2(outputFs, outputD, name,threads)


if __name__ == '__main__':
    import os
    import sys
    import time
    import commands
    from subprocess import Popen, PIPE, STDOUT
    import multiprocessing as mp

    source = sys.argv[1]
    thread = '10'
    sonNumber = 1
    if source == "Catholic":
        inputD = "/home/seokju/New/LiverCancer/Catholic/2.mRNA/0.fastq/"
        outputD = "/home/seokju/New/LiverCancer/Catholic/2.mRNA/"
    elif source == "Tsinghua":
        inputD = "/home/seokju/New/LiverCancer/Tsinghua/2.mRNA/0.fastq/"
        outputD = "/home/seokju/New/LiverCancer/Tsinghua/2.mRNA/"
    build_index = '/home/seokju/Data/bowtie_hg19/hg19' 
    mmr_threshold = 10
    quality_threshold = 20
    length_threshold = 20
    samples = filter(lambda x: not x.endswith('txt'), os.listdir(inputD))[:1]
    #print samples
    #exit()
    fastaF = '/home/seokju/Data/gencode.v27lift37.transcripts.novel.lncRNA.fa'
    #fastaF = '/home/seokju/Data/gencode_v19_3PseqUpdated.fa'
    lists = []
    for i in xrange(sonNumber):
        lists.append([])
    for i in xrange(len(samples)):
        lists[i%sonNumber].append(samples[i])

    threads = []
    for i in xrange(len(lists)):
        tempSamples = lists[i]
        threads.append(mp.Process(target=main,args=(inputD,tempSamples,outputD,thread,fastaF)))
        threads[-1].deamon=True
    for i in xrange(0,len(threads)):
        threads[i].start()
    for i in xrange(0,len(threads)):
        threads[i].join()
        
    #main(inputD, samples, outputD, bsub_host, threads, fastaF)
