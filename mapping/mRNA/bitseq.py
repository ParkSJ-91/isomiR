def build_bowtie2Index(fastaF, outputD):
    buildDir = '/'.join(outputD.split('/')[:-2]) + '/bowtie2_index/'
    if not os.path.exists(buildDir): os.makedirs(buildDir)
    build_cmd = '/home/seokju/build/bowtie2/bowtie2-2.3.4.1/bowtie2-build -f ' + fastaF + ' ' + buildDir + fastaF.split('/')[-1].rstrip('.fa')
    print build_cmd
    p = Popen(build_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

    time.sleep(10)
    while commands.getoutput('bjobs -w').find('build' + fastaF.split('/')[-1])>=0: time.sleep(10)
    time.sleep(10)
    
    return buildDir + fastaF.split('/')[-1].rstrip('.fa')

def make_dir(outputD, name):
    if not os.path.exists(outputD + '/2.bitseq/'+name): os.makedirs(outputD + '/2.bitseq/'+name)

def bitseq(outputFs, bitseq_outputD, bowtie2_index, fastaF, name, threads):
    if not os.path.exists(bitseq_outputD): os.makedirs(bitseq_outputD)
    cmd = '/home/seokju/New/LiverCancer/src/git/mapping/mRNA/bitseq.sh ' + outputFs[0] + ' '+ outputFs[1] + ' ' + bitseq_outputD + ' ' + bowtie2_index + ' ' + fastaF + ' ' + threads + ' ' + name
    print cmd

    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def main(inputD, samples, outputD, threads, fastaF):
    bowtie2_index = '/home/seokju/Project/LiverCancer/Catholic/2.mRNA/1.Bitseq/2.bitseq/bowtie2_index/gencode_v19_3PseqUpdated'
    #bowtie2_index = '/home/seokju/Project/LiverCancer/Catholic/NovelLncRNA/index/gencode.v27lift37.transcripts.novel.lncRNA'

    for sample in samples:
        print sample
        inputFs = sorted(filter(lambda x: 'fastq' in x or 'fq' in x, os.listdir(inputD + sample)))
        inputFs = filter(lambda x: not x.endswith('s.trimmed.filtered.fastq'), inputFs)
        inputFs = map(lambda x: inputD + sample + '/' + x, inputFs)
        bitseq(outputFs, outputD + '/2.bitseq/' + sample + '/', bowtie2_index, fastaF, sample, threads)

if __name__ == '__main__':
    import os
    import sys
    import time
    import commands
    from subprocess import Popen, PIPE, STDOUT
    import multiprocessing as mp

    source = sys.argv[1]
    thread = '10'
    sonNumber = 2
    if source == "Catholic":
        inputD = "/home/seokju/New/LiverCancer/Catholic/2.mRNA/1.preprocess/sickle/"
        outputD = "/home/seokju/New/LiverCancer/Catholic/2.mRNA/"
    build_index = '/home/seokju/Data/bowtie_hg19/hg19' #'/Data_Set/Genome/human/hg19/bowtie/hg19'
    mmr_threshold = 10
    quality_threshold = 20
    length_threshold = 20
    samples = filter(lambda x: not x.endswith('txt'), os.listdir(inputD))
    #fastaF = '/home/seokju/Data/gencode.v27lift37.transcripts.novel.lncRNA.fa'
    fastaF = '/home/seokju/Data/gencode_v19_3PseqUpdated.fa'
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
        
