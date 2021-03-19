####################
##   2016.08.11   ##
##   made by SJ   ##
##    miRDeep2    ##
####################

def unzip(inputF, inputD, name):
    unzip_cmd = 'gzip -d ' + inputD + inputF
    print unzip_cmd
    p = Popen(unzip_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def fastqc(inputF, inputD, outputD, threads, name):
    fastqcD = inputD + '/fastqc/' + name + '/'
    if not os.path.exists(inputD + '/fastqc/' + name + '/'): os.makedirs(inputD + '/fastqc/' + name + '/')
    fastqc_cmd = 'fastqc -o ' + fastqcD + ' -f fastq -t ' + threads + ' -q ' + inputD + inputF
    print fastqc_cmd
    p = Popen(fastqc_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def cut_adapt(inputF, inputD, outputD, threads, name, adaptorSeq):
    outputDir = outputD + '/1.adapter_trimmed_fastq/'
    if not os.path.exists(outputD + '/1.adapter_trimmed_fastq/log/'): os.makedirs(outputD + '/1.adapter_trimmed_fastq/log/')
    if '_' in adaptorSeq:
        adaptor_fw, adaptor_rv = adaptorSeq.split('_')
        cut_command = '/home/seokju/.local/bin/cutadapt --overlap=6 -f fastq -a ' + adaptor_rv + ' ' + inputD + inputF + ' | /home/seokju/.local/bin/cutadapt --overlap=6 -g ' + adaptor_fw + ' -m 18 -M 26 - > ' + outputDir + name + '.fastq 2> ' + outputDir + '/log/e.log'
    else:
        cut_command = '/home/seokju/.local/bin/cutadapt --overlap=6 -f fastq -a ' + adaptorSeq + ' -m 18 -M 26 -o ' + outputDir + name + '.fastq ' + inputD + inputF + ' 2> ' + outputDir + '/log/e.log'
    print cut_command
    p=Popen(cut_command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()

def sickle(inputF, inputD, outputD, name, quality_type='sanger'):
    quality = int(20); length = int(18)
    outputDir = outputD + '/2.sickled_fastq/'
    if not os.path.exists(outputDir + '/log/'): os.makedirs(outputDir + '/log/')
    sc_command = 'sickle se -x -f ' + inputD + inputF + ' -t ' + quality_type + ' -q ' + str(quality) + ' -l ' + str(length) + ' -o ' + outputDir + name + '.fastq 2> ' + outputDir + '/log/e.log'
    print sc_command

    p = Popen(sc_command, shell=True)
    p.wait()

def fastqc2(inputF, inputD, outputD, name):
    fastqc_cmd = '/home/seokju/New/build/FastQC/fastqc -o ' + outputD + '/1.adapter_trimmed_fastq/ -f fastq -t ' + threads + ' ' + outputD + '/1.adapter_trimmed_fastq/' + name + '.fastq'
    print fastqc_cmd

    p = Popen(fastqc_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.wait()

def mapper(inputF, inputD, outputD, genome_index, threads, name):
    outputDir = outputD + '/2.miRDeep2/'
    if not os.path.exists(outputDir+name): os.makedirs(outputDir+name)
    os.chdir(outputDir + name)
    link_command = 'ln -s ' + outputD + '/1.adapter_trimmed_fastq/' + name + '.fastq .'
    print link_command
    p=Popen(link_command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()

    mapper_command = '/home/seokju/.anyenv/envs/plenv/shims/perl /home/seokju/New/LiverCancer/src/mapper.pl ' + outputDir + name + '/' + name + '.fastq -e -h -j -l 18 -m -s ' + outputDir + name + '/' + name + '_reads_collapsed.fa > ' + outputDir + '/' + name + '/mapper_o.log 2> ' + outputDir + '/' + name + '/mapper_e.log'
    print mapper_command
    p=Popen(mapper_command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    
def quantifier(inputF, inputD, outputD, name, species, precursorD, matureD):
    outputDir = outputD + '/2.miRDeep2/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    os.chdir(outputDir + name)
    quantifier_command = '/home/seokju/.anyenv/envs/plenv/shims/perl /home/seokju/New/LiverCancer/src/quantifier.pl -p ' + precursorD + ' -m ' + matureD + ' -r ' + outputDir + name + '/' + name + '_reads_collapsed.fa -t '+species+' -d -g 2 > ' + outputDir + '/' + name + '/quantifier_o.log 2> ' + outputDir + '/' + name + '/quantifier_e.log'
    print quantifier_command
    p=Popen(quantifier_command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()

def main(inputF, inputD, outputD, genome_index, threads, name, adaptorSeq, species, precursorD, matureD):
    if inputF.endswith('gz'):
        unzip(inputF, inputD, name, bsub_host)
        inputF = inputF.split('.gz')[0]
    #fastqc(inputF, inputD, outputD, threads, name)
    #cut_adapt(inputF, inputD, outputD, threads, name, adaptorSeq)
    #sickle(inputF, inputD, outputD, name,'sanger')
    #fastqc2(inputF, inputD, outputD, name)
    mapper(inputF, inputD, outputD, genome_index, threads, name)
    quantifier(inputF, inputD, outputD, name, species, precursorD, matureD)

if __name__ == '__main__':
    import os
    import sys
    import time
    import commands
    from subprocess import *

    analysisName = sys.argv[1]
    threads = sys.argv[2]

    precursorD = '/home/seokju/New/LiverCancer/data/hairpin.fa'
    matureD = '/home/seokju/New/LiverCancer/data/mature.fa'
    if analysisName == 'Liver_mmu':
        species = 'mmu'
        #genomeIndex = '/Data_Set/Genome/mouse/mm9/bowtie/mm9'
        inputD = '/home/seokju/Project/LiverCancer/cellLine2/Mouse_Liver_GSE113598/0.fastq/'  
        outputD = '/home/seokju/Project/LiverCancer/cellLine2/Mouse_Liver_GSE113598/'
        adaptorSeq = 'TGGAATTCTCGGGTGCCAAGG'
        inputFs = sorted(filter(lambda x: '.fastq' in x, os.listdir(inputD)))
    elif analysisName == "Liver_sihnRNPC":
        species = 'hsa'
        inputD = '/home/seokju/New/LiverCancer/CellLine/KD_hnRNPC/1.fastq/'
        outputD = '/home/seokju/New/LiverCancer/CellLine/KD_hnRNPC/'
        adaptorSeq = 'TGGAATTCTCGGGTGCCAAGG'
        inputFs = sorted(filter(lambda x: x.endswith('.fastq'), os.listdir(inputD)))
    elif analysisName == "Liver_sihnRNPC_U2AF2":
        species = 'hsa'
        inputD = '/home/seokju/New/LiverCancer/CellLine/KD_hnRNPC_U2AF2/1.fastq/'
        outputD = '/home/seokju/New/LiverCancer/CellLine/KD_hnRNPC_U2AF2/'
        adaptorSeq = 'TGGAATTCTCGGGTGCCAAGG'
        inputFs = sorted(filter(lambda x: x.endswith('.fastq'), os.listdir(inputD)))


    if species == 'mmu':
        genomeIndex = '/Data_Set/Genome/mouse/mm9/bowtie/mm9'
    elif species == 'hsa':
        genomeIndex = '/home/seokju/New/LiverCancer/data/bowtie/hg19'

    for inputF in inputFs:
        #print inputF
        name = inputF.split('.fastq')[0]
        #if 'Hiseq' in name: continue
        main(inputF, inputD, outputD, genomeIndex, threads, name, adaptorSeq, species, precursorD, matureD)
