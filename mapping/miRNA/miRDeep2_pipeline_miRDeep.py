def mapper(inputF, inputD, outputD, genome_index, threads, name):
    outputDir = outputD + '/2.miRDeep2/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    if not os.path.exists(outputDir+name): os.makedirs(outputDir+name)
    os.chdir(outputDir + name)
    #link_command = 'ln -s ' + outputD + '/1.adapter_trimmed_fastq/' + name + '.fastq .'
    link_command = 'ln -s ' + inputD + name + '.fastq .'

    print(link_command)
    p=Popen(link_command, shell=True)
    p.wait()

    mapper_command = '/home/seokju/.anyenv/envs/plenv/shims/perl /home/seokju/New/LiverCancer/src/mapper.pl ' + outputDir + name + '/' + name + '.fastq -e -h -j -l 18 -m -s ' + outputDir + name + '/' + name + '_reads_collapsed.fa > ' + outputDir + name + '/mapper_o.log 2> ' + outputDir + name + '/mapper_e.log'
    print(mapper_command)
    p=Popen(mapper_command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    
def quantifier(inputF, inputD, outputD, name, species, precursorD, matureD):
    outputDir = outputD + '/2.miRDeep2/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    os.chdir(outputDir + name)
    quantifier_command = '/home/seokju/.anyenv/envs/plenv/shims/perl /home/seokju/New/LiverCancer/src/quantifier.pl -p ' + precursorD + ' -m ' + matureD + ' -r ' + outputDir + name + '/' + name + '_reads_collapsed.fa -t '+species+' -d -g 2 > ' + outputDir + name + '/quantifier_o.log 2> ' + outputDir + name + '/quantifier_e.log'
    print(quantifier_command)
    p=Popen(quantifier_command, shell=True)
    p.wait()

def main(inputF, inputD, outputD, genome_index, threads, name, species, precursorD, matureD):
    mapper(inputF, inputD, outputD, genome_index, threads, name)
    quantifier(inputF, inputD, outputD, name, species, precursorD, matureD)

if __name__ == '__main__':
    import os
    import sys
    import time
    from subprocess import *

    source = sys.argv[1]
    threads = '10'

    precursorD = '/home/seokju/New/LiverCancer/data/hairpin.fa'
    matureD = '/home/seokju/New/LiverCancer/data/mature.fa'

    if source == 'Catholic':
        species = 'hsa'
        inputD = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/1.adapter_trimmed_fastq/'
        outputD = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/'
        inputFs = sorted(filter(lambda x: '.fastq' in x, os.listdir(inputD)))
    elif source == 'TCGA':
        species = 'hsa'
        inputD = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/2.sickled_fastq/'
        outputD = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/'
        inputFs = sorted(filter(lambda x: '.fastq' in x, os.listdir(inputD)))
    elif source == 'Tsinghua':
        species = 'hsa'
        inputD = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/1.adapter_trimmed_fastq/'
        outputD = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/'
        inputFs = sorted(filter(lambda x: '.fastq' in x, os.listdir(inputD)))

    if species == 'mmu':
        genomeIndex = '/Data_Set/Genome/mouse/mm9/bowtie/mm9'
    elif species == 'hsa':
        genomeIndex = '/home/seokju/New/LiverCancer/data/bowtie/hg19'

    for inputF in inputFs:
        name = inputF.split('.fastq')[0]
        main(inputF, inputD, outputD, genomeIndex, threads, name, species, precursorD, matureD)
