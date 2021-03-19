
def sickle(inputF, inputD, outputD, name, quality_type='sanger'):
    quality = int(20); length = int(18)
    outputDir = outputD + '/2.sickled_fastq/'
    if not os.path.exists(outputDir + '/log/'): os.makedirs(outputDir + '/log/')
    sc_command = 'sickle se -x -f ' + inputD + inputF + ' -t ' + quality_type + ' -q ' + str(quality) + ' -l ' + str(length) + ' -o ' + outputDir + name + '.fastq 2> ' + outputDir + '/log/e.log'
    print(sc_command)

    p = Popen(sc_command, shell=True)
    p.wait()

def fastqc2(inputF, inputD, outputD, name):
    fastqc_cmd = '/home/seokju/New/build/FastQC/fastqc -o ' + outputD + '/1.adapter_trimmed_fastq/ -f fastq -t ' + threads + ' ' + outputD + '/1.adapter_trimmed_fastq/' + name + '.fastq'
    print(fastqc_cmd)

    p = Popen(fastqc_cmd, shell=True)
    p.wait()

def main(inputF, inputD, outputD, threads, name):
    sickle(inputF, inputD, outputD, name,'sanger')
    fastqc2(inputF, inputD, outputD, name)

if __name__ == '__main__':
    import os
    import sys
    import time
    from subprocess import *

    source = sys.argv[1]
    threads = '10'

    if source == 'TCGA':
        inputD = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/'

    inputFs = sorted(filter(lambda x: '.fastq' in x, os.listdir(inputD+'/1.adapter_trimmed_fastq/')))

    for inputF in inputFs:
        name = inputF.split('.fastq')[0]
        main(inputF, inputD+'/1.adapter_trimmed_fastq/', inputD, threads, name)
