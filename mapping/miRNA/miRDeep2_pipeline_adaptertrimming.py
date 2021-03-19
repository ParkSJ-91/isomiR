def fastqc(inputF, inputD, outputD, threads, name):
    fastqcD = inputD + '/fastqc/' + name + '/'
    if not os.path.exists(inputD + '/fastqc/' + name + '/'): os.makedirs(inputD + '/fastqc/' + name + '/')
    fastqc_cmd = 'fastqc -o ' + fastqcD + ' -f fastq -t ' + threads + ' -q ' + inputD + inputF
    print(fastqc_cmd)
    p = Popen(fastqc_cmd, shell=True)
    p.wait()

def cut_adapt(inputF, inputD, outputD, threads, name, adaptorSeq):
    outputDir = outputD + '/1.adapter_trimmed_fastq/'
    if not os.path.exists(outputD + '/1.adapter_trimmed_fastq/log/'): os.makedirs(outputD + '/1.adapter_trimmed_fastq/log/')
    if '_' in adaptorSeq:
        adaptor_fw, adaptor_rv = adaptorSeq.split('_')
        cut_command = '/home/seokju/.local/bin/cutadapt --overlap=6 -f fastq -a ' + adaptor_rv + ' ' + inputD + inputF + ' | /home/seokju/.local/bin/cutadapt --overlap=6 -g ' + adaptor_fw + ' -m 18 -M 26 - > ' + outputDir + name + '.fastq 2> ' + outputDir + '/log/e.log'
    else:
        cut_command = '/home/seokju/.local/bin/cutadapt --overlap=6 -f fastq -a ' + adaptorSeq + ' -m 18 -M 26 -o ' + outputDir + name + '.fastq ' + inputD + inputF + ' 2> ' + outputDir + '/log/e.log'
    print(cut_command)
    p=Popen(cut_command, shell=True)
    p.wait()

def main(inputF, inputD, outputD, threads, name, adaptorSeq):
    fastqc(inputF, inputD, outputD, threads, name)
    cut_adapt(inputF, inputD, outputD, threads, name, adaptorSeq)

if __name__ == '__main__':
    import os
    import sys
    import time
    from subprocess import *

    source = sys.argv[1]
    threads = '10'
    adaptorSeq = 'TGGAATTCTCGGGTGCCAAGG'
    
    if source == 'Catholic':
        inputD = '/home/seokju/New/LiverCancer/Catholic/1.miRNA/'
    elif source == 'TCGA':
        inputD = '/home/seokju/New/LiverCancer/TCGA/1.miRNA/'
    elif source == 'Tsinghua':
        inputD = '/home/seokju/New/LiverCancer/Tsinghua/1.miRNA/'

    inputFs = sorted(filter(lambda x: '.fastq' in x, os.listdir(inputD + '/0.fastq/')))

    for inputF in inputFs:
        name = inputF.split('.fastq')[0]
        main(inputF, inputD+'/0.fastq/', inputD, threads, name, adaptorSeq)
