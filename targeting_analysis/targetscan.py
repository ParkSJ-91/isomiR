####################################
###                              ###
###        made by seokju        ###
###                              ###
###          2016.09.07          ###
###                              ###
####################################

# All files required for running this scripts do not have header lines. 
from threading import Thread
from threading import Semaphore
sem = Semaphore(50)
#import threading
#sem_print = Semaphore(1)
def seperate_files(inputD, sequenceF, geneID_list, number_of_seperation,now):
    print "Seperating", sequenceF, " files..."
    file_name = sequenceF.split('.')[0]
    if not os.path.exists(inputD + '/temp/'): os.makedirs(inputD + '/temp/')
    os.chdir(inputD + '/temp/')
    seperation_point = []
    for i in xrange(number_of_seperation):
        tempLength = len(geneID_list) * (i+1) / number_of_seperation
        seperation_point.append(tempLength)
    #print len(geneID_list)
    #print seperation_point

    now_point = 0
    now_length = 1 
    temp_geneID = []
    hit = 1
    f = open(now + '/'  + sequenceF, 'r'); lines = f.readlines(); lines.sort(); f.close()
    #print len(lines)
    writer = open(file_name + '_1.txt', 'w')
    writer.write(lines[0])
    temp_geneID.append(lines[0].split('\t')[0])
    for line in lines[1:]:
        line = line.strip()
        #if line.startswith('Ensembl') or line.startswith('CDR1as'): continue
        geneID = line.split('\t')[0]
        if temp_geneID[-1] != geneID: 
            temp_geneID.append(geneID)
            now_length += 1
        
        if now_length - 2  == seperation_point[now_point]:
            #print len(temp_geneID), temp_geneID[-1], lines[hit].split('\t')[0]
            writer.close()
            now_point += 1
            writer = open(file_name + '_' + str(now_point+1) + '.txt', 'w')

        writer.write(line + '\n')
        hit +=1

    writer.close()
    print "done"
    return

def get_geneID(sequenceF):
    geneID_list = []
    f = open(sequenceF, 'r'); lines = f.readlines(); f.close()
    lines.sort()
    geneID_list.append(lines[0].split('\t')[0])
    for line in lines[1:]:
         geneID = line.split('\t')[0]
         if geneID_list[-1] != geneID: geneID_list.append(geneID)
    print "Total number of transcripts :", len(geneID_list)
    return geneID_list

def check_files(utrF, orfF, airF):
    print "Check required files..."
    geneID_utr = get_geneID(utrF)
    geneID_orf = get_geneID(orfF)
    geneID_air = get_geneID(airF)
    if geneID_utr != geneID_orf and geneID_utr != geneID_air and geneID_orf != geneID_air:
        return False
    else: return True

class Auto(Thread):
#class Auto:
    def bsub(self, command_line, job_name, bsub_host):
        #print command_line
        #p = Popen(['bsub', '-I', '-J', job_name, '-q', bsub_host, command_line], stdout=PIPE, stderr=PIPE)
        p = Popen(command_line, shell=True)
        p.wait()
        return

    #def bsub2(self, command_line, job_name, bsub_host, outputF):
    #    writer = open(outputF, 'w')
    #    p = Popen(['bsub', '-I', '-J', job_name, '-q', bsub_host, command_line], stdout=PIPE, stderr=PIPE)
    #    output = p.stdout.read().split('\n')
    #    #writer.write(output[output.index('>.')+1:])
    #    #dohun debug#
    #    #print p.stdout.read()#dohun added"
    #    print output[0]
    #    writer.write('\n'.join(output[1:])) #; err = p.stderr.read()#dohun deleted
    #    writer.close()
    #    p.wait()
    #def bsub3(self, command_line, job_name, bsub_host, outputF):
    #    writer = open(outputF, 'w')
    #    p = Popen(['bsub','-I','-J',job_name, '-q', bsub_host, command_line], stdout=PIPE, stderr=PIPE)
    #    output = p.stdout.read()
    #    print output.split('\n')[:3]
    #    writer.write('\n'.join(output.split('\n')[3:]))
    #    writer.close()
    #    p.wait()

class AutoTargetScan(Auto):
    def __init__(self, inputD, miR_famF, miR_seqF, seperated_utrF, seperated_orfF, bsub_host, job_name):
        Thread.__init__(self)
        self._inputD = inputD
        self._miR_famF = miR_famF
        self._miR_seqF = miR_seqF
        self._seperated_utrF = seperated_utrF
        self._seperated_orfF = seperated_orfF
        self._bsub_host = bsub_host
        self._job_name = job_name
    
    def inputD(self): return self._inputD
    def miR_famF(self): return self._inputD + '/' + self._miR_famF
    def miR_seqF(self): return self._inputD + '/' + self._miR_seqF
    def utrF(self): return self._seperated_utrF
    def orfF(self): return self._seperated_orfF
    def bsub_host(self): return self._bsub_host
    def job_name(self): return self._job_name
    def job_number(self): return self._job_name.split('_')[1]

    def targetscan_70(self): #miR_famF, utrF, outputD, bsub_host):
        #os.chdir(inputD())
        command_line = 'perl /home/seokju/bin/targetscan_70.pl ' + self.miR_famF() + ' ' \
                    + self.utrF() + ' ' + self.inputD() + '/temp/predicted_targets' + self.job_number() + '.txt'
        print command_line
        self.bsub(command_line, self.job_name() + '_target' , self.bsub_host())
        return
    def targetscan_70_BL_bins(self): #utrF, outputD, bsub_host):
        #os.chdir(inputD)
        command_line = 'perl /home/seokju/bin/targetscan_70_BL_bins.pl ' + self.utrF() 
        print command_line
        self.bsub2(command_line, self.job_name() + '_bin' , self.bsub_host(), self.inputD() + '/temp/UTRs_median_BLs_bins' + self.job_number() + '.txt')
        return 
    def targetscan_70_BL_PCT(self): #miR_famF, outputD, bsub_host):
        #os.chdir(inputD)
        command_line = 'perl /home/seokju/bin/targetscan_70_BL_PCT.pl ' + self.miR_famF() + ' ' \
                    + self.inputD() + '/temp/predicted_targets' + self.job_number() + '.txt ' + self.inputD() + '/temp/UTRs_median_BLs_bins' \
                    + self.job_number() + '.txt'
        print command_line
        self.bsub3(command_line, self.job_name() + '_pct', self.bsub_host(), self.inputD() + '/temp/targetscan_70_output.BL_PCT' + self.job_number() + '.txt')
        return
    def targetscan_count_8mers(self): #miR_famF, orfF, outputD, bsub_host):
        #os.chdir(inputD)
#        command_line = '"perl /home/seokju/bin/targetscan_count_8mers.pl ' + self.miR_famF() + ' '\ #dohun deleted
#                    + self.orfF() + '"'# dohun deleted
        command_line = 'perl /home/seokju/bin/targetscan_count_8mers.pl ' + self.miR_famF() + ' '\
                    + self.orfF()# dohun edit
        print command_line
        self.bsub2(command_line, self.job_name() + '_8mer' , self.bsub_host(), self.inputD() + '/temp/ORF_8mer_counts' + self.job_number() + '.txt')
        return
    def targetscan_70_contexts_scores(self): #miR_seqF, utrF, outputD, job_name, bsub_host):
        #os.chdir(inputD)
        command_line = 'perl /home/seokju/bin/targetscan_70_context_scores.pl ' + self.miR_seqF()\
                    + ' ' + self.utrF() + ' ' + self.inputD() + '/temp/targetscan_70_output.BL_PCT' + self.job_number() \
                    + '.txt ' + '.'.join(self.orfF().split('.')[:-1]) + '.lengths.txt ' + self.inputD() + '/temp/ORF_8mer_counts' \
                    + self.job_number() + '.txt' + self.inputD() + '/result/context_scores' + self.job_number() + '.txt'
        print command_line
        self.bsub(command_line, self.job_name() + '_context', self.bsub_host())
        return
    def targetscan(self):
        cmd = '/home/seokju/Project/LiverCancer/src_rev1/ContextScore/targetscan.sh ' + self.miR_famF() + ' ' + self.miR_seqF() + ' ' + \
            self.utrF() + ' ' + self.orfF() + ' ' + self.job_number() + ' ' + self.inputD()
        print cmd
        self.bsub(cmd, self.job_name(), self.bsub_host())
    #def run(self):
    #    #os.chdir(self.inputD() + '/temp/')
    #    self.targetscan_70()
    #    self.targetscan_count_8mers()
    #    self.targetscan_70_BL_bins()
    #    self.targetscan_70_BL_PCT()
    #    #os.chdir(self.inputD())
    #    self.targetscan_70_contexts_scores()
    #    return
    def run(self):
        sem.acquire()
        self.targetscan()
        sem.release()

def main(miR_famF, miR_seqF, utrF, orfF, airF, number_of_seperation, inputD, bsub_host,now):
    if not check_files:
        print "Your geneID sets in ORF file and UTR file, AIR file are different"
        sys.exit()
    geneID_list = get_geneID(now+'/'+utrF)
    print len(geneID_list)
    #geneID_list2 = get_geneID(orfF)
    #if geneID_list != geneID_list2: print "fuck"; print geneID_list; print geneID_list2; exit()
    #seperate_files(inputD, utrF, geneID_list, number_of_seperation,now)
    #seperate_files(inputD, orfF, geneID_list, number_of_seperation,now)
    #exit()
    #time.sleep(500)
    for i in xrange(number_of_seperation):
        #if i!=0:continue
        job_name = 'targetscan_' + str(i+1)
        seperated_utrF = inputD + '/temp/' + utrF.split('/')[-1].split('.')[0] + '_' + str(i+1) + '.txt'
        seperated_orfF = inputD + '/temp/' + orfF.split('/')[-1].split('.')[0] + '_' + str(i+1) + '.txt'
        process = AutoTargetScan(inputD, miR_famF, miR_seqF, seperated_utrF, seperated_orfF, bsub_host, job_name)
        process.start()
    #for i in xrange(number_of_seperation)[:50]:
    #m = 0
    #for n in xrange(0,number_of_seperation,10):
    #    m += 10
    #    for i in xrange(n, m):
    #        if i == number_of_seperation: break
    #        job_name = 'targetscan_' + str(n) + '_' + str(i+1)
    #        #print job_name
    #        seperated_utrF = inputD + '/temp/' + utrF.split('/')[-1].split('.')[0] + '_' + str(i+1) + '.txt'
    #        seperated_orfF = inputD + '/temp/' + orfF.split('/')[-1].split('.')[0] + '_' + str(i+1) + '.txt'
    #
    #        process = AutoTargetScan(inputD, miR_famF, miR_seqF, seperated_utrF, seperated_orfF, bsub_host, job_name)
    #        process.start()
    #        
    #
    #    time.sleep(10)
    #    while commands.getoutput('bjobs -w').find('targetscan_' + str(n))>=0: time.sleep(10)
    #    time.sleep(10)
    #    #process.run()
    #return
    
if __name__=='__main__':
    import sys
    import os
    import time
    import commands
    from subprocess import *
    #from threading import Thread
    #import threading
    #Thread.__init__(self)
    if len(sys.argv) != 4:
        print 'usage : targetcan.py input_directory number_of_seperation bsub_host'
        sys.exit()
    print "Wait !! Do you change the file directory of TA_SPS, Agarwal parameters and RNA fold in targetscan_70_context_scores.pl script ??"
    time.sleep(5)
    inputD = sys.argv[1]
    number_of_seperation = int(sys.argv[2])
    bsub_host = sys.argv[3]
    miR_famF = filter(lambda x: 'miRNA' in x and 'family' in x, os.listdir(inputD))[0]
    miR_seqF = filter(lambda x: 'miRNA' in x and 'sequence' in x, os.listdir(inputD))[0]
    #utrF = filter(lambda x: 'UTR' in x and 'Sequences' in x, os.listdir(inputD))[0]
    #orfF = filter(lambda x: 'ORF' in x and 'Sequences' in x, os.listdir(inputD))[0]
    #airF = filter(lambda x: 'AIRs' in x, os.listdir(inputD))[0]
   
    now = '/home/seokju/Project/LiverCancer/src_rev1/ContextScore/'
    #utrF = filter(lambda x: 'UTR' in x and 'Sequences' in x, os.listdir(now))[0]
    #orfF = filter(lambda x: 'ORF' in x and 'Sequences' in x, os.listdir(now))[0]
    #airF = filter(lambda x: 'AIRs' in x, os.listdir(now))[0]
    utrF = 'UTR_Sequences.txt' #'UTR_Sequences_IsoSelected.txt'
    orfF = 'ORF_Sequences.txt' #'ORF_Sequences_IsoSelected.txt'
    airF = 'AIRs_Huh7.txt'
    # exist file check
    #if not 'Agarwal_2015_parameters.txt' in os.listdir(now):
    #    print 'Agarwal parameter file does not exist in input directory'
    #    sys.exit()
    #elif not 'TA_SPS_by_seed_region.txt' in os.listdir(now):
    #    print 'TA_SPS_by_seed_region file does not exist in input directory'
    #    sys.exit()
    #elif not 'RNAplfold_in_out' in os.listdir(now):
    #    print 'RNAplfold_in_out directory does not exist in input directory'
    #    sys.exit()
    #else:
    main(miR_famF, miR_seqF, utrF, orfF, airF, number_of_seperation, inputD, bsub_host,now)
    print "End"

