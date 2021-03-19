import sys, os
inputD = sys.argv[1]

writer = open(inputD + '/context_scores.txt','w')
#context_scores_82.txt

files = sorted(filter(lambda x: x.startswith('context_scores_'), os.listdir(inputD)))
#print files[0]
#print open(inputD + files[0])
for line in open(inputD + files[0]):
    #print files[0]
    writer.write(line)

for file in files[1:]:
    print file
    for line in open(inputD + file):
        if line.startswith('Gene'): continue
        writer.write(line)
writer.close()

