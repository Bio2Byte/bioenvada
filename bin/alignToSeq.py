import sys


oldfile = sys.argv[1]
newfile =oldfile.split('.')[0]+ '.fasta'
print (newfile)

nf = open(newfile,'w') 

with open (oldfile, 'r') as of:
    for line in of:
        if line.startswith('>') == False:
            line = line.replace('-','').replace('*','')
        nf.write(line)

nf.close()