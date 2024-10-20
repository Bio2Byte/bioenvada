import json
from b2bTools.multipleSeq.Predictor import MineSuiteMSA
import sys

from Bio import SeqIO

sequences=sys.argv[1]
sequences_name=sys.argv[2]
predtypes=sys.argv[3].replace(' ','').split(',')


print(predtypes)
print ('import complete')



#check if all seqs are the same, add trailing -
seq_dict = {rec.id : str(rec.seq) for rec in SeqIO.parse(sequences, "fasta")}
sameseqs=list(set(seq_dict.values()))
print(sameseqs)
if len(sameseqs) ==1:
    print('this file contains only identical sequences')

    with open (sequences, 'w') as m:
            for key in seq_dict.keys():
                row =  key+ '\n'+seq_dict[key]+'-'
                m.write( '>'+row + '\n')   




msaSuite = MineSuiteMSA()
msaSuite.predictAndMapSeqsFromMSA(sequences, predTypes=predtypes)
print("predictions done")

predictions_single_seq = msaSuite.allAlignedPredictions
json.dump(predictions_single_seq, open(sequences_name+'_b2b.json', 'w'), indent=2)

predictions=msaSuite.getDistributions()
json.dump(predictions, open(sequences_name+'b2b_stats.json', 'w'), indent=2)
