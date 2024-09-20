import json
from b2bTools.multipleSeq.Predictor import MineSuiteMSA
import sys


sequences=sys.argv[1]
sequences_name=sys.argv[2]
predtypes=sys.argv[3].replace(' ','').split(',')


print(predtypes)
print ('import complete')

msaSuite = MineSuiteMSA()
msaSuite.predictAndMapSeqsFromMSA(sequences, predTypes=predtypes)
print("predictions done")

predictions_single_seq = msaSuite.allAlignedPredictions
json.dump(predictions_single_seq, open(sequences_name+'_b2b.json', 'w'), indent=2)

predictions=msaSuite.getDistributions()
json.dump(predictions, open(sequences_name+'b2b_stats.json', 'w'), indent=2)
