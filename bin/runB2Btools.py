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
json.dump(predictions_single_seq, open('b2b_msa_results_'+sequences_name+'.json', 'w'), indent=2)

predictions=msaSuite.getDistributions()
json.dump(predictions, open('b2b_msa_stats_'+sequences_name+'.json', 'w'), indent=2)
