from ete4 import Tree
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster,linkage
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt


mytree=sys.argv[1] #"results_testing/example_filtered_NT_checked.afasta.treefile" 

inputname=mytree.split('.')[0]

#read tree
t = Tree(open(mytree),parser=1)
#cophenic matrix of tree
data = t.cophenetic_matrix()

df=pd.DataFrame(data=data[0],     columns=data[1], index = data[1]) 
df.reset_index(inplace=True)
print (df)
#probably should have an alert if distance >10 //10x the other


a=np.array([np.array(xi) for xi in data[0]])
#distlist is a list of distances [AB,AC,AD,BC,BD,CD]
distlist=a[np.triu_indices(a.shape[0],k=1)]

Z=linkage(distlist, method='ward', metric='euclidean')
print('ward')
print(Z)
#google if so else had cluster dist > node dist

Z=linkage(distlist, method='single', metric='euclidean')
print('single')
print(Z)

#Z=ward(distlist) == Z=linkage(distlist, method='ward', metric='euclidean')
Z=linkage(distlist, method='median', metric='euclidean')
print('median')
print(Z)


coeffs={}

dfmax=round(df.max(numeric_only=True).max()*10)
dfmin=round(df.mask(df==0).min(numeric_only=True).min()*10)

print('Testing distance tresholds ',dfmin/10, ' to ',dfmax/10 )
for i in range(dfmin,dfmax,2):
    i=i/10

    df['clusters']=fcluster(Z, i, criterion='distance')
    clusterlist=list(set(df['clusters'].to_list()))
    if clusterlist != [1]:
        sil=silhouette_score(a,df.clusters )
        coeffs[i]=sil
    else:
        print('only one cluster at', i)
        break


#print (coeffs)

#get max score
ideal_split_c=max(list(coeffs.values()))
ideal_split=list(coeffs.keys())[list(coeffs.values()).index(ideal_split_c)]
print("Ideal treshold: ", ideal_split)

df['clusters']=fcluster(Z, ideal_split, criterion='distance')
df['clusters'] = 'clade_'+df['clusters'].astype(str)

#df['clusters_max4_cl']=fcluster(Z, 4, criterion='maxclust')


outname=inputname+'_dist_thr_'+str(ideal_split)+'.csv'
df.to_csv(outname, sep='\t')


plt.plot(list(coeffs.keys()), list(coeffs.values()))
plt.xlabel('Distance threshold')
plt.ylabel('Silhouette score')

plt.savefig(inputname+"_dist_thres.pdf")





