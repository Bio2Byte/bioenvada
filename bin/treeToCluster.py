from ete4 import Tree
from ete4 import smNodeStyle

#from ete4.smartview.ete.layouts import TreeStyle
#from ete4 import Tree,NodeStyle#,TreeStyle #ImportError: cannot import name 'TreeStyle' from 'ete4' (/home/sheidig/miniconda3/envs/bioenvada_base/lib/python3.12/site-packages/ete4/__init__.py). 

import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster,linkage
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt


mytree=sys.argv[1] #"results_testing/example_filtered_NT_checked.afasta.treefile" 

inputname=mytree.split('/')[-1].split('.')[0]

#read tree
t = Tree(open(mytree),parser=1)
#cophenic matrix of tree
data = t.cophenetic_matrix()

df=pd.DataFrame(data=data[0],     columns=data[1], index = data[1]) 
df.reset_index(inplace=True)
#print (df)
#probably should have an alert if distance >10 //10x the other


a=np.array([np.array(xi) for xi in data[0]])
#distlist is a list of distances [AB,AC,AD,BC,BD,CD]
distlist=a[np.triu_indices(a.shape[0],k=1)]

Z=linkage(distlist, method='ward', metric='euclidean')
#print('ward')
#print(Z)
#google if so else had cluster dist > node dist

Z=linkage(distlist, method='single', metric='euclidean')
#print('single')
#print(Z)

#Z=ward(distlist) == Z=linkage(distlist, method='ward', metric='euclidean')
Z=linkage(distlist, method='median', metric='euclidean')
#print('median')
#print(Z)


coeffs={}

dfmax=round(df.max(numeric_only=True).max()*100)
dfmin=round(df.mask(df==0).min(numeric_only=True).min()*100)


stepsize=3
print('Testing distance tresholds ',dfmin/100, ' to ',dfmax/100 )
for i in range(dfmin+1,dfmax,stepsize):#range(dfmin,dfmax,20): ##dfmin: ValueError: Number of labels is 8. Valid values are 2 to n_samples - 1 (inclusive)-->dfmin+1
    i=i/100

    df['clusters']=fcluster(Z, i, criterion='distance')
    clusterlist=list(set(df['clusters'].to_list()))
    print(i, clusterlist)

    if 5 in clusterlist:
        clust5hres=i

    if 4 in clusterlist:
        clust4hres=i
    
    if 3 in clusterlist:
        clust3hres=i
    
    if len(clusterlist) > 1:
        sil=silhouette_score(a,df.clusters )
        coeffs[i]=sil
    else:
        #sil=silhouette_score(a,df.clusters )
        #coeffs[i]=sil
        print('only 2 cluster before', i)
        break

print (coeffs)

#get max score

print("max threshold for 3 clusters:",clust3hres)
print("max threshold for 4 clusters:",clust4hres)
print("max threshold for 5 clusters:",clust5hres)

ideal_split_c=max(list(coeffs.values()))
print("global max thres: ",ideal_split_c)
sil=[k for k, v in coeffs.items() if v == ideal_split_c]
print(sil)
ideal_split=sil[-1]

localmax_split=''
scores=list(coeffs.values())
for i in range(stepsize,len(scores)-stepsize):
    coef=scores[i]
    if coef >  scores[i-stepsize] and coef > scores[i+stepsize]:
        localmax_split=list(coeffs.keys())[i]
print('local max sil:', localmax_split)
#sil=[k for k, v in coeffs.items() if v == t]

thresholds=[clust3hres ,clust4hres,clust5hres]
chosenSil=0
for t in thresholds:
    s=coeffs[t]
    print(s)
    if s >= chosenSil:
        chosenThres = t
        chosenSil = s
print('clusters_set_clades', chosenThres,chosenSil)



gmax='c_globalMax_thres'+str(ideal_split)
df[gmax]=fcluster(Z, ideal_split, criterion='distance')
df[gmax] = 'clade_'+df[gmax].astype(str)

try:
    lmax='c_localMax_thres'+str(localmax_split)
    df[lmax]=fcluster(Z, localmax_split, criterion='distance')
    df[lmax] = 'clade_'+df[lmax].astype(str)
except:
    print("No local maxima found")

df["clusters_set_clades"]=fcluster(Z, chosenThres, criterion='distance')
df['clusters_set_clades'] = 'clade_'+df['clusters_set_clades'].astype(str)

#df['clusters_max4_cl']=fcluster(Z, 4, criterion='maxclust')


outname=inputname+'_dist_thr_'+str(chosenThres)+'.csv'
df.to_csv(outname, sep='\t')


plt.plot(list(coeffs.keys()), list(coeffs.values()))
plt.xlabel('Distance threshold')
plt.ylabel('Silhouette score')

plt.savefig(inputname+"_dist_thres.pdf")

