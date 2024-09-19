from ete4 import Tree,NodeStyle#,TreeStyle #ImportError: cannot import name 'TreeStyle' from 'ete4' (/home/sheidig/miniconda3/envs/bioenvada_base/lib/python3.12/site-packages/ete4/__init__.py). 
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
    if 4 in clusterlist:
        clust4hres=i
    if len(clusterlist) > 2:
        sil=silhouette_score(a,df.clusters )
        coeffs[i]=sil
    else:
        sil=silhouette_score(a,df.clusters )
        coeffs[i]=sil
        print('only two cluster at', i)
        break


#print (coeffs)

#get max score
ideal_split_c=max(list(coeffs.values()))
print(ideal_split_c)
ideal_split=list(coeffs.keys())[list(coeffs.values()).index(ideal_split_c)]
print("Ideal treshold: ", ideal_split)


localmax=[]
scores=list(coeffs.values())
for i in range(stepsize,len(scores)-stepsize):
    coef=scores[i]
    if coef >  scores[i-stepsize] and coef > scores[i+stepsize]:
        localmax.append(list(coeffs.keys())[i])
        localmax_split=list(coeffs.keys())[i]
print(localmax)



print("max threshold for 4 clusters:",clust4hres)


df['clusters_globalMax']=fcluster(Z, ideal_split, criterion='distance')
df['clusters_globalMax'] = 'clade_'+df['clusters_globalMax'].astype(str)

df['clusters_4clades']=fcluster(Z, clust4hres, criterion='distance')
df['clusters_4clades'] = 'clade_'+df['clusters_4clades'].astype(str)

try:
    df['clusters_localMax']=fcluster(Z, localmax_split, criterion='distance')
    df['clusters_localMax'] = 'clade_'+df['clusters_localMax'].astype(str)
except:
    print("No local maxima found")

#df['clusters_max4_cl']=fcluster(Z, 4, criterion='maxclust')


outname=inputname+'_dist_thr_'+str(ideal_split)+'.csv'
df.to_csv(outname, sep='\t')


plt.plot(list(coeffs.keys()), list(coeffs.values()))
plt.xlabel('Distance threshold')
plt.ylabel('Silhouette score')

plt.savefig(inputname+"_dist_thres.pdf")



"""
clusteroptions=['clusters_4clades','clusters_localMax','clusters_globalMax']


colors=['LightCoral','LightGrey','PeachPuff','LightCyan','AntiqueWhite','LightSalmon','Khaki','PaleVioletRed','Lavender','PaleGreen','Beige','PowderBlue']

for option in clusteroptions:
    clade_df=df.filter(items=['index',option])
    cladel=clade_df.groupby(by=option)['index'].apply(list)
    claded=cladel.to_dict()

    colmap={}
    for i in range(len(claded.keys())):
        clade=list(claded.keys())[i]
        if i >= len(colors):
            if i == len(colors):
                x=0
            else:
                x=x+1
        elif i ==0:
            x=len(colors)-1
        else:
            x=i
        colmap[clade]=colors[x]

    #print(colmap)

    nt=t
    for i in range(len(claded.keys())):

        clade=list(claded.keys())[i]

        leaves=list(claded.values())[i]

        

        print(clade,colmap[clade])
        nst = NodeStyle()
        nst["bgcolor"] = colmap[clade]
        n = nt.common_ancestor(leaves)
        n.set_style(nst)

    ts = TreeStyle()
    ts.mode = "c"
    nt.show(tree_style=ts)

"""
