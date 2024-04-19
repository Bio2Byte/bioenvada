import pandas as pd 
import sys


branchtab = sys.argv[1]

df = pd.read_csv(branchtab, sep='\t')

properties=['omegaCany2spe','dNCany2spe','OCNCoD','OCNany2spe','OCSany2spe','omegaCany2dif','dNCany2dif','OCNany2dif','OCSany2dif']
numbranches = int(branchtab.split('.')[0].split('cb_')[1])

n=0
branchonly = []
for num in range(0,numbranches):
    n +=1
    properties.append("branch_id_"+str(n))
    branchonly.append("branch_id_"+str(n))


filterdf = df[((df.OCNany2spe >= 3) &(df.OCSany2spe >= 3))] #| ((df.OCNany2dif >= 3) &(df.OCSany2dif >= 3))]
filterdf =filterdf.nlargest(10,'OCNany2spe')


if filterdf.empty ==False:
    coldf = filterdf.filter(items=properties)
    print (coldf)

    branchonlydf = filterdf.filter(items=branchonly)
    x = branchonlydf.to_string(header=False,
                    index=False,
                    index_names=False).split('\n')
else:
    filterdf = df.loc[df.OCNany2spe.idxmax()]
    
    branchonlydf = filterdf.filter(items=branchonly)
    print(filterdf)
    
    x = [branchonlydf.to_string(header=False,index=False).replace('\n',',')]


vals = [','.join(ele.split()) for ele in x]
print(vals)

with open('selectedbranches.txt','a')as sel:
    for val in vals:
        sel.write(val +'\n')






