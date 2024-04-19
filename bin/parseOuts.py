
import pandas as pd
import glob
import json
import numpy as np
from matplotlib import pyplot as plt
#import seaborn as sns
import sys



etem8 =  glob.glob('pamlwd/*M8*/rst')[0] #'pamlwd/M8.CK_00001561_mini_filtered_clustered_NT_checked.afasta/rst'
csubstbranchfiles =  glob.glob('csubst_branch*/csubst_site.branch_id*/csubst_site.tsv')
b2btoolsjson =  sys.argv[1]#'b2b_msa_stats_CK_00001561_mini_filtered_clustered_AA_checked.json'

dataname = b2btoolsjson.split('stats_')[1].split('_filtered')[0] 

def ete_to_df(etefile,dataname):
    with open(etefile,'r') as ete:
        extractStart=False
        getdata = False
        dataset=[]
        for line in ete:
            if 'Bayes Empirical Bayes (BEB) probabilities' in line:
                extractStart = True
            if extractStart:
                if line.startswith('   1 '):
                    getdata = True
                if 'Positively selected sites' in line:
                    getdata = False
                    extractStart = False
                    dataset.pop()
                if getdata:
                    line = line.replace('(',"").replace(')',"")
                    dataset.append(line.rsplit())
    
    df = pd.DataFrame(dataset)
    print (df)
    cols = df.columns[[0,1,-3,-1, -5]]
    df = df.filter(items=cols)

    labels = {cols[0]:'msa_pos', cols[1]:'outgroup_seq', cols[2]:'dnds', cols[3]:'SD_dnds',cols[4]:'P_dnds_larger_1'}
    df.rename(columns=labels,inplace=True)

    df.to_csv(dataname+'_ete_M8.csv')
    print (df)
    return (df)




def find_best_branch(csubstbranchfiles,dataname):
    allthebest=[]
    for pair in csubstbranchfiles:
        print (pair)
        branchtup = pair.split('/')[1].split('_id')[1].split(',')

        df = pd.read_csv(pair, delimiter='\t', usecols=['codon_site_alignment','OCNany2spe','OCNany2dif'])
        #ana_df = df[(df.OCNany2spe >= 0.5) | (df.OCNany2dif >= 0.5) ] 
        spe =  df[df.OCNany2spe >= 0.5]
        dif =  df[df.OCNany2dif >= 0.5]

        findbest =[branchtup, spe.OCNany2spe.mean(), spe.OCNany2spe.sum(), dif.OCNany2dif.mean() , dif.OCNany2dif.sum(), pair]

        allthebest.append(findbest)
    
    print(allthebest)

    if allthebest!=[]:
        bestdf = pd.DataFrame(allthebest, columns=['branchtup','P50+_OCNany2spe_mean', 'P50+_OCNany2spe_sum','P50+_OCNany2dif_mean', 'P50+_OCNany2dif_sum','filename'])
        bestdf =bestdf.fillna(0)

        print(bestdf)

        bestdf.to_csv(dataname+'_summary_of_branch_combination.csv')

        best_spe = bestdf.loc[bestdf['P50+_OCNany2spe_sum'].idxmax(),'filename']
        best_dif = bestdf.loc[bestdf['P50+_OCNany2dif_sum'].idxmax(),'filename']

        if best_spe == best_dif:
            bestbranch = best_spe
        else:
            if bestdf['P50+_OCNany2spe_sum'].max() > bestdf['P50+_OCNany2dif_sum'].max():
                bestbranch = best_spe
            else:
                bestbranch = best_dif
        
        bestbranchdf =  pd.read_csv(bestbranch, delimiter='\t', usecols=['codon_site_alignment','OCNany2spe','OCNany2dif','gap_rate_all'])
        #print  (bestbranchdf)
    else:
        bestbranchdf=pd.DataFrame()
        print("noone meat the 50% criterion")
    return (bestbranchdf)



def b2btools_json_to_df(b2btoolsjson, etedf,bestbranchdf):
    with open(b2btoolsjson) as data_file:    
        data = json.load(data_file).get('results') 
    
    df = pd.json_normalize(data, sep='_')
    df = df.apply(lambda col: col.explode())
    df = df.reset_index(drop=True)


    # Extract the unique tool names from column names
    tools  = list(set(tool.split('_')[0] for tool in df.columns))

    # Calculate the spread for each tool
    for tool in tools:
        tool_columns = [col for col in df.columns if col.startswith(tool)]
        df[f'{tool}_spread'] = df[tool_columns].apply(lambda row: row[tool + '_thirdQuartile'] - row[tool + '_firstQuartile'], axis=1)
        df[f'{tool}_spread']=df[f'{tool}_spread'].round(5)

        tool_columns = [col for col in df.columns if col.startswith(tool)]
        df_tool = df.filter(items=tool_columns)

        #print (df_tool)
        cooldf = pd.concat([etedf,bestbranchdf,df_tool ],axis=1)
        allplot(cooldf, tool)
        
    print (df)
    return df

def allplot(df, tool):
    
    #conditions = [(df['OCNany2spe'] > 0.5), (df['OCNany2dif'] > 0.5)]
    #choices = [1, -1]
    #df['csubst_cat'] = np.select(conditions, choices, default=0)
    #condition1 = (df['dnds'].values.astype(float) > 1)
    #condition2= (df['P_dnds_larger_1'].values.astype(float) > 0.5)
    #df['dnds_cat'] = np.select([condition1 & condition2,condition1 ], [2,1], default=0)

    keywords = ['OCNany2spe', 'OCNany2dif', 'median','dnds' ,'spread']#'Outlier',

    # Use regex to filter columns based on partial matches
    filtered_df = df.filter(regex='|'.join(keywords), axis=1)

    print (filtered_df)

    #plt.figure(figsize=(10, 7.5))
    #sns.pairplot(filtered_df)

    #plname = dataname+'_'+tool+'_pairplot.png'
    #plt.savefig(plname, dpi=600) 

    #print (plname +' ist done')



print (etem8)
etedf = ete_to_df(etem8, dataname)
bestbranchdf = find_best_branch(csubstbranchfiles,dataname)

b2bdf = b2btools_json_to_df(b2btoolsjson,etedf,bestbranchdf)

cooldf = pd.concat([etedf,bestbranchdf,b2bdf],axis=1)
print (cooldf)

cooldf.to_csv(dataname+'_msa_combined_outputs.csv')

#allplot(cooldf)