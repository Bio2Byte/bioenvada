import matplotlib.pyplot as plt
import pandas as pd
import glob
import itertools
import sys
import re


b2bfile=sys.argv[1]# 'results/b2b_msa_results_example_filtered_AA_checked.json' #
clade_csv=sys.argv[4] #glob.glob('*dist_thr*.csv')[0]#sys.argv[2]# 'results/example_filtered_NT_checked_dist_thr_1.5.csv' #
env_info=sys.argv[5]

figwid=sys.argv[2]
#figwid is figwidth, in cm, optional
min_occ=sys.argv[3]
#occupancy filter, percent, optional

outname=b2bfile.split('.')[0]


env_info_df=pd.read_csv(env_info,sep='\t')

try:
    clade_df=pd.read_csv(clade_csv,sep='\t',usecols=['strain','category'])
except:
    clade_df=pd.read_csv(clade_csv,sep='\t',usecols=['index','clusters_globalMax']) ##this is not supposed to be hardcoded !
    clade_df.rename(columns={"index": "strain", "clusters_globalMax": "category"}, inplace=True)

sgroup = clade_df['category'].unique()

if 'Cool' in sgroup:
    colors = {'Cold':'blue','Cool':"green", 'Warm':'orange', "Hot":'red'}

    #colors = {"Pro":'grey', 'Cold':'blue','Cool':"green", 'Warm':'orange', "Hot":'red'}
else:
    colors={}
    for i in range(len(sgroup)):
        colors[sgroup[i]]=plt.cm.tab10(i)

def json_to_df(infile):
    df = pd.read_json(infile, orient='index')
    dfseq = df.loc['sequence',:]

    df = df.drop(labels='sequence',axis=0)
    df=df.dropna(axis=1, how='all')
    dfseq=dfseq.dropna(how='all')

    df['sequence'] = dfseq
    #print(df)
    return df

def find_strain_in_index(index_value, strain_list):
    # Create a function that matches strains from df2 to df1 partially
    for strain in strain_list:
        if re.search(strain, index_value):
            return strain
    return None

def df_to_csv(df,clade_df):
    #add clades to df
    df.reset_index(inplace=True)
    cat_species=clade_df['strain'].tolist()
    
    df['cat_species']=df['index'].apply(lambda x: ''.join([part for part in cat_species if part in x ]))

    
    #qc step
    qc_df = df.merge(clade_df, left_on='cat_species', right_on='strain', how='outer')  

    filtered_df = qc_df[qc_df['index'].notnull()]

    print(filtered_df)
    filtered_df2=filtered_df[~filtered_df['category'].notnull()]

    unmatched_cols=filtered_df2['index'].tolist()
    for elem in unmatched_cols:
        if 'Syn_' in elem: #or 'Cya_' in elem:
                print ("Missing temp for strain", elem)
 
    qc_df[['index','cat_species','category','strain']].to_csv('testerdf.csv', sep='\t')

    df = df.merge(clade_df, left_on='cat_species', right_on='strain', how='inner') 


    #add evironmental info, this has to be partial match!
    # Step 2: Apply the function to df1['index'] to find the corresponding strain
    strain_list = env_info_df['strain'].tolist()
    df['short_strain'] = df['index'].apply(lambda x: find_strain_in_index(x, strain_list))

    # Step 3: Merge df1 and df2 based on the 'strain' column
    df = pd.merge(df, env_info_df,left_on='short_strain', right_on='strain', how='left')

    statdfs=[]
    cols=['backbone', 'sidechain', 'ppII', 'coil', 'sheet', 'helix', 'earlyFolding', 'disoMine']
    for col in cols:
        dfcols=df.columns.values.tolist()
        if col not in dfcols:
            print('Property', col, 'not calculated')
            continue
            
        dfname=outname+'_'+col+'.tsv'
        outdf = df[col].apply(pd.Series)
        #outdf['category'] = outdf.index.str.split('_',n=1,).str[0]
        outdf['category'] = df['category']
        outdf['species'] = df['index']
        outdf['sequence'] = [''.join(map(str, l)) for l in df['sequence']]

        outdf.to_csv(dfname)

        #df for pca
        pca_df=outdf.drop(columns=['sequence'])
        
        pca_df['temp'] = df['temp'] #should make this generic for all env properties
        pca_df=pca_df.add_prefix('R')#, axis =1)
        pca_df.rename(columns={'Rspecies':'species','Rcategory': 'category', 'Rtemp':'temp'}, inplace=True)
        
        
        pca_df.dropna(axis=1, how='any', inplace=True)
        
        
        nrow,ncol=pca_df.shape
        #pca_df.insert(ncol-2,'','')
        
        pca_df=pca_df.transpose()

        pcaname=outname+'_'+col+'_pca.tsv'
        pca_df.to_csv(pcaname)

        #get stats df
        statdf=outdf.groupby('category').describe().stack()
        statname=outname+'_'+col+'_stats.tsv'
        statdf.to_csv(statname)

        statdfs.append((col,statname))
    
    return statdfs


def plots(df, statname, prop,stretch):
    fig, (ax1,ax2) = plt.subplots(2,sharex=True, gridspec_kw={'height_ratios': [10, 1]})
    fig.set_figwidth(stretch)
    fig.set_figheight(8)
    # Plot raw data for each group

    print(statname)
    
    for i in prop:
        try:
            ax1.plot(df[i,'50%'], '-',  markersize=3 ,color=colors[i], label=f'{i} median')
        except:
            print("No member of group ", i, " found!")
            continue
       # ax1.plot(df[i,'50%'], '-',  markersize=3 ,color=colors[i], label=f'{i} median')
        #ax1.plot(df[i,'mean'], '-',  markersize=2 ,color=colors[i],alpha=0.3, label=f'{i} mean')

        ax1.fill_between(range(0,rows), df[i,'25%'].tolist(), df[i, '75%'].tolist(), alpha=0.15, color=colors[i], label=f'{i} 1st-3rd quartiles')
        #ax1.fill_between(range(0,rows), df[i,'min'].tolist(), df[i,'max'].tolist(), alpha=0.05, color=colors[i], label=f'{i} outliers')
        
        ax2.bar(range(0,rows), df[i,'count'],  width=1, color=colors[i],alpha=0.15,)

    fig.subplots_adjust(top=0.9, hspace = 0.0001)
    ax1.set_title(f'Prediction of biophysical propensity for {statname} for {str(rows)} residues')
    ax2.set_xlabel('Residue position')
    
    ax2.set_ylabel('Occupancy')
    ax1.set_ylabel('Propensity')
    ax1.legend(bbox_to_anchor=(1.001,0.5), loc="center left")
    
    plt.tight_layout()

    ax1.margins(x=0)
    ax2.margins(x=0)
    ax1.set_xticks(ax1.get_xticks()[::50])


    props=''.join(prop)

    if len(prop) == len(sgroup):
        plt.savefig(f'{outname}_{statname}_all.pdf')
    else:
        plt.savefig(f'{outname}_{statname}_{props}.pdf')
    plt.close()
    #plt.show()

b2bfile_df=json_to_df(b2bfile)
b2bfile_df=json_to_df(b2bfile)

stat_csv= df_to_csv(b2bfile_df,clade_df)
#fillis = glob.glob("*_stats.tsv")

for tup in stat_csv:
    statname=tup[0]
    f=tup[1]

    df = pd.read_csv(f,index_col=[0,1])
    df=df.transpose()

    rows,col = df.shape
    if figwid == "":
        figwid = rows*0.2
    else:
        figwid=float(figwid)

    plots(df, statname, sgroup, figwid)

    #occupancy filter
    if min_occ != '':
        min_occ=float(min_occ)
        mof=min_occ/100
        idx =[]
        for elem in sgroup:
            idx.append(list(df[df[elem,'count'] < df[elem,'count'].max()*mof].index))
        idx= list(itertools.chain(*idx))
        sidx = list(set(idx))
        ffdf = df.drop(index=sidx)

        ffstatname=statname+'_occupancy_'+str(min_occ)
        plots(ffdf, ffstatname, sgroup)


        combis=list(itertools.combinations(sgroup, 2))
        for elem in combis:
            fdf = pd.concat([df[elem[0]],df[elem[1]]],axis=1,keys=[elem[0],elem[1]])
            idx =[]
            for cl in elem:
                idx.append(list(df[df[cl,'count'] < df[cl,'count'].max()*mof].index))
            idx= list(itertools.chain(*idx))
            sidx = list(set(idx))
            print(sidx)
            fdf = df.drop(index=sidx)
            fstatname=statname+'_filtered_'+str(min_occ)
            plots(fdf, fstatname, elem)