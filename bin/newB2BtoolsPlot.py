import matplotlib.pyplot as plt
import pandas as pd
import glob
import itertools
import sys


b2bfile=sys.argv[1]# 'results/b2b_msa_results_example_filtered_AA_checked.json' #
clade_csv=sys.argv[4] #glob.glob('*dist_thr*.csv')[0]#sys.argv[2]# 'results/example_filtered_NT_checked_dist_thr_1.5.csv' #

figwid=sys.argv[2]
#figwid is figwidth, in cm, optional
min_occ=sys.argv[3]
#occupancy filter, percent, optional

outname=b2bfile.split('.')[0]
clade_df=pd.read_csv(clade_csv,sep='\t',usecols=['species','clusters'])

sgroup = clade_df['clusters'].unique()

if 'Cool' in sgroup:
    colors = {"Pro":'grey', 'Cold':'blue','Cool':"green", 'Warm':'orange', "Hot":'red'}
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

def df_to_csv(df,clade_df):
    #add clades to df
    df.reset_index(inplace=True)
    cat_species=clade_df['species'].tolist()
    
    df['cat_species']=df['index'].apply(lambda x: ''.join([part for part in cat_species if part in x ]))

    
    #qc step
    qc_df = df.merge(clade_df, left_on='cat_species', right_on='species', how='outer')  

    filtered_df = qc_df[qc_df['index'].notnull()]
    filtered_df2=filtered_df[~filtered_df['clusters'].notnull()]

    unmatched_cols=filtered_df2['index'].tolist()
    for elem in unmatched_cols:
        if 'Syn_' in elem: #or 'Cya_' in elem:
                print ("Missing temp for strain", elem)
 
    qc_df[['index','cat_species','clusters','species']].to_csv('testerdf.csv', sep='\t')

    df = df.merge(clade_df, left_on='cat_species', right_on='species', how='inner') 


    statdfs=[]
    cols = df.columns.values.tolist()

    for col in cols:
        if 'sequence' in col:
            continue
        if 'clusters'in col:
            continue
        if 'species'in col:
            continue
        if 'index'in col:
            continue

        dfname=outname+'_'+col+'.csv'
        outdf = df[col].apply(pd.Series)
        #outdf['clusters'] = outdf.index.str.split('_',n=1,).str[0]
        outdf['sgroup'] = df['clusters']
        outdf['sequence'] = [''.join(map(str, l)) for l in df['sequence']]

        outdf.to_csv(dfname)

        statdf=outdf.groupby('sgroup').describe().stack()
        statname=outname+'_'+col+'_stats.csv'
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
        ax1.plot(df[i,'50%'], '-',  markersize=3 ,color=colors[i], label=f'{i} median')
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
    if len(props) == len(sgroup):
        plt.savefig(f'{outname}_{statname}_all.pdf')
    else:
        plt.savefig(f'{outname}_{statname}_{props}.pdf')
    plt.close()
    #plt.show()

b2bfile_df=json_to_df(b2bfile)
b2bfile_df=json_to_df(b2bfile)

stat_csv= df_to_csv(b2bfile_df,clade_df)
#fillis = glob.glob("*_stats.csv")

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