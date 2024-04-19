import pandas as pd
import glob
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches

from matplotlib import pyplot as plt

files = glob.glob('*msa_combined_outputs.csv')
first = True

for f in files:
    minidf = pd.read_csv(f)
    fname = f.split('_msa')[0]
    minidf.loc[:, 'name'] = fname
    minidf.loc[:, 'gene_name'] = fname.split('_')[2]

    if first == False:
        df = pd.concat([df,minidf],axis=0)
    else:
        first = False
        df = minidf

print (df)
        
conditions = [(df['OCNany2spe'] > 0.5),(df['OCNany2dif'] > 0.5)]
choices = ['P(Convergence)>0.5', 'P(Divergence)>0.5']
df['csubst_cat'] = np.select(conditions, choices, default='P(Adaptation)<0.5')
choices = [1, -1]
df['csubst_cat_num'] = np.select(conditions, choices, default=0)

condition1 = (df['dnds'].values.astype(float) > 1)
#condition2= (df['P_dnds_larger_1'].values.astype(float) > 0.5)
#df['dnds_cat'] = np.select([condition1 & condition2,condition1 ], [2,1], default=0)
df['dnds_cat'] = np.select([condition1], ['dN/dS>1'], default='dN/dS<=1')



tools  = list(set(tool.split('_')[0] for tool in df.columns))
interesting_tools = ['backbone', 'sidechain', 'sheet' , 'helix', 'earlyFolding', 'disoMine', 'agmata']
mytools=[]
for elem in tools:
    if elem in interesting_tools:
        mytools.append(elem)


def get_n(colx, tool, ax):
    #get median from plot
    #ax = box_plot.axes
    lines = ax.get_lines()
    categories = ax.get_xticks()

    medians=[]
    for cat in categories: 
        # every 4th line at the interval of 6 is median line
        # 0 -> p25 1 -> p75 2 -> lower whisker 3 -> upper whisker 4 -> p50 5 -> upper extreme value
        y = lines[4+cat*6].get_ydata()[0]
        #print(y)
        medians.append(y)

    #medians = df.groupby([colx])[f'{tool}_median'].median().values
    print(medians)
    nobs = df[colx].value_counts().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n=" + i for i in nobs]
    
    pos = range(len(nobs))
    for tick,label in zip(pos,ax.get_xticklabels()):
        ax.text(pos[tick],
                medians[tick] - 0.025,
                nobs[tick],
                horizontalalignment='center',
                size='x-small')
        ax.text(pos[tick],
                medians[tick] + 0.005,
                round(medians[tick],3),
                horizontalalignment='center',
                size='x-small')

#plots collated
for tool in mytools:
    tool_columns =  ['name','msa_pos']+[col for col in df.columns if col.startswith(tool)] +['gene_name','gap_rate_all','csubst_cat_num','OCNany2spe','OCNany2dif','csubst_cat','dnds','dnds_cat']
    df_tool = df.filter(items=tool_columns)

    #plt.figure(figsize=(10, 7.5))
    f, axs = plt.subplots(1,2,figsize=(9,5),sharey=True, gridspec_kw=dict(width_ratios=[2,3]))
    sns.boxplot(df_tool, x='dnds_cat', y=f'{tool}_median',fliersize=1, ax=axs[0])#showcaps=False,
    get_n('dnds_cat', tool, axs[0])
    
    csubstplot=sns.boxplot(df_tool, x='csubst_cat', y=f'{tool}_median',fliersize=1,ax=axs[1])
    get_n('csubst_cat', tool, axs[1])

    f.tight_layout()
    axs[1].set_ylabel('')
    #axs[1].tick_params(labelleft=False) 
    #axs[1].set(ylabel=None)#, yticks=[])
    f.subplots_adjust(wspace=0.05)  
    plt.savefig(tool+'_collated.png', dpi=600)
    plt.close()

    sns.stripplot(df_tool, x='gene_name', y=f'{tool}_median', size=4, hue= 'csubst_cat' ,dodge=True, palette="deep",legend=True)
    #sns.stripplot(df_tool, x='gene_name', y=f'{tool}_median', size=4, hue='dnds_cat', dodge=True, palette="pastel", legend=True)

    dndsdf = df_tool[df_tool.dnds_cat != 'dN/dS<=1']
    sns.stripplot(dndsdf, x='gene_name', y=f'{tool}_median', size=4,dodge=True,  color='indianred',  label='dN/dS>1')

    # Create a Custom Legend Handler
    unique_labels = list(set(df_tool['csubst_cat'])) + ['dN/dS>1']
    palette = sns.color_palette("deep")
    legend_handles = [mpatches.Patch(color=palette[i], label=label) for i, label in enumerate(unique_labels)]
    # Create the Legend outside the plot on the left
    plt.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(0, 1.1), ncol=2,fontsize='x-small')



    plt.savefig(tool + '_per_prot.png', dpi=600)
    plt.close()


'''
#plots per protein
for tool in mytools:
    tool_columns =  ['name','msa_pos']+[col for col in df.columns if col.startswith(tool)] +['gene_name','gap_rate_all','csubst_cat_num','OCNany2spe','OCNany2dif','csubst_cat','dnds','dnds_cat']
    df_tool = df.filter(items=tool_columns)

    plt.figure(figsize=(10, 7.5))

    sns.boxplot(df_tool, y=f'{tool}_median', x='gene_name', showcaps=False, color='lightgrey', saturation=0.2, fliersize=1, linewidth=1)
    dndsdf = df_tool[df_tool.dnds_cat > 0]
    sns.stripplot(dndsdf, x='gene_name', y=f'{tool}_median', size=4,dodge=True,  color='green', legend=True, label='dN/dS>1')
    #plt.savefig(tool+'_dnds.png', dpi=600)
    #plt.close()

    csubstdf = df_tool[df_tool.csubst_cat != 'P(Adaptation)<0.5']
    ax=sns.stripplot(csubstdf, x='gene_name', y=f'{tool}_median', size=4,  hue= 'csubst_cat' ,dodge=True, palette="deep")
    #plt.savefig(tool+'_csubst.png', dpi=600)

    # Create a Custom Legend Handler
    unique_labels = list(set(csubstdf['csubst_cat'])) + ['dN/dS>1']
    palette = sns.color_palette("deep")
    legend_handles = [mpatches.Patch(color=palette[i], label=label) for i, label in enumerate(unique_labels)]
    # Create the Legend outside the plot on the left
    plt.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(0, 1.1), ncol=3)

    # Save the plot
    plt.savefig(tool + '.png', dpi=600)
    plt.close()

    #Write csv
    outdf = pd.concat([csubstdf,dndsdf],axis=0)
    print(outdf)
    outdf.sort_values(by=['msa_pos'])
    outdf.to_csv(tool+'_selected.csv')
'''

    



