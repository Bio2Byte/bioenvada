import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

inpath = '' #results/CK_00001561_null_2023_03_17_16_36_03/csubst_CK_00001561_null_filtered.fasta_new.msa_out/'
results =  sys.argv[1] #inpath + 'csubst_cb_2.tsv' 
wdata =  sys.argv[2] #inpath+'csubst_s.tsv'


interesting_cols = {'omegaCany2dif' :  'omegaCany2dif \n Rate of combinatorial substitutions from any ancestral codons to a different codon \n divergence',
                    'omegaCany2spe': 'omegaCany2spe \n Rate of combinatorial substitutions from any ancestral codons to the specific descendant codons \n convergence', 
                    'OCNany2spe':  'OCNany2spe \n Number of nonsynonymous combinatorial substitutions \n convergence'  }



def matrixify(df, cols, rows):
    matrix = np.zeros(( rows, cols))
    list_of_df = df.values.tolist()
    for elem in list_of_df:
        col,row,val = int(elem[0])-1, int(elem[1])-1, elem[2]
        matrix[row,col] = val

    return matrix


def heatmap(matrix, name, max, cols, rows):
    #correlation matrix
    """plt.xticks(rotation=90) 
    ax.set_xticks(range(cols))
    ax.set_yticks(range(rows))
    ax.set_xticklabels(range(cols))
    ax.set_yticklabels(range(rows))"""

    plt.figure()
    # create a mask to hide the upper triangle of the matrix
    mask = np.zeros_like(matrix, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    matrix[mask] = np.nan
    
    # create a heatmap using matshow
    plt.matshow(matrix, cmap='YlGnBu', vmin=0.05, vmax=max)
    cbar = plt.colorbar(aspect=40, fraction=0.02)
    cbar.ax.tick_params(labelsize=8)
    
    plt.title(interesting_cols[name], fontsize=6)
    plt.savefig(inpath+name, dpi=300)

def dnds_bar(data):
    data['dN/dS'] = data['N_sub'] / data['S_sub']
    w = data['dN/dS'].values.tolist()

    plt.figure()
    plt.bar(range(0,len(w)), w,  width=.8)
    plt.ylim(0,5)
    plt.xlim(0,len(w))
    plt.savefig(inpath+'cSubst_dN_over_dS.png', dpi=300)




#plots for different cols in csubst_cb_2.tsv
res_df = pd.read_csv(results,sep='\t')
for datcol in interesting_cols.keys():
    wad = res_df[['branch_id_1','branch_id_2',datcol]]
    max= wad.max()
    cols, rows =int(max[0]),int(max[1])

    data= matrixify (wad, cols, rows)

    if 'omega' in datcol:
        max =10
    else:
        max =1
    heatmap(data, datcol, max, cols, rows)


#dnds Bar plot
w_df = pd.read_csv(wdata, sep='\t')
dnds_bar(w_df)

    