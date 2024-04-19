import sys
from Bio import SeqIO
import pandas as pd
import glob
import os
import time
import numpy as np
from statistics import mean

folder= "data/wholeGenomeSyn/"
outfolder = folder+'clusters/'

outfolder_sel = folder +'clusters_sel/'
cdhitfolder = folder + "cdhit_sel/"

globaltic = time.perf_counter()

def createFullGenomeDataFrame (genome_nuc, genome_info, outname):

    #parse fna file (nucleotide sequences)

    seq_dict = {rec.id : str(rec.seq) for rec in SeqIO.parse(genome_nuc, "fasta")}
    seq_df = pd.DataFrame.from_dict(seq_dict, orient =  'index' , columns=['seq'])

    seq_df['id'] = seq_df.index
    seq_df[['id_short','start-stop','direction']] = seq_df['id'].str.split(":", expand=True)
    seq_df[['direction','product']] = seq_df['direction'].str.split("|", expand=True)

    seq_df = seq_df.filter(items=['id_short','start-stop','direction', 'seq']).reset_index(drop=True)
    #print (seq_df)


    #IMPROVE:probably better to just keep the dictionary until the end and make it then?
    #parse gff file (all other info)
    info_df = pd.read_csv(genome_info, sep='\t',header=1)
    info_df = info_df.attributes.str.split(';')

    rows = len(info_df.index)
    attribute_df = pd.DataFrame()
    first =True

    for row in range(0,rows):
        myd = {}
        for elem in info_df.loc[row]:
            key,val = elem.split('=')
            myd[key] = val

        mys = pd.Series(myd)

        if attribute_df.empty:
            attribute_df = mys
            firstname = myd['ID']
        elif first == True:
            attribute_df = pd.merge(attribute_df.rename('ID'), mys.rename(myd['ID']), how='outer', left_index=True, right_index=True)
            first = False
        else:
            attribute_df= pd.merge(attribute_df, mys.rename(myd['ID']), how='outer', left_index=True, right_index=True)


    attribute_df= attribute_df.transpose().reset_index(drop=True).drop(0)
    attribute_df['id_short'] = attribute_df['ID']
    attribute_df= attribute_df.drop("ID", axis='columns')
    #print (attribute_df)


    #merge into full genome info df
    all_df = pd.merge(seq_df,attribute_df , on='id_short')
    #print (all_df)
    all_df.to_csv(outname)
    return all_df

def treatAllGenomes (folder):

    genome_nuc_list = glob.glob(folder+"*.fna")
    counter =0

    for genome_nuc in genome_nuc_list:
        genome_name = genome_nuc.split('/')[-1].split('.')[0]
        
        csv = time.perf_counter()

        genome_info = genome_nuc[:-3]+'gff'
        outname = folder + genome_name + '_full_annotation.csv'

        createFullGenomeDataFrame (genome_nuc, genome_info, outname)

        tic = time.perf_counter()
        counter +=1
        print (f"Parsed genome {counter} of {len(genome_nuc_list)}: {genome_name} in {tic-csv:0.4f} seconds")

def collectAllHomologs (folder, outfolder,outfolder_sel):

    all_genomes_list =  glob.glob(folder+"fullGenomes/*.csv")
    all_genomes_df = pd.DataFrame()

    for genome in all_genomes_list:
        genome_df = pd.read_csv(genome)

        if all_genomes_df.empty:
            all_genomes_df = genome_df
        else:
            all_genomes_df = pd.concat([all_genomes_df,genome_df ], axis =0, ignore_index=True)
    
    #for i, x in all_genomes_df.groupby(by = 'cluster_number'): x.to_csv(f'{outfolder}{i}.csv', index=False)

    #count occurences of the cluster number to get most common //unique genomes per cluster

    #get stats of sequence length per cluster
    all_genomes_df['seq_len'] = all_genomes_df['seq'].str.len()
    occurrence_df = pd.DataFrame()
    occurrence_df['occurence'] = all_genomes_df.groupby('cluster_number')['cyanorak_Role'].count()
    occurrence_df['nuc_seq_len_min'] = all_genomes_df.groupby('cluster_number')['seq_len'].min()
    occurrence_df['nuc_seq_len_max'] = all_genomes_df.groupby('cluster_number')['seq_len'].max()
    occurrence_df['nuc_seq_len_spread'] = occurrence_df['nuc_seq_len_max'] - occurrence_df['nuc_seq_len_min']
    occurrence_df['nuc_seq_len_mean'] = all_genomes_df.groupby('cluster_number')['seq_len'].mean()
    occurrence_df['nuc_seq_len_median'] = all_genomes_df.groupby('cluster_number')['seq_len'].median()
    occurrence_df['nuc_seq_len_std'] = all_genomes_df.groupby('cluster_number')['seq_len'].std()
    print (occurrence_df)

    #get list of all cluster number and their meaning
    all_functions_df = all_genomes_df.filter(['cluster_number','cyanorak_Role','cyanorak_Role_description','eggNOG','eggNOG_description','kegg','kegg_description','ontology_term_description','product','protein_domains','protein_domains_description','tIGR_Role','tIGR_Role_description'])
    print ('genes in all genomes:', len(all_functions_df.index))
    
    #get only one entry per cluster
    all_functions = all_functions_df.drop_duplicates('cluster_number').sort_values('cluster_number')
    all_functions = all_functions.merge(occurrence_df , on = 'cluster_number' )
    print ('unique cyanorak cluster numbers:', len(all_functions.index))

    outfilename = folder + 'all_cluster_desc.csv'
    all_functions.to_csv(outfilename, index=False)

    #remove hypothetical proteins and genes without protein_domains annotation
    print ('clusters annotated as "hypothetical protein":', all_functions['product'].str.count('hypothetical protein').sum())
    
    all_functions['protein_domains'].replace('',np.nan, inplace=True)
    all_functions.dropna(subset=['protein_domains'], inplace=True)

    all_functions=  all_functions[all_functions['product'].str.contains('hypothetical protein')==False]
    print ('cyanorak cluster with protein annot:', len(all_functions.index))

    outfilename = folder + 'all_protein_cluster_desc.csv'
    all_functions.to_csv(outfilename, index=False)

    # check for ideal occurence cutoff 
    cutoff =[25,40,50]
    for cut in cutoff:
        all_functions = all_functions[all_functions['occurence']>cut]
        outfilename = folder + f'common_{cut}_protein_cluster_desc.csv'
        all_functions.to_csv(outfilename, index=False)
        print (f'cyanorak clusters in at least {cut}:', len(all_functions.index))
    selc = all_functions.cluster_number.to_list()

    filtered_all_functions = all_genomes_df[all_functions_df.cluster_number.isin(selc)]
    print ('selected clusters:', len(filtered_all_functions.groupby(by = 'cluster_number')))

    for i, x in filtered_all_functions.groupby(by = 'cluster_number'): x.to_csv(f'{outfolder_sel}{i}.csv', index=False)

def fetchGenomeInfo (folder):

    toc_file = folder + 'CyanorakOrganismTable.csv'
    toc_df = pd.read_csv(toc_file , sep=',', header=0)

    print (toc_df)


    all_genomes_list = toc_df.Name.values.tolist()

    print (all_genomes_list)

    for genome in all_genomes_list:

        print ('fetching: ' , genome)
        genome_nuc_url = 'http://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/genes/fna/' + genome
        genome_info_url ='http://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/gff/' + genome


        os.system(f'curl {genome_nuc_url} -o {folder + genome}.fna')
        os.system(f'curl {genome_info_url} -o {folder + genome}.gff')

        size = os.stat(folder +genome+'.fna').st_size
        if size <=127:
            with open(f'{folder}missing_genomes.txt', 'a') as lost:
                lost.write(genome + ': ' + genome_nuc_url + '\n')

def homologToFasta (outfolder):
    all_cluster_list =  glob.glob(outfolder+"*.csv")

    for cluster  in all_cluster_list:
        cluster_df = pd.read_csv(cluster)

        seq_df = cluster_df.filter(['id_short','translation'])

        #prep and write outfile
        seq_df['id_short']= ">" + seq_df['id_short'] + '\n'
        rows = seq_df.to_string(header=False,index=False,index_names=False).split('\n')

        out_name = cluster[: -4]+'.faa'

        with open (out_name, 'w') as m:
            for row in rows:
                row = row.replace('\\n','\n').replace(' ','').replace('+',"*").replace('#',"*")
                m.write( row + '\n')

def cdhit(outfolder_sel, cdhitfolder):
    all_cluster_seqs_list =  glob.glob(outfolder_sel+"*.faa")

    print ('clustering: ', len(all_cluster_seqs_list), 'files')

    for ortho in all_cluster_seqs_list:

        cd_tic = time.perf_counter()
        name = ortho[:-4].split('/')[-1]
        outname = cdhitfolder + name

        """treshl = [0.8, 1]
        for tresh in treshl:
            if tresh <=0.8:
                word = 4
            elif tresh <=0.85:
                word =5
            elif tresh <=0.88:
                word =6
            elif tresh <=0.9:
                word =7
            else:
                word = 10
            
            cdhitc = f'cd-hit-est -i {ortho} -o {outname}_{tresh*100}_clustered  -c {tresh} -n {word} -d 40'
            os.system(cdhitc)
        """
        
        cdhitc = f'cd-hit -i {ortho} -o {outname}_{50}_clustered_AA  -c 0.5 -n 3 -d 40'
        os.system(cdhitc)
        cd_toc = time.perf_counter()
        print (f"ran CDhit for {name} in {cd_toc-cd_tic:0.4f} seconds")

def getSeqIdentity (cdhitfolder, folder):
    cluster_list =  glob.glob(cdhitfolder+"*50_clustered_AA.clstr")
    print (len(cluster_list))

    cluster_dic = {}
    perc_list = []

    for cluster in cluster_list:
        with open (cluster, 'r') as cf:
            lines = cf.readlines()
            cluster_name = cluster.split('/')[-1].split('_')[1]

            for l in lines:
                if l.startswith(">"):
                    if perc_list!= []:
                        cluster_dic[key] = perc_list 
                    key = 'CK_'+ cluster_name + '-' +l.split(' ')[-1].replace('\n','')
                    perc_list = []
                else:
                    perc = l.replace('%\n','').split(' ')[-1]
                    if perc == '*\n':
                        perc = 99.0
                    perc_list.append(float(perc))

    cluster_df = pd.DataFrame.from_dict(cluster_dic, orient='index')
    cluster_df = cluster_df.apply(pd.DataFrame.describe, axis=1)
    cluster_df.insert(0, 'cluster_number', cluster_df.index)

    cluster_df = cluster_df.sort_values('cluster_number')

    outfilename = folder + 'sequence_variation_per_cluster.csv'
    #cluster_df.to_csv(outfilename, index=False)

    print ('cdhit clusters in all ck_cluster:', len(cluster_df.index))
    #make global approximation


    cluster_df[['cluster_number','cdhit_clstr']] =  cluster_df.cluster_number.str.split("-", expand=True)

    idx = cluster_df.groupby('cluster_number')['count'].transform(max) == cluster_df['count']    
    l = cluster_df[idx].index.to_list()
    
    print ('max cluster:', len(l))
    #if count is the same, it does not deciede, both are max

    filterl = []
    for fl in l:
        filterl.append(fl.split('-')[0])
    
    filterl = list(set(filterl))
    result = []
    for elem in filterl:
        result.append(list(filter(lambda x: x.startswith(elem), l))[0])
    
    print ('max cluster:', len(result))

    cluster_df = cluster_df.filter(['cluster_number','cdhit_clstr','count','mean'])
    cluster_df.loc[~cluster_df.index.isin(l), 'mean'] = 45.0
    print ('replaced perc: ', cluster_df['mean'].value_counts()[45])


    cluster_df = cluster_df.sort_values('cluster_number')
    outfilename = folder + 'sequence_variation_per_clusterONE.csv'
    #cluster_df.to_csv(outfilename, index=False)

    #print(cluster_df)
    #now average percertage:  ( x *xperc + y*yperc + ** ) / ( x+y+**)
    cluster_df['weight'] = cluster_df['count'] *cluster_df['mean']

    gdf= cluster_df.groupby('cluster_number')[['weight','count']].sum()
    gdf['approx_total_similarity'] = gdf['weight'] / gdf['count']
    gdf.drop("weight", axis='columns', inplace=True)

    print(gdf)

    outfilename = folder + 'approx_total_similarity_per_cluster.csv'
    gdf.to_csv(outfilename)


    #merge with common 50 df
    c50 = folder + 'common_50_protein_cluster_desc.csv'
    common50_df = pd.read_csv(c50, header =0)

    n_common50_df = pd.merge(common50_df,gdf,on='cluster_number')

    print (n_common50_df)
    n_common50_df.to_csv(c50, index=False)


    


#fetchGenomeInfo (folder)
#treatAllGenomes(folder)
#collectAllHomologs(folder, outfolder, outfolder_sel)
#homologToFasta (outfolder_sel)
#cdhit(outfolder_sel, cdhitfolder)
getSeqIdentity (cdhitfolder, folder)


globalcsv = time.perf_counter()
print (f"collated all homologs  in {globalcsv-globaltic:0.4f} seconds")

