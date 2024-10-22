import pandas as pd
import sys

aa_file=sys.argv[1]
nuc_file=sys.argv[2]
outname=sys.argv[3]


aa_df=pd.read_csv(aa_file, header=None)
aa_df = pd.DataFrame({'label':aa_df[0].iloc[::2].values, 'aseq':aa_df[0].iloc[1::2].values})
aa_df.label= aa_df.label.str.lstrip(' ')
aa_df.label = aa_df.label.str.replace('-', '_')
aa_df.label = aa_df.label.str.replace('.', '_')

nuc_df=pd.read_csv(nuc_file, header=None)
nuc_df = pd.DataFrame({'label':nuc_df[0].iloc[::2].values, 'nseq':nuc_df[0].iloc[1::2].values})
nuc_df.label= nuc_df.label.str.lstrip(' ')
aa_df.label = aa_df.label.str.replace('.', '_')
nuc_df.label = nuc_df.label.str.replace('-', '_')

seqs_df=aa_df.merge(nuc_df)

print (seqs_df)

gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }


def tripletify(label,nuc_seq):
    nuc_seq=nuc_seq.upper()
    triplets=[nuc_seq[i:i+3] for i in range(0, len(nuc_seq), 3)]

    if len(triplets[-1]) != 3:
        print(triplets)
        raise ValueError("The nucleotide sequence of ",label," can not be mapped to the amino acid alignment, as it is not divisable by 3")

    return triplets


def map_gaps(triplets,aseq):
    mapped=''
    aseq=[*aseq]

    j=0
    print(len(aseq), len(triplets))
    for i in range(len(aseq)):
        if aseq[i] == '-':  #if aa == '-' add '---'
            mapped+='---'
        else: #for every aa add next 3 nuc. 
            #check if triplet translation matches nuccode for aa
            codon=triplets[j]
            if 'N' in codon.upper():
                print('Undefined codons detected!',codon)
            else:
                try:
                    gencode[codon]==aseq[j]
                except:
                    raise ValueError('Missmatch between nucleotide triplet and amino acid detected!')
            mapped+=codon
            j+=1
    
    #check if anuc=3xaa. 
    try:
        len(mapped) == 3*len(aseq)
    except:
        raise ValueError("Alignment mapped to nucleotide sequence does not have the correct length")
    
    return mapped


seqs_df=seqs_df.assign(triplets='',mapped_nucs='')
triplets=[]
mapped_nucs=[]
for i in range (len(seqs_df)):
    triplet=tripletify( seqs_df.iloc[i]['label'], seqs_df.iloc[i]['nseq'])
    triplets.append(triplet)
    
    mapped_gap=map_gaps(triplet, seqs_df.iloc[i]['aseq'])
    mapped_nucs.append(mapped_gap)



seqs_df['triplets']=triplets
seqs_df['mapped_nucs']=mapped_nucs

print(seqs_df)


def write_fasta_from_df(df,labelcol,seqcol,outname):
    out_name = outname + '.anuc'
    out_df = df.filter(items=[labelcol, seqcol])
    out_df[labelcol] =  out_df[labelcol]+"\n"
    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')
    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ',"")
            m.write( row + '\n')



write_fasta_from_df(seqs_df, 'label', 'mapped_nucs',outname)





