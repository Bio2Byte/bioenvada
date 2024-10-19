
import sys
from ete3 import Tree


#allow root to be group of seqs! > get common anc of cya &vulcanococcus as root
#ancestor = t.get_common_ancestor("E","D")
#t.set_outgroup(ancestor)


phylogeneticTree = sys.argv[1] 
outGroup =sys.argv[2]   #'Cya_NS01_5_2B_1_CK_Cya_NS01_01838_1666461_1666877_1_CK_00001561_null'



t= Tree(phylogeneticTree, format=1)


def set_outGroup(outGroup,t):
    #find full name node
    setancestor = ''
    for node in t.traverse():
        if outGroup in node.name:
            setancestor = node.name
    
    if setancestor == '':
        raise ValueError("Outgroup not found in tree")
    
    t.set_outgroup(setancestor)
    return (t) 

def find_common_anc(outGroupS,t):
    full_ogs=[]

    for og in outGroupS:
        for node in t.traverse():
            if og in node.name:
                full_ogs.append(node.name)
    
    
    if full_ogs == []:
        raise ValueError("Outgroups not found in tree")

    print("full names of outgroups",full_ogs)

    ancestor = t.get_common_ancestor(full_ogs)
    print("lca of outgroups:",ancestor)
    
    t.set_outgroup(ancestor.name)

    return (t) 




if outGroup != '':
    outGroupS=outGroup.split(',')
    if len(outGroupS)>1:
        print('find common ancestor of outgroups and root')
        find_common_anc(outGroupS,t)
    else:
        print('set outgroup as root')
        set_outGroup(outGroup,t)
else:
    print('no outgroup found, use midpoint rooting')
    ancestor = t.get_midpoint_outgroup()
    t.set_outgroup(ancestor)


file_extension = len(phylogeneticTree.split('.')[-1]) 
out_name = phylogeneticTree[: -file_extension]+ 'rooted.treefile'

t.write(format=1, outfile=out_name)
