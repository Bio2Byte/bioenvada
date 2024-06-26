process findRoot{
    tag "${outGroup}"

    input:
    path phylogeneticTree
    val outGroup

    output:
    path '*rooted.treefile', emit: rootedTree

    script:
    """
    python  $projectDir/bin/rootTree.py $phylogeneticTree $outGroup 
    """
}

process runCsubst{

    label  'tinysnail'
    tag "${multipleSequenceAlignmentNuc.name}"

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree
    path iqtreefiles
    val outGroup

    output:
    path "csubst_*" , emit: csubstOut

    script:
    """
   
    echo "1	$outGroup" > foreground.txt
    cat foreground.txt

    csubst analyze --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --iqtree_redo no \
        --foreground foreground.txt \
        --fg_exclude_wg no \
        --cutoff_stat 'OCNany2spe,2.0|omegaCany2spe,5.0' \
        --max_arity 10 \
        --exhaustive_until 3
    
    """
} // python $projectDir/bin/cSubstAnalyzer.py csubst_cb_2.tsv csubst_s.tsv
//  echo '1	Syn_KORDI_100.*' > foreground.txt
//    echo '2	Syn_TAK9802.*' >> foreground.txt
//    echo '3	Syn_RS9916.*' >> foreground.txt
//    echo '4	Syn_BOUM118.*' >> foreground.txt
//    echo '5	Syn_WH7.*' >> foreground.txt
 //   echo '4	Syn_WH8103.*' >> foreground.txt
  //  echo '6	Syn_A15.*' >> foreground.txt


process runCsubstBranch{
    label  'tinysnail'
    tag "${multipleSequenceAlignmentNuc.name}"

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree
    path csubstOutZip
    val branchIds

    output:
    path "csubst_*" , emit: csubstBranchOut

    script:
    """

    if [ '$branchIds' = 'all' ]; then
        echo "starting scan"

        for file in csubst_cb_[0-9].tsv ; do
            python3 $projectDir/bin/filterBranchpairs.py \$file
        done

        pairs=\$(sort selectedbranches.txt | uniq)

        for pair in \$pairs; do
            echo "Pair: \$pair"
            csubst site --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --branch_id \$pair
        done

    else
        csubst site --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --branch_id ${branchIds} 
    fi
    
    """
}


process runEteEvol{
    errorStrategy 'ignore'
    label  'tinysnail'
    tag "${multipleSequenceAlignmentNuc.name}"
    conda '/Users/sophie/miniconda3/envs/ete3'

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree
    each eteEvol 

    output:
    path "ete_${multipleSequenceAlignmentNuc}_*_out.tar.gz" , emit: eteOut

    script:

    """
    python $projectDir/bin/EteEvol.py ${multipleSequenceAlignmentNuc.name} $rootedTree $multipleSequenceAlignmentNuc $eteEvol
    
    tar -cvzhf ete_${multipleSequenceAlignmentNuc}_${eteEvol}_out.tar.gz pamlwd/ plots/
    """



}
