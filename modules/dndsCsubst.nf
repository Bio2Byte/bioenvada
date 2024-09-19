process findRoot{
    tag "${outGroup}"

    publishDir "$params.outFolder/tree", mode: "copy"

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

    publishDir "$params.outFolder/csubst", mode: "copy"

    input:
    tuple val(id), path (multipleSequenceAlignmentNuc), path (rootedTree), path (iqtreefiles)

    output:
    path "csubst_*" , emit: csubstOut

    script:
    """
    outgroup=${params.outGroup}
    echo "1	\$outGroup" > foreground.txt
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

    publishDir "$params.outFolder/csubst", mode: "copy"
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

    publishDir "$params.outFolder/dnds", mode: "copy"

    input:
    path multipleSequenceAlignmentNuc
    path rootedTree
    each eteEvol 

    output:
    //path "ete_${multipleSequenceAlignmentNuc}_*_out.tar.gz" , emit: eteOut
    path "ete_out/*", emit: eteOut

    script:

    """

    ete3 evol -t $rootedTree --alg $multipleSequenceAlignmentNuc --models M7 M8  --leaves --tests M7 M8 -o ete_out/  --cpu 10 -i ete_out/mypic.png &| tee etelog.log

    python $projectDir/bin/EteEvol.py ${multipleSequenceAlignmentNuc.name} $rootedTree $multipleSequenceAlignmentNuc $eteEvol
    
    """
    //python $projectDir/bin/EteEvol.py ${multipleSequenceAlignmentNuc.name} $rootedTree $multipleSequenceAlignmentNuc $eteEvol
    //ete3 evol -t $rootedTree --alg $multipleSequenceAlignmentNuc --models  b_free b_neut --leaves --tests b_free,b_neut -o ete_out/  --cpu 10
      //tar -cvzhf ete_${multipleSequenceAlignmentNuc}_${eteEvol}_out.tar.gz pamlwd/ plots/
  



}
