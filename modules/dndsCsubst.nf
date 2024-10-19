process findRoot{
    tag "${outGroup}"

    publishDir "$params.outFolder/${phylogeneticTree.baseName.substring(0, 11)}/tree", mode: "copy"

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

    publishDir "$params.outFolder/${multipleSequenceAlignmentNuc.baseName.substring(0, 11)}/csubst", mode: "copy"

    input:
    tuple val(id), path(multipleSequenceAlignmentNuc), path(rootedTree), path(iqtreefiles), path (foregroundFile)
    output:
    path "csubst_*" , emit: csubstOut

    script:
    """
    csubst analyze --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --iqtree_redo no \
        --foreground  ${foregroundFile} \
        --fg_exclude_wg no \
        --cutoff_stat 'OCNany2spe,2.0|omegaCany2spe,5.0' \
        --max_arity 10 \
        --exhaustive_until 3
    
    """
} 


process runCsubstBranch{

    publishDir "$params.outFolder/${multipleSequenceAlignmentNuc.baseName.substring(0, 11)}/csubst/$oid", mode: "copy"
    label  'tinysnail'
    tag "${oid}"

    input:

    tuple val(oid), path(multipleSequenceAlignmentNuc), path(rootedTree), path(iqtreefiles), val(branchIds), path (foregroundFile)

    output:
    path "csubst_*" , emit: csubstBranchOut

    script:
    """
    cat foreground.txt

    csubst analyze --alignment_file ${multipleSequenceAlignmentNuc}  --rooted_tree_file ${rootedTree} --iqtree_redo no \
        --foreground ${foregroundFile} \
        --fg_exclude_wg no \
        --cutoff_stat 'OCNany2spe,2.0|omegaCany2spe,5.0' \
        --max_arity 10 \
        --exhaustive_until 3

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

   // echo '1	Syn_RS9902' > foreground.txt
   // echo '1	Syn_TAK9802' >> foreground.txt
    //echo '1	Syn_A15_127' >> foreground.txt
}


process runEteEvol{

    tag "${oid}"

    publishDir "$params.outFolder/${multipleSequenceAlignmentNuc.baseName.substring(0, 11)}/", mode: "copy"

    input:
    tuple val(oid), path(multipleSequenceAlignmentNuc), path(rootedTree), val(evolModel)

    output:
    //path "ete_${multipleSequenceAlignmentNuc}_*_out.tar.gz" , emit: eteOut
    path "ete_out/*", emit: eteOut

    script:

    """
    mkdir ete_out
    ete3 evol -t $rootedTree --alg $multipleSequenceAlignmentNuc --models $evolModel  --leaves --tests $evolModel -o ete_out/ --cpu 10  >> ete_out/etelog_${oid}.log
    
    """
    //ete3 evol -t $rootedTree --alg $multipleSequenceAlignmentNuc --models $evolModel  --leaves --tests $evolModel -o ete_out/ --cpu 10  -i ete_out/mypic.png >> ete_out/etelog_${oid}.log
        // error: end of tree file.qt.qpa.xcb: could not connect to display 
        //qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
    
    
    //python $projectDir/bin/EteEvol.py ${multipleSequenceAlignmentNuc.name} $rootedTree $multipleSequenceAlignmentNuc $eteEvol
    //tar -cvzhf ete_${multipleSequenceAlignmentNuc}_${eteEvol}_out.tar.gz pamlwd/ plots/
}
