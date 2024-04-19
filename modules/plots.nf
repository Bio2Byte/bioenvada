
process plotBiophysicalFeaturesOverview {
    tag "${msa.name}"
    debug true
    errorStrategy "ignore"
    //publishDir "results/${msa.name}", mode: 'symlink'
    input:
    path msa
    val efoldmine
    val disomine
    val agmata
    val selected_proteins
    path predictions
    path stats
    val split

    output:
    path "*.png", emit: plots

    script:
    """
    python3 $projectDir/bin/MsaPlotB2btoolsBar.py $msa $predictions "${efoldmine ? 'efoldmine,' : ''} ${disomine ? 'disomine,' : ''} ${agmata ? 'agmata,' : ''}" $selected_proteins $stats $split
    """
}


process plotPhylogeneticTree {
    tag "${tree.name}"
    debug true

    // /home/sophie/miniconda3/envs/ete/bin/python
    //conda '/Users/sophie/miniconda3/envs/ete3'

    input:
    path tree
    val plotTree

    output:
    path "*.png", emit: treePlot

    script:
    """

    if [ $plotTree = 'phylo' ]; then
        plotstyle=0
        fi
    
    if [ $plotTree = 'evo' ]; then
        plotstyle=1
        fi

    python $projectDir/bin/treePlot.py $tree \$plotstyle
    """

    //PATH=/home/sophie/miniconda3/envs/ete/bin/:$PATH
    // /home/sophie/miniconda3/envs/ete/bin/python
}


process plotEvoVsPhys {
    input:
    path etetar
    path branchtar
    path b2bjson

    output:
    path "*.csv"
    path "*.png"


    script:
    """

    for file in *.tar.gz; do
        tar -xvzf \$file
    done

    python3 $projectDir/bin/parseOuts.py $b2bjson
    """
}