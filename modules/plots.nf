
process plotBiophysicalFeaturesOverview {
    
    publishDir "$params.outFolder/plots/b2b", mode: "copy"

    tag "${msa.name}"
    debug true
    
    input:
    path msa
    val dynamine
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
    python3 $projectDir/bin/MsaPlotB2btoolsBar.py $msa $predictions "${dynamine ? 'dynamine' : ' '} ${efoldmine ? ',efoldmine' : ' '} ${disomine ? ',disomine' : ' '} ${agmata ? ',agmata' : ' '}" $selected_proteins $stats $split
    """
}


process plotPhylogeneticTree {


    publishDir "$params.outFolder/plots/", mode: "copy"
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
    echo $plotTree
    if [ $plotTree == 'phylo' ]; then
        plotstyle=0
        fi
    
    if [ $plotTree == 'evo' ]; then
        plotstyle=1
        fi

    python $projectDir/bin/treePlot.py $tree \$plotstyle
    """

    //PATH=/home/sophie/miniconda3/envs/ete/bin/:$PATH
    // /home/sophie/miniconda3/envs/ete/bin/python
}


process plotEvoVsPhys {

    publishDir "$params.outFolder/plots/evo_v_b2b", mode: "copy"
    errorStrategy 'ignore'
    debug true
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


process cladePlots {

    publishDir "$params.outFolder/plots/clade_plots", mode: "copy"
    input:
    path b2bjson
    val b2bfigwidth
    val b2boccupancy
    path cladeTab
    path envInfoFile

    output:
    path '*.pdf', optional: true 
    path '*.png', optional: true 
    path 'b2b*.csv' , emit: b2bPerTool

    script:

    """
    python3 $projectDir/bin/newB2BtoolsPlot.py "$b2bjson" "$b2bfigwidth" "$b2boccupancy" $cladeTab $envInfoFile
    
    """




}