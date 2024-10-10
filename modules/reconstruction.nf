process tempRec{
    tag "${phylogeneticTree.baseName}"
    publishDir "$params.outFolder/reconstructions", mode: "copy"

    input:
    path phylogeneticTree
    path tempFile

    output:
    path '*.csv', emit: recTemps
    path '*.png'

    script:
    """
    Rscript $projectDir/bin/ancestral_state_rec.R $phylogeneticTree $tempFile ${phylogeneticTree.baseName}_rec_temps
    """
}

process pic{
    tag "${phylogeneticTree.baseName}"
    publishDir "$params.outFolder/pic", mode: "copy"

    input:
    path phylogeneticTree
    path traitFilelist

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    for traitFile in $traitFilelist ; do
        Rscript $projectDir/bin/run_pic.R $phylogeneticTree \$traitFile 
        done
    """
}