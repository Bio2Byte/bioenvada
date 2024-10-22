process tempRec{
    tag "${phylogeneticTree.baseName}"
    publishDir "$params.outFolder/${phylogeneticTree.baseName.substring(0, 11)}/reconstructions", mode: "copy"

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
    publishDir "$params.outFolder/${phylogeneticTree.baseName.substring(0, 11)}/pic", mode: "copy"

    input:

    tuple val(oid), path(traitFilelist), path(phylogeneticTree)
    //path phylogeneticTree
    //path traitFilelist

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