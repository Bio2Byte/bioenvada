process fetchEsmAtlasStructure {

    publishDir "$params.outFolder/structures", mode: "copy"
    tag "${header}"
    errorStrategy 'ignore'
    debug true

    input:
    tuple val(header), val(sequence)

    output:
    path "${header}.pdb", emit: esmStructures

    script:
    """
    echo Folding sequence using ESM Atlas
    echo To fetch: ${sequence.replaceAll('-', '').replaceAll('_', '')}
    sleep \$((RANDOM % 30))

    curl -X POST  -k -connect-timeout 8 --data "${sequence.replaceAll('-', '').replaceAll('_', '')}" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${header}.pdb

    find . -type f -name "*.pdb" -size -2k -exec bash -c 'mv "\$1" "\${1%.pdb}".fail' - '{}' +
    """
}


process aaaToaaSeq {

    publishDir "$params.outFolder/msa", mode: "copy"
    tag "${multipleSequenceAlignment}"

    input:
    path multipleSequenceAlignment

    output:
    path "*.fasta", emit: aaSeq

    script:
    """
    python3  $projectDir/bin/alignToSeq.py $multipleSequenceAlignment
    """
}