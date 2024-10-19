process predictBiophysicalFeatures {
    publishDir "$params.outFolder/${sequences.baseName.substring(0, 11)}/msa", mode: "copy"

    tag "${sequences.name}"
    errorStrategy 'finish'
    time '10h'
    input:
    path sequences

    val dynamine
    val efoldmine
    val disomine
    val agmata
    
    output:
    path '*b2b.json', emit: predictions
    path '*stats.json', emit: stats

    script:
    """
   
    python3 $projectDir/bin/runB2Btools.py $sequences $sequences.baseName  "${dynamine ? 'dynamine' : ' '} ${efoldmine ? ',efoldmine' : ' '} ${disomine ? ',disomine' : ' '} ${agmata ? ',agmata' : ' '}"
    """
    // echo "${dynamine ? 'dynamine' : ' '} ${efoldmine ? ',efoldmine' : ' '} ${disomine ? ',disomine' : ' '} ${agmata ? ',agmata' : ' '}"
    
}

process buildMultipleSequenceAlignmentNuc {

    publishDir "${params.outFolder}/${sequences.baseName.substring(0, 11)}/msa", mode: "copy"

    label  'bigboy'
    tag "${sequences.baseName}"

    input:
    path sequences

    output:
    path '*NT.afasta', emit: msaNuc
    path "*AA.afasta", emit: msaAA

    script:
    """
    macse -prog alignSequences -seq ${sequences} -out_NT ${sequences.baseName}.afasta -out_AA ${sequences.baseName}.aafasta 
    macse -prog exportAlignment -align ${sequences.baseName}.afasta -codonForFinalStop --- -codonForInternalStop NNN -codonForInternalFS NNN -codonForExternalFS --- -charForRemainingFS ---
    """
}

process buildMultipleSequenceAlignmentAA {

    publishDir "$params.outFolder/${sequences.baseName.substring(0, 11)}/msa", mode: "copy"

    label  'bigboy'
    tag "${sequences.name}"

    input:
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    script:
    """
    clustalo -i $sequences -o ${sequences}.msa
    """
}

process takeMultipleSequenceAlignment {

    publishDir "$params.outFolder/${sequences.baseName.substring(0, 11)}/msa", mode: "copy"
    tag "${sequences.name}"

    input:
    path sequences
    val buildTreeEvo 
    val qc

    output:
    path "${sequences.baseName}_checked*" , emit: multipleSequenceAlignment

    script:
    """
    python3 $projectDir/bin/MsaChecker.py $sequences $buildTreeEvo $qc
    """

}

process takeMultipleSequenceAlignmentNuc {

    publishDir "$params.outFolder/${sequences.baseName.substring(0, 11)}/msa", mode: "copy"
    tag "${sequences.name}"

    input:
    path sequences
    val buildTreeEvo 
    val qc

    output:
    path "${sequences.baseName}_checked*" , emit: multipleSequenceAlignmentNuc

    script:
    """
    python3 $projectDir/bin/MsaChecker.py $sequences $buildTreeEvo $qc
    """

}

process buildPhylogeneticTree {

    publishDir "$params.outFolder/${multipleSequenceAlignment.baseName.substring(0, 11)}/tree", mode: "copy"

    label  'bigboy'
    tag "${multipleSequenceAlignment.name}"

    input:
    path multipleSequenceAlignment

    output:
    path "*.treefile", emit: tree
    path "${multipleSequenceAlignment}.*" , emit:iqtreefiles

    script:
    """
    iqtree -s $multipleSequenceAlignment  -nt AUTO -B 10000 
    """
    //FastTree $multipleSequenceAlignment > ${multipleSequenceAlignment}.tree
    //add ultrafast bootstrap -B 10000
    //iqtree -s $multipleSequenceAlignment  -nt AUTO -B 10000 -m LG+R10

}

process buildPhylogeneticTreeEvol {

    publishDir "$params.outFolder/${multipleSequenceAlignmentNuc.baseName.substring(0, 11)}/tree", mode: "copy"

    label  'bigboy'
    tag "${multipleSequenceAlignmentNuc.name}"

    input:
    path multipleSequenceAlignmentNuc

    output:
    path "*.treefile", emit: tree
    path "${multipleSequenceAlignmentNuc}.*" , emit:iqtreefiles

    script:
    """
    iqtree -s $multipleSequenceAlignmentNuc -m ECMK07+F+R4 -B 1000 --seqtype CODON11 -nt AUTO --ancestral --rate 
    """
    //iqtree -s $multipleSequenceAlignmentNuc -m ECMK07+F+R4 -B 1000 --seqtype CODON11 --threads-max 1 -T AUTO --ancestral --rate 

}

process buildLogo {

    publishDir "$params.outFolder/${multipleSequenceAlignment.baseName.substring(0, 11)}/msa", mode: "copy"
    tag "${multipleSequenceAlignment.name}"

    input:
    path multipleSequenceAlignment

    output:
    path "*_logo.png", emit: logo

    script:
    """
    weblogo --sequence-type protein --title "MSA logo" --first-index 0 --size large --format png_print < $multipleSequenceAlignment > ${multipleSequenceAlignment.simpleName}_logo.png
    """
}


process treeToClade {

    publishDir "$params.outFolder/${tree.baseName.substring(0, 11)}/tree", mode: "copy"
    tag "${tree.baseName}"

    input:
    path tree

    output:

    path '*.csv' , emit: cladeTab
    path "*.pdf"

    script:
    """
    
    python3 $projectDir/bin/treeToCluster.py $tree
    
    """

}

process mapMSA {
    publishDir "$params.outFolder/${alignment.baseName.substring(0, 11)}/msa", mode: "copy"

    tag "${alignment.name}"

    input:
    path alignment
    path sequencesNuc

    output:
    path "*.anuc", emit: msaNuc

    script:
    """
    python3 $projectDir/bin/mapAAtoNUC.py $alignment $sequencesNuc $sequencesNuc.baseName

    """

    //macse -prog reportGapsAA2NT -align_AA $sequences -seq $sequencesNuc -out_NT ${sequencesNuc.baseName}.anuc
    
}