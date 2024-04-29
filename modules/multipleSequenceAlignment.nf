process predictBiophysicalFeatures {
    publishDir "$params.outFolder/msa", mode: "copy"

    tag "${sequences.name}"
    errorStrategy 'finish'
    time '10h'
    input:
    path sequences

    output:
    path 'b2b_msa_results_*.json', emit: predictions
    path 'b2b_msa_stats_*.json', emit: stats

    // Dynamine runs always; psp has no msa function

    script:
    """
    #!/usr/local/bin/python

    import json
    from b2bTools.multipleSeq.Predictor import MineSuiteMSA
    
    print ('import complete')

    msaSuite = MineSuiteMSA()
    msaSuite.predictAndMapSeqsFromMSA('$sequences', predTypes = ("dynamine", ${params.efoldmine ? '"eFoldMine",' : ''} ${params.disomine ? '"disoMine",' : ''} ${params.agmata ? '"agmata",' : ''}))
    print("predictions done")

    predictions_single_seq = msaSuite.allAlignedPredictions
    json.dump(predictions_single_seq, open('b2b_msa_results_${sequences.baseName}.json', 'w'), indent=2)

    predictions=msaSuite.getDistributions()
    json.dump(predictions, open('b2b_msa_stats_${sequences.baseName}.json', 'w'), indent=2)
    """

    
}

process buildMultipleSequenceAlignmentNuc {

    publishDir "$params.outFolder/msa", mode: "copy"

    label  'bigboy'
    tag "${sequences.name}"

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

    publishDir "$params.outFolder/msa", mode: "copy"

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

    publishDir "$params.outFolder/msa", mode: "copy"
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

    publishDir "$params.outFolder/msa", mode: "copy"
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

    publishDir "$params.outFolder/tree", mode: "copy"

    label  'bigboy'
    tag "${multipleSequenceAlignment.name}"

    input:
    path multipleSequenceAlignment

    output:
    path "*.treefile", emit: tree
    path "${multipleSequenceAlignment}.*" , emit:iqtreefiles

    script:
    """
    iqtree -s $multipleSequenceAlignment  -nt AUTO -B 10000 -m LG+R10
    """
    //FastTree $multipleSequenceAlignment > ${multipleSequenceAlignment}.tree
    //add ultrafast bootstrap -B 10000
    //iqtree -s $multipleSequenceAlignment  -nt AUTO -B 10000 -m LG+R10

}

process buildPhylogeneticTreeEvol {

    publishDir "$params.outFolder/tree", mode: "copy"

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

    publishDir "$params.outFolder/msa", mode: "copy"
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
