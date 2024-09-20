#!/usr/bin/env nextflow


log.info """\

Usage:

\$ nextflow run main.nf \\
    -profile standard,withdocker \\
    --targetSequences ../input_example.fasta \\
    --type 'aa' or 'nuc' \\
    --qc \\ 
    --clustering 0.85\\
    --relabel \\
    --alignSequences \\
    --efoldmine \\
    --disomine \\
    --agmata \\
    --fetchStructures \\
    --buildTreeEvo \\   
    --outGroup 'Species name to root your tree on' \\
    --csubst \\
    --branchIds '1,2,3'\\
    --eteEvol 'M7,M8' \\
    --selectedProteins 'your,proteins,as,str' \\
    --plotBiophysicalFeatures \\
    --buildLogo \\
    --plotTree

================================================================================
                                LIST OF PARAMETERS
================================================================================
                                GENERAL

Launch dir      : $launchDir
Project dir     : $projectDir
Execution time  : $params.executionTimestamp
Compressed file : ${params.compressedFile}.tar.gz
================================================================================
                                INPUT FILES

Input file (--targetSequences) : $params.targetSequences
Input file type (--type )      : $params.type  
================================================================================
                                FILTERING

Set minimal ooccupancy of position in MSA  (--qc 0.85)                            : $params.qc     
Clustering with CD-Hit (Siminarity percentage [0.85], word length 6)(--clustering): $params.clustering
Adapt labels to clustering (--relabel)                                            : $params.relabel
================================================================================
                                PREDICTORS

Align sequences for AA with Clustal, 
    for Nuc with MACSE(--alignSequences)        : $params.alignSequences

DynaMine                             : ALWAYS
DisoMine (--disomine)                : $params.disomine
EFoldMine (--efoldmine)              : $params.efoldmine
AgMata (--agmata)                    : $params.agmata
PSP                                  : NOT IMPLEMENTED

Fetch structures (--fetchStructures) : $params.fetchStructures

Phylo. Tree (--buildTreeEvo)         : $params.buildTreeEvo
Species name to root your tree on (--outGroup)   : $params.outGroup
Csubst (--csubst)                    : $params.csubst
CsubstSite (--branchIds)             : $params.branchIds
EteEvol (--eteEvol)                  : $params.eteEvol
================================================================================
                                OUTPUT FILES


Proteins to be highlighted in the plots (--selectedProteins)    : $params.selectedProteins
Plot B2btools (--plotBiophysicalFeatures)                       : $params.plotBiophysicalFeatures
Logo (--buildLogo)                : $params.buildLogo
Phylo. Tree plot (--plotTree)     : $params.plotTree
================================================================================
"""


// Modules
include {
    cdHitClustering;
    postClusteringLabels
    } from "${projectDir}/modules/cdhit"

include {
    predictBiophysicalFeatures;
    buildMultipleSequenceAlignmentAA;
    buildMultipleSequenceAlignmentNuc;
    takeMultipleSequenceAlignment;
    takeMultipleSequenceAlignmentNuc;
    buildPhylogeneticTree;
    buildPhylogeneticTreeEvol;
    buildLogo;
    treeToClade;
    mapMSA;
} from "${projectDir}/modules/multipleSequenceAlignment"

include {
    fetchEsmAtlasStructure;
    aaaToaaSeq
} from "${projectDir}/modules/structures"

include {
    plotBiophysicalFeaturesOverview as plotBiophysicalFeatures;
    plotPhylogeneticTree;
    plotEvoVsPhys;
    cladePlots;
} from "${projectDir}/modules/plots"

include {
    findRoot;
    runCsubst;
    runCsubstBranch;
    runEteEvol;
} from "${projectDir}/modules/dndsCsubst"

include { compressDirectory } from "${projectDir}/modules/utils"


workflow {
    
    params.compressedFile   = params.outFolder

    if (params.preprocessing == 'protein'){
        targetSequencesFile     = file(params.targetSequences)
        allSequences            = Channel.fromPath(params.targetSequences)
        
        //allSequences.view()
        seqsFiltered = allSequences
                .splitFasta( record: [header: true,  sequence: true ])
                .map { record -> [header: record.header.replaceAll("[^a-zA-Z0-9]", "_"),
                        sequence: record.sequence.replaceAll("\n","").replaceAll("[^ARNDCQEGHILKMFPOSUTWYVarndcqeghilkmfposutwvy-]", "X")] }

        seqsQC = seqsFiltered
            .branch{
                valid: it.sequence.count('X') / it.sequence.size() < 0.5
                invalid: it.sequence.count('X') / it.sequence.size() >= 0.5
            }.set { result}
        /*
                invalid:  it.header =~ /_Pro_/
                valid: true
        result.invalid.last().view{ "INVALID >${it.header}" }
        result.invalid.view{ "INVALID >${it.header}" }
        */
        sequencesSanitized = result.valid
        sequencesValid = sequencesSanitized.collectFile(name: "${targetSequencesFile.baseName}_filtered.fasta", newLine: true) {
            item -> '>' + item.header + '\n' + item.sequence + '\n'
        }
        sequencesRemoved = result.invalid.collectFile(name: "${targetSequencesFile.baseName}_sequences_ignored.fasta", newLine: true) {
            item -> '>' + item.header + '\n' + item.sequence + '\n'
        }
    }else if (params.preprocessing == 'proteome'){
       allSequences            = Channel.fromPath(params.targetSequences) //should be folder here!
        
        //allSequences.view()
        seqsFiltered = allSequences
                .splitFasta( record: [header: true,  sequence: true ])
                .map { record -> [header: record.header.replaceAll("[^a-zA-Z0-9]", "_"),
                        sequence: record.sequence.replaceAll("\n","").replaceAll("[^ARNDCQEGHILKMFPOSUTWYVarndcqeghilkmfposutwvy-]", "X"), orthologID: record.header[params.orthologIDLen] ]} //or record.header[0..params.orthologIDLen] 

        seqsQC = seqsFiltered
            .branch{
                valid: it.sequence.count('X') / it.sequence.size() < 0.5
                invalid: it.sequence.count('X') / it.sequence.size() >= 0.5
            }.set { result}

        sequencesSanitized = result.valid
        sequencesValid = sequencesSanitized.collectFile( newLine: true) {
                    item -> [ "${item.orthologID}_filtered.fasta", '>' + item.header + '\n' + item.sequence]
        }
        sequencesValid.view()
        sequencesRemoved = result.invalid.collectFile(newLine: true) {
                    item -> [ "${item.orthologID}_ignored.fasta", '>' + item.header + '\n' + item.sequence]
        }

    }


    if (params.clustering){
        cdHitClustering(sequencesValid, params.clustering)
        clusters =  cdHitClustering.out.clusters

        postClusteringLabels(sequencesValid,clusters, params.relabel)
        sequencesFiltered = postClusteringLabels.out.repSeqs
        representativeRepresented =  postClusteringLabels.out.representativeRepresented
    } else {
        sequencesFiltered = sequencesValid
        clusters = Channel.empty()
        representativeRepresented = Channel.empty()
    }

    if (params.alignSequences){
        if (params.type == 'nuc') {
            buildMultipleSequenceAlignmentNuc(sequencesFiltered)

            msaNuc = buildMultipleSequenceAlignmentNuc.out.msaNuc
            takeMultipleSequenceAlignmentNuc(msaNuc, params.buildTreeEvo, params.qc)
            multipleSequenceAlignmentNuc = takeMultipleSequenceAlignmentNuc.out.multipleSequenceAlignmentNuc

            msaAA = buildMultipleSequenceAlignmentNuc.out.msaAA
            takeMultipleSequenceAlignment(msaAA, params.buildTreeEvo, params.qc)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment
            
        } 
        if (params.type == 'aa') {
            multipleSequenceAlignmentNuc = Channel.empty()

            buildMultipleSequenceAlignmentAA(sequencesFiltered)
            msaAA = buildMultipleSequenceAlignmentAA.out.multipleSequenceAlignment
            takeMultipleSequenceAlignment(msaAA, params.buildTreeEvo, params.qc)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment

        } 
    } else {
        if (params.type == 'nuc') {
            multipleSequenceAlignment = Channel.empty()

            msaNuc =  sequencesFiltered
            takeMultipleSequenceAlignmentNuc(msaNuc, params.buildTreeEvo, params.qc)
            multipleSequenceAlignmentNuc = takeMultipleSequenceAlignmentNuc.out.multipleSequenceAlignmentNuc

        }
        else if (params.type == 'aa') {
            multipleSequenceAlignmentNuc = Channel.empty()
            msaAA = sequencesFiltered
            takeMultipleSequenceAlignment(msaAA, params.buildTreeEvo, params.qc)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment

        }
        else if (params.type == 'both') {
            multipleSequenceAlignmentNuc = Channel.empty()
            msaAA = sequencesFiltered
            takeMultipleSequenceAlignment(msaAA, params.buildTreeEvo, params.qc)
            multipleSequenceAlignment = takeMultipleSequenceAlignment.out.multipleSequenceAlignment

            mapMSA(multipleSequenceAlignment,params.nucToMap)
            multipleSequenceAlignmentNuc=mapMSA.out.msaNuc

        }else {println "Please select the input data type with --type 'nuc'"}
        
    }


    if (params.buildTree) {
        buildPhylogeneticTree(multipleSequenceAlignment)
        phylogeneticTree = buildPhylogeneticTree.out.tree
        iqtreeFiles = buildPhylogeneticTree.out.iqtreefiles

        if (params.plotTree) {
            plotPhylogeneticTree(phylogeneticTree, params.plotTree )
            plottedPhylogeneticTree = plotPhylogeneticTree.out.treePlot
        } else {
            //println "Skipping phylogenetic tree plot from MSA: ${targetSequencesFile}"

            plottedPhylogeneticTree = Channel.empty()
        }
    } else if (params.buildTreeEvo) {
        buildPhylogeneticTreeEvol(multipleSequenceAlignmentNuc)
        phylogeneticTree = buildPhylogeneticTreeEvol.out.tree
        iqtreeFiles = buildPhylogeneticTreeEvol.out.iqtreefiles

        if (params.plotTree) {
            plotPhylogeneticTree(phylogeneticTree, params.plotTree)
            plottedPhylogeneticTree = plotPhylogeneticTree.out.treePlot
        } else {
            //println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

            plottedPhylogeneticTree = Channel.empty()
        }
    } else {
        //println "Skipping Phylogenetic tree from MSA: ${targetSequencesFile}"
        //println "Skipping Phylogenetic tree plot from MSA: ${targetSequencesFile}"

        phylogeneticTree = Channel.empty()
        iqtreeFiles = Channel.empty()
        plottedPhylogeneticTree = Channel.empty()
    }

    if (params.buildLogo) {
        buildLogo(multipleSequenceAlignment)
        logo = buildLogo.out.logo
    } else {
        //println "Skipping Logo from MSA: ${targetSequencesFile}"

        logo = Channel.empty()
    }

    predictBiophysicalFeatures(multipleSequenceAlignment,params.dynamine, params.efoldmine, params.disomine, params.agmata)
    if (params.plotBiophysicalFeatures) {
        plotBiophysicalFeatures(
            multipleSequenceAlignment,
            params.dynamine,
            params.efoldmine,
            params.disomine,
            params.agmata,
            params.selectedProteins,
            predictBiophysicalFeatures.out.predictions,
            predictBiophysicalFeatures.out.stats,
            params.plotBiophysicalFeatures
        )
        plottedBiophysicalFeaturesInPNG = plotBiophysicalFeatures.out.plots
        
    } else {
        //println "Skipping Biophysical features plots from MSA: ${targetSequencesFile}"
        plottedBiophysicalFeaturesInPNG = Channel.empty()
    }

    if (params.cladePlots){ 
        if (params.groupInfo){
            //this needs testing
                      
            mapB2Bjson = predictBiophysicalFeatures
                .out.predictions
                .map(it -> tuple(it.baseName[params.orthologIDLen] , it))
                .combine(Channel.fromPath(params.groupInfo))


            mapB2Bjson.combine(Channel.fromPath(params.envInfoFile)).view()
            mapB2Bjson.combine(Channel.fromPath(params.envInfoFile)) | cladePlots

        }
        else{
            treeToClade(phylogeneticTree)
            cladeTab = treeToClade.out.cladeTab


            mapCladeTab = cladeTab.map(it -> tuple(it.baseName[params.orthologIDLen] , it))
            mapB2Bjson = predictBiophysicalFeatures.out.predictions.map(it -> tuple(it.baseName[params.orthologIDLen] , it))
            mapB2Bjson.join(mapCladeTab, remainder:true)
                .branch {
                    cladeNotFound:  it[1] == null
                    b2bNotFound:    it[2] == null
                    mapped:         true 
                }.set{matchCladesB2B}

            matchCladesB2B.mapped.combine(Channel.fromPath(params.envInfoFile))| cladePlots 
        }
        
    }

    if (params.fetchStructures) {
        // Sequences need to be shorter then 400 residues!)

        fetchEsmAtlasStructure(multipleSequenceAlignment.splitFasta(record: [header: true, seqString: true]).map { record -> [header: record.header, seqString: record.seqString.take(400)] })
        structures = fetchEsmAtlasStructure.out.esmStructures
    } else {
        //sequencesFiltered.view{ "Skipping fetching structures from EsmAtlas for sequences: "+ it}
        structures = Channel.empty()
    }


    if (params.csubst) {
        findRoot(phylogeneticTree, params.outGroup)
        rootedTree = findRoot.out.rootedTree

        mapAli = multipleSequenceAlignmentNuc.map(it -> tuple(it.baseName[params.orthologIDLen], it))
        mapTree = rootedTree.map(it -> tuple(it.baseName[params.orthologIDLen] , it))

        mapAli.join(mapTree, remainder:true)
            .branch {
                aliNotFound:        it[1] == null
                treeNotFound:       it[2] == null
                mapped:             true
            }.set{matchAliTree}
        
        mapIqtreeFiles =iqtreeFiles.map(it -> tuple(it[0].baseName[params.orthologIDLen] ,it))

        matchAliTree.mapped.join(mapIqtreeFiles, remainder:true)
            .branch {
                    aliNotFound:        it[1] == null
                    treeNotFound:       it[2] == null
                    mapped:         true 
                }.set{matchAliTreeIq}
        //matchAliTreeIq.mapped.view()

        if (params.branchIds) {
            matchAliTreeIq.mapped.combine(Channel.of(params.branchIds)) | runCsubstBranch
            //runCsubstBranch(multipleSequenceAlignmentNuc, rootedTree, csubstOutZip, params.branchIds)
            csubstOut = runCsubstBranch.out.csubstBranchOut
        } else{
            matchAliTreeIq.mapped | runCsubst
            csubstOut = runCsubst.out.csubstOut
        }

    } else{
        rootedTree = Channel.empty()
        csubstOut = Channel.empty()
    }


    if (params.eteEvol) {
        if (params.csubst == false){
            findRoot(phylogeneticTree, params.outGroup )
            rootedTree = findRoot.out.rootedTree

            mapTree = rootedTree.map(it -> tuple(it.baseName[params.orthologIDLen] , it))

            mapAli.join(mapTree, remainder:true)
                .branch {
                    aliNotFound:        it[1] == null
                    treeNotFound:       it[2] == null
                    mapped:             true
                }.set{matchAliTree}
        }

        //modelList = params.eteEvol?.split(',') as List
        //modelChannel = Channel.fromList(modelList)

        matchAliTree.mapped.combine(Channel.of(params.eteEvol)) | runEteEvol
        eteOutZip = runEteEvol.out.eteOut.toList()

    } else{
        rootedTree = Channel.empty()
        eteOutZip = Channel.empty()
    }

    //this step needs matching of msa and ete/csubst out
    //plotEvoVsPhys(eteOutZip,csubstBranchOutZip,predictBiophysicalFeatures.out.stats)

   

}

workflow.onComplete {
    println "Pipeline completed at               : $workflow.complete"
    println "Time to complete workflow execution : $workflow.duration"
    println "Execution status                    : ${workflow.success ? 'Success' : 'Failed' }"
    println "Compressed file                     : $params.outFolder/${params.compressedFile}.tar.gz"
    println "${workflow}.commandLine"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    println "Details: \n ${workflow.errorReport}"
}
