#!/bin/bash
now=`date +"%s"`
data=example 
nextflow run main.nf \
    -profile standard,withdocker \
    --targetSequences "/home/sheidig/bioenvada/sequences/*.fasta" \
    --preprocessing "proteome" \
    --alignSequences \
    --buildTreeEvo \
    --efoldmine \
    --disomine \
    --plotBiophysicalFeatures \
    --fetchStructures false \
    --csubst \
    --branchIds '1,5' \
    --eteEvol 'M8' \
    --outGroup 'Cya_NS01' \
    --selectedProteins  'AncNode14,Syn_BIOS_U3' \
    --buildLogo \
    --plotTree 'evo' \
    --cladePlots ${pwd}/results/testerClades.tsv


#-resume \
#-profile standard, withdocker, withsingularity, withconda \
#General Data Preparation
#    --type aa \
#    --clustering 1 \
#    --relabel \
#    --alignSequences \
#    --buildTreeEvo \
#B2Btools selection
#    --efoldmine \
#    --disomine \
#    --agmata \
#Evol predictions
#    --csubst \
#    --branchIds '1,2,3' \
#    --eteEvol\
#    --outGroup 'Cya_NS01_5_2B_1' \
#Generate Plots
#    --plotBiophysicalFeatures \
#    --selectedProteins  'Syn_WH5701 ,Syn_RS9917' \
#    --buildLogo \
#    --plotTree \
#    --fetchStructures \