#!/usr/bin/env bash

if [ $# != 6 ]
    then
        echo 'input1 : miRNA family information'\
             'input2 : miRNA sequence information'\
             'input3 : UTR sequence file'\
             'input4 : ORF sequence file'\
             'input5 : file number'
             'input6 : input directory'
        exit 1
fi
nowDir=$6
outputDir=$nowDir'/temp/'
logDir=$outputDir'/log/'
if [ ! -d "$logDir" ]
    then
        mkdir $logDir
fi
resultDir=$nowDir'/result/'
if [ ! -d "$resultDir" ]
    then
        mkdir $resultDir
fi

# find target
targetscan_70.pl $1 $3 $outputDir'predicted_targets_'$5'.txt' 2> $logDir'targetscan_70_err_'$5'.txt'

# calculate ORF length and number of 8mer
targetscan_count_8mers.pl $1 $4 1> $outputDir'ORF_8mer_counts_'$5'.txt' 2> $logDir'targetscan_count_8mers_err_'$5'.txt'

# calculate branch length (BL) and probability of conserbed targeting (PCT) of 3UTR sequence
targetscan_70_BL_bins.pl $3 1> $outputDir'UTRs_median_BLs_bins_'$5'.txt' 2> $logDir'targetscan_70_BL_bins_err_'$5'.txt'

targetscan_70_BL_PCT.pl $1 $outputDir'predicted_targets_'$5'.txt' $outputDir'UTRs_median_BLs_bins_'$5'.txt' 1> $outputDir'targetscan_70_output.BL_PCT_'$5'.txt' 2> $logDir'targetscan_70_BL_PCT_err_'$5'.txt'

# calculate context scores
temp=$4
echo ${temp/%.txt/.lengths.txt}
targetscan_70_context_scores.pl $2 $3 $outputDir'targetscan_70_output.BL_PCT_'$5'.txt' ${temp/%.txt/.lengths.txt} $outputDir'ORF_8mer_counts_'$5'.txt' $resultDir'/context_scores_'$5'.txt' 2> $logDir'targetscan_70_context_scores_err_'$5'.txt'
