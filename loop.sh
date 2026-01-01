#!/bin/bash
#$ -cwd
#$ error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=18:00:00,h_data=20G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load plink

source ~/.bashrc


list=$1

while read line; do
        #pref=${line::-7} #if starting from bcf
        #basename="${line##*/}"
	#echo "$basename"
	pref=${line::-11} #if starting from fastq, remember to un-comment the bwa mem and mpileup lines
        echo $pref
        #qsub proc1.sh ${pref}
        qsub makebcf.sh ${pref}
	#bash makebcf2.sh $pref
done < $list
