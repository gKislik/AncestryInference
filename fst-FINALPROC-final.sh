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
	#pref=${line::-7} #if starting from vcf
	pref=${line::-11} #if starting from fastq
	echo $pref
	#qsub proc1-pca-mergefirst-final.sh ${pref}
	#qsub newapproach-updated-final.sh $pref 
        qsub newapproach-mergefirst-abr-final.sh $pref
        #bash na-Ronly.sh $pref
done < $list
