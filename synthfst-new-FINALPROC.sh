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
	#pref=${line::-7} #if starting from vcf (likely Garrett's script)
	pref=${line::-4} #if starting from bed file/plink output (GK's new script) 
	echo $pref
	#qsub proc1-pca-synthvcf-fst.sh ${pref} #if starting from vcf
	qsub proc1-pca-newsynthvcf-fst.sh ${pref} #if starting from bed file/plink output (GK's new script)
done < $list
