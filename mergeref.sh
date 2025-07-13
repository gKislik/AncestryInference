#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=70G

source ~/.bashrc
source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load plink

list=$1

ref=/u/project/pellegrini/gkislik/all_geno/final/LABR-POOD-CKSP/ref-SRR8614044

while read line; do
	pref="${line%%_*}"
	echo $pref
	#plink -bfile $line -chr-set 38 -allow-extra-chr -snps-only -make-bed -set-missing-var-ids @_# -out $pref
	plink -bfile $ref -bmerge $line -chr-set 38 -allow-extra-chr -make-bed -out test
	plink -bfile $ref -exclude test-merge.missnp -chr-set 38 -allow-extra-chr -make-bed -out $ref-semi
	plink -bfile $line -exclude test-merge.missnp -chr-set 38 -allow-extra-chr -make-bed -out $pref-semi
	plink -bfile $ref-semi -bmerge $pref-semi -chr-set 38 -allow-extra-chr -make-bed -out ref-$pref
	ref=ref-$pref
done < $list
		
