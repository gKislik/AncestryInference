#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=1:00:00,h_data=70G

source ~/.bashrc
source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load plink
module load R

N=$1
clust=$2
id=$3 #can be anything pretty much, just make sure you don't use duplicate ones between runs
ref=/u/project/pellegrini/gkislik/all_geno/final/LABR-POOD-CKSP/filt-ref-SRR8614076
#shuffle SNPs
shuf finSNPs.txt > shuffledfinSNPs.txt
nlines=$(wc -l shuffledfinSNPs.txt | cut -d' ' -f1)

#pick random lines from cluster
shuf $clust > shufclust.txt
head -n $N shufclust.txt > chosen.txt
awk '{print $1"\t"$2}' chosen.txt > samples.txt
awk '{print $3}' chosen.txt > breeds.txt

echo "$nlines"
echo "$N"
#make random cutoffs 
shuf -i 1-"${nlines}" -n "$N" | sort -n > cutoffs.txt

iter=1
while read cutoff ; do
	#echo "$iter"
	#echo "$cutoff"
	if [ $iter == 1 ]; then #for first entry
		#stuff=$(sed -n "${iter}p" < chosen.txt)
		echo $cutoff
		len=$(head -n ${cutoff} shuffledfinSNPs.txt | wc -l)
		breed=$(head -n 1 breeds.txt)
		prop="$((1000 * len / nlines))"
		echo "$prop"
		prop="${prop:0:3}"
		sed -n "${iter}p" < samples.txt > ${breed}-$prop.txt
		sed -n "1,${cutoff}p" < shuffledfinSNPs.txt > ${breed}-$prop-chosenSNPs.txt
		plink -bfile $ref -chr-set 38 -allow-extra-chr -keep ${breed}-$prop.txt -extract ${breed}-$prop-chosenSNPs.txt -out test-${breed}-${prop}-${id} -make-bed
		echo -e "$id\t$id\t0\t0\t0\t-9" > temp.fam #do this to merge later
		mv temp.fam test-${breed}-$prop-$id.fam
		name=${breed}${prop}
		oldcutoff=$cutoff
		iter=$(($iter+1))
		echo test-${breed}-$prop-$id > listoffiles.txt
	elif [ $iter == $N ]; then #for last entry
		echo "$oldcutoff"
		echo "$cutoff"
		len=$(sed -n "${oldcutoff},${nlines}p" < shuffledfinSNPs.txt | wc -l)
                breed=$(sed -n "${iter}p" < breeds.txt)
                prop=$((1000 * len / nlines))
                echo "$prop"
		prop="${prop:0:3}"
                sed -n "${iter}p" < samples.txt > ${breed}-$prop.txt
                sed -n "${oldcutoff},${nlines}p" < shuffledfinSNPs.txt > ${breed}-$prop-chosenSNPs.txt
                plink -bfile $ref -chr-set 38 -allow-extra-chr -keep ${breed}-$prop.txt -extract ${breed}-$prop-chosenSNPs.txt -out test-${breed}-$prop-$id -make-bed
		echo -e "$id\t$id\t0\t0\t0\t-9" > temp.fam #do this to merge later
                mv temp.fam test-${breed}-$prop-$id.fam
		name=${name}"X"${breed}${prop}
		echo test-${breed}-$prop-$id >> listoffiles.txt
	else
		iter=$(($iter+1))
		len=$(sed -n "${oldcutoff},${cutoff}p" < shuffledfinSNPs.txt | wc -l)
                breed=$(sed -n "${iter}p" < breeds.txt)
                prop=$((1000 * len / nlines))
                echo "$prop"
		prop="${prop:0:3}"
                sed -n "${iter}p" < samples.txt > ${breed}-$prop.txt
                sed -n "${oldcutoff},${cutoff}p" < shuffledfinSNPs.txt > ${breed}-$prop-chosenSNPs.txt
                plink -bfile $ref -chr-set 38 -allow-extra-chr -keep ${breed}-$prop.txt -extract ${breed}-$prop-chosenSNPs.txt -out test-${breed}-$prop-$id -make-bed
		echo -e "$id\t$id\t0\t0\t0\t-9" > temp.fam #do this to merge later
                mv temp.fam test-${breed}-$prop-$id.fam
		name=${name}"X"${breed}${prop}
		echo test-${breed}-$prop-$id >> listoffiles.txt
		oldcutoff=$cutoff
	fi
done < cutoffs.txt
echo $name
plink -merge-list listoffiles.txt -chr-set 38 -allow-extra-chr -make-bed -out newcompletedsynth/$name 

