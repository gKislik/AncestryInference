#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=4:00:00,h_data=70G

source ~/.bashrc
source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load plink
module load R
module load bcftools

conda activate hap_env

N=$1
breeds=$2

## generate random proportions & breeds


shuf $breeds > shufbreeds.txt
head -n $N shufbreeds.txt > chosen.txt

shuf -i 1-100 -n $N > numbers.txt 

sum=$(awk '{s+=$1} END {print s}' numbers.txt)

echo $sum

props=()
breeds=()
for i in $(seq 1 $N); do
	b=$(head -n $i chosen.txt | tail -1)
	n=$(head -n $i numbers.txt | tail -1) 
	echo $n
	n2=$(echo "scale=4; $n / $sum" | bc)
	echo $n2
	props+=($n2)
	breeds+=($b)
	echo $i
done

###make entries sum to 1 (fix floating point issues)
sum2=0
for i in "${props[@]}"; do
  sum2=$(echo "scale=4; $sum2 + $i" | bc)
done
echo "$sum2"

diff=$(echo "1 - $sum2" | bc)
echo "$diff"
le=${props[-1]}
echo "$le"
up_le=$(echo "scale=4; $le + $diff" | bc)


last_index=$(( ${#props[@]} - 1 ))
props[last_index]="$up_le"

#echo ${props[@]}
echo -e "${breeds[@]}\n${props[@]}"
fname="${breeds[@]}_${props[@]}"
echo "${fname}"
new_fname="${fname// /_}"
echo "$new_fname"

echo -e "${breeds[@]}\n${props[@]}" > MISC/${new_fname}.txt

paste  p1_modfile_forhaptools.txt MISC/${new_fname}.txt >  MISC/${new_fname}_forhaptools.txt


base="/u/home/g/gkislik/project-pellegrini/all_geno/final/ALL_VAR_REF"
haptools simgenotype \
--model MISC/${new_fname}_forhaptools.txt \
--mapdir $base/EXP_REF4_0/SYNTH/MAP \
--ref_vcf $base/EXP_REF4_0/SYNTH/QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8.vcf \
--sample_info two-col.txt \
--seed 123 \
--chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38 \
--out SYNTH_FILES/${new_fname}_simgenotype.vcf.gz

/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools view SYNTH_FILES/${new_fname}_simgenotype.vcf.gz -Ob -o SYNTH_FILES/${new_fname}_simgenotype.bcf.gz
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools index SYNTH_FILES/${new_fname}_simgenotype.bcf.gz
