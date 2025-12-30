#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=6:00:00,h_data=75G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load plink
module load bwa
module load samtools
module load matlab
module load R
inp=$1
scope_dir=/u/home/g/gkislik/project-pellegrini/cat-pipeline/SCOPE/build
seed=12345

refname="${ref##*/}"
sname="${inp##*/}"

cp /u/project/pellegrini/gkislik/all_geno/final/newsynth/newcompletedsynth/${sname}.bcf.gz /u/project/pellegrini/gkislik/all_geno/final/REF-REDO/
cd /u/project/pellegrini/gkislik/all_geno/final/REF-REDO/
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools index ${sname}.bcf.gz


echo " beginning merge"
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools merge filt-bcftools-ref-redone.imputed.filtered.bcf.gz ${sname}.bcf.gz -Ob -o ${sname}_merge.bcf.gz
plink -bcf ${sname}_merge.bcf.gz -chr-set 38 -allow-extra-chr -make-bed -out ${sname}_merge -double-id -snps-only
awk '{print $1}' ${sname}_merge.fam > ${sname}_tempcol.txt
cat ${sname}_tempcol.txt | xargs -I {} basename {} > ${sname}_tempcol2.txt
awk -F "_" '{print $1" "$1" 0 0 0 -9"}' ${sname}_tempcol2.txt > ${sname}_tempcol3.txt
mv ${sname}_merge.fam old_${sname}_merge.fam
mv ${sname}_tempcol3.txt ${sname}_merge.fam

plink -bfile ${sname}_merge -chr-set 38 -allow-extra-chr -extract imputed_filtered_finSNPs.txt -make-bed -out filt-${sname}


echo "merge done, converting to PLINK"
plink -bfile filt-${sname} -keep-allele-order -chr-set 38 -allow-extra-chr -freq -within redone-clust.txt -out ${sname}-freq
head -n 49 ${sname}-freq.frq.strat > hf_${sname}.txt
#echo "Running SCOPE"
$scope_dir/scope -g filt-${sname} -k 48 -seed $seed -nt 64 -o ${sname}_ -freq ${sname}-freq.frq.strat

#generate output
echo $sname
Rscript procQhat-synth.R ${sname}_Qhat.txt ${sname}_merge.fam hf_${sname}.txt redone-clust.txt ${sname}

echo "Done"
