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
synth_dir=/u/home/g/gkislik/project-pellegrini/all_geno/final/ALL_VAR_REF/EXP_REF4_0/SYNTH/SYNTH_FILES
main_dir=/u/home/g/gkislik/project-pellegrini/all_geno/final

variants=$main_dir/ALL_VAR_REF/EXP_REF4_0/10kSNPs_FST0.350.txt #path to variants to filter for
cluster_file=$main_dir/ALL_VAR_REF/EXP_REF3_FILT_SEP/POST_CORRECTION/BSPEC/EXP-G-corrected-clust.txt #path to 3 column cluster file
seed=12345

refname="${ref##*/}"
sname="${inp##*/}"


echo " beginning merge"
##$main_dir/bcftools-1.18/bcftools merge $main_dir/ALL_VAR_REF/EXP_REF4_0/QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-FST0.375_2500.bcf.gz $synth_dir/${sname}.bcf.gz -Ob -o $synth_dir/merged/${sname}_merge.bcf.gz
##plink -bcf $synth_dir/merged/${sname}_merge.bcf.gz -chr-set 38 -allow-extra-chr -make-bed -out $synth_dir/merged/${sname}_merge -double-id -snps-only

plink -bcf $synth_dir/${sname}.bcf.gz -chr-set 38 -allow-extra-chr -make-bed -out $synth_dir/${sname} -double-id -snps-only
plink -bfile $main_dir/ALL_VAR_REF/EXP_REF4_0/QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350 -bmerge $synth_dir/${sname} -chr-set 38 -allow-extra-chr -out $synth_dir/merged/${sname}_merge

awk '{print $1}' $synth_dir/merged/${sname}_merge.fam > $synth_dir/merged/tempcol.txt
cat $synth_dir/merged/tempcol.txt | xargs -I {} basename {} > $synth_dir/merged/tempcol2.txt
awk -F "_" '{print $1" "$1" 0 0 0 -9"}' $synth_dir/merged/tempcol2.txt > $synth_dir/merged/tempcol3.txt
mv $synth_dir/merged/${sname}_merge.fam $synth_dir/merged/old_${sname}_merge.fam
mv $synth_dir/merged/tempcol3.txt $synth_dir/merged/${sname}_merge.fam


plink -bfile $synth_dir/merged/${sname}_merge -chr-set 38 -allow-extra-chr -extract $variants -make-bed -out $synth_dir/merged/filt-${sname}

plink -bfile $synth_dir/merged/filt-${sname} -chr-set 38 -allow-extra-chr -freq -within $cluster_file -out $synth_dir/merged/filt-${sname}

echo "merge done, converting to PLINK"
head -n 66 $synth_dir/merged/filt-${sname}.frq.strat > $synth_dir/outputs/hf_${sname}.txt
#echo "Running SCOPE"
$scope_dir/scope -g $synth_dir/merged/filt-${sname} -k 65 -seed $seed -nt 64 -o $synth_dir/outputs/${sname}_ -freq $synth_dir/merged/filt-${sname}.frq.strat

#generate output
echo $sname
Rscript $main_dir/REF-REDO/procQhat-synth.R $synth_dir/outputs/${sname}_Qhat.txt $synth_dir/merged/${sname}_merge.fam $synth_dir/outputs/hf_${sname}.txt $cluster_file ${sname}

echo "Done"
