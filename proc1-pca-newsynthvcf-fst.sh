#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=2:00:00,h_data=100G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load plink
module load bwa
module load samtools
module load matlab

inp=$1
scope_dir=/u/home/g/gkislik/project-pellegrini/cat-pipeline/SCOPE/build
seed=12345
#ref=/u/project/pellegrini/gkislik/all_geno/final/LABR-POOD-CKSP/ref-SRR8614076
#ref=fst-0p1-ref-SRR8614076
ref=/u/project/pellegrini/gkislik/all_geno/final/LABR-POOD-CKSP/filt-ref-SRR8614076-0p03fst

#bwa mem -t 64 /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa ${inp}_1.fastq.gz ${inp}_2.fastq.gz | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools view -u -  | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools sort -l 9  - -o ${inp}_test.bam
#/u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools index -c ${inp}_test.bam
#/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools mpileup -f /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa ${inp}_test.bam -R /u/project/pellegrini/gkislik/all_geno/final/LABR-POOD-CKSP/ref-SRR8614076-pos.txt | /u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools call -mA -Ob -o ${inp}_test.bcf.gz
#samtools coverage ${inp}_test.bam > ${inp}_stats.txt

echo " beginning merge"
#plink -bcf ${inp}_test.bcf.gz -make-bed -chr-set 38 -allow-extra-chr -snps-only -out ${inp}_test -double-id
#plink -vcf ${inp}.vcf.gz -make-bed -chr-set 38 -allow-extra-chr -snps-only -out ${inp}_test -double-id
#plink -bfile ${inp}_test -chr-set 38 -allow-extra-chr -double-id -make-bed -snps-only -set-missing-var-ids @_# -out ${inp}_test2

#awk '$2=$1"_"$4' ${inp}_test.bim > ${inp}_test_temp.txt
#mv ${inp}_test_temp.txt ${inp}_test.bim

#no need to intentionally nuke plink to make a missnps file for synth data
#plink -bfile $ref -bmerge ${inp}_test2 -make-bed -chr-set 38 -allow-extra-chr -out ${inp}_test_merged
plink -bfile ${inp} -bmerge $ref -make-bed -chr-set 38 -allow-extra-chr -out ${inp}_test_merged -indiv-sort 0
plink -bfile ${inp}_test_merged -chr-set 38 -allow-extra-chr --fill-missing-a2 -make-bed -out ${inp}_merge  

#plink -bfile ${inp}_test2 -exclude ${inp}_test_merged-merge.missnp -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_nm
#plink -bfile $ref -exclude ${inp}_test_merged-merge.missnp -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_ref
#plink -bfile ${inp}_ref -bmerge ${inp}_nm -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_merge

#### old filtering, replace with fst #### 
#plink -bfile ${inp}_merge -maf 0.03 -geno 0.1 -make-bed -chr-set 38 -allow-extra-chr -out ${inp}_prefin
#plink -bfile ${inp}_prefin -indep-pairwise 50 5 0.5 -chr-set 38 -allow-extra-chr -out ${inp}_prefin_LD
#plink -bfile ${inp}_prefin -extract ${inp}_prefin_LD.prune.in -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_fin 
####

#plink -bfile ${inp}_merge -indep-pairwise 50 5 0.5 -chr-set 38 -allow-extra-chr -out ${inp}_prefin_LD
####

#plink -bfile ${inp}_test_merged -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_prefin -maf 0.03 #-extract ${inp}_prefin_LD.prune.in -maf 0.05 -geno 0.1
#plink -bfile ${inp}_prefin -within sortsixclust.txt -fst -chr-set 38 -allow-extra-chr -out ${inp}_fst
#awk '$5>0.1 && $5 != "nan" {print $2}' ${inp}_fst.fst > ${inp}_var.txt
#plink -bfile ${inp}_prefin -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_fin -extract ${inp}_var.txt

####try just the selected SNPs####
#cat ${inp}_var.txt uSNPs.txt > ${inp}_intSNPs.txt
#sort -u ${inp}_intSNPs.txt > ${inp}_ufinSNPs.txt
#plink -bfile ${inp}_test_merged -chr-set 38 -allow-extra-chr -make-bed -out ${inp}_fin -extract ${inp}_var.txt #_ufinSNPs.txt 
#only extract informative SNPs, no need to FST filter twice
plink -bfile ${inp}_merge -keep-allele-order -chr-set 38 -allow-extra-chr -freq -within updsixclust.txt -out ${inp}_merge-freq
#echo "Running SCOPE"
$scope_dir/scope -g ${inp}_merge -k 50 -seed $seed -nt 64 -o ${inp}_test_SUP -freq ${inp}_merge-freq.frq.strat

#generate output

awk '{print $3}' ${inp}_merge-freq.frq.strat > ${inp}_p1.txt
head -n 51 ${inp}_p1.txt > ${inp}_p2.txt
tail -n +2 ${inp}_p2.txt > ${inp}_clust.txt
awk '{print $1}' ${inp}_test_SUPQhat.txt > ${inp}_res.txt

paste ${inp}_clust.txt ${inp}_res.txt > ${inp}_out.txt

#### old pipeline ####
#echo "Try IBS"
#mkdir -p $inp

# Process SNP blocks and calculate IBS distances
#awk '{print $2}' ${inp}_fin.bim | split -l 100 - $inp/snp_block_


#i=1
#for file in $inp/snp_block_*; do
#  echo "* $(cat $file)" > $inp/snp_block.txt
#  start_pos=$(head -n 1 $file)
#  end_pos=$(tail -n 1 $file)

#  echo $i
  #wc -l $merged_dir/$sample/snp_block.txt
  #echo $start_pos
  #echo $end_pos

#  plink --bfile ${inp}_fin \
#    --distance ibs \
#    --out $inp/ibs_distance \
#    --chr-set 38 \
#    --allow-extra-chr \
#    --extract $inp/snp_block.txt >> /u/project/pellegrini/gkislik/COMPUTE_DISTANCE_DEBUG.log
#  col1=$(cut -f1 $inp/ibs_distance.mibs)
#  col1_line=$(echo "$col1" | tr '\n' ' ')

  # echo $col1_line

#  if [[ $i -eq 1 ]]; then
#    echo $col1_line > $inp.DISTANCE.txt
#    echo -e "${start_pos}\t${end_pos}" > $inp.snp_block_coords.txt
#  else
#    echo $col1_line >> $inp.DISTANCE.txt
#    echo -e "${start_pos}\t${end_pos}" >> $inp.snp_block_coords.txt
#  fi
#  rm $file
#  ((i++))
#done
#rm -r $inp

#get JSON
#matlab -sd /u/project/pellegrini/gkislik/scripts/ -batch "newG_breed_inference_ACTIVE $inp /u/home/g/gkislik/project-pellegrini/scripts"


echo "Done"
