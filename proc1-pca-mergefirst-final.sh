#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=75G

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

bwa mem -t 64 /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa ${inp}_1.fastq* ${inp}_2.fastq* | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools view -u -  | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools sort -l 9  - -o ${inp}_test.bam
/u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools index -c ${inp}_test.bam
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools mpileup -f /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa ${inp}_test.bam -R /u/project/pellegrini/gkislik/all_geno/final/LABR-POOD-CKSP/ref-SRR8614076-pos.txt | /u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools call -mA -Ob -o ${inp}_test.bcf.gz

/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools norm -f /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa -m -both ${inp}_test.bcf.gz -Ob \
        -o ${inp}.norm.bcf.gz



bcftools annotate ${inp}.norm.bcf.gz -Ob \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -o ${inp}.norm.id.bcf.gz
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools index ${inp}.norm.id.bcf.gz

refname="${ref##*/}"
sname="${inp##*/}"
echo $sname
cp ${inp}.norm.id.bcf.gz /u/project/pellegrini/gkislik/all_geno/final/REF-REDO/
cd /u/project/pellegrini/gkislik/all_geno/final/REF-REDO/


#echo " beginning merge"
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools merge bcftools-ref-redone.imputed.filtered.bcf.gz ${inp}.norm.id.bcf.gz -Ob -o ${sname}_merge.bcf.gz

plink -bcf ${sname}_merge.bcf.gz -chr-set 38 -allow-extra-chr -make-bed -out ${sname}_merge -double-id -snps-only
awk '{print $1}' ${sname}_merge.fam > tempcol.txt
cat tempcol.txt | xargs -I {} basename {} > tempcol2.txt
awk -F "_" '{print $1" "$1" 0 0 0 -9"}' tempcol2.txt > tempcol3.txt
mv ${sname}_merge.fam old_${sname}_merge.fam
mv tempcol3.txt ${sname}_merge.fam

#plink -bfile ${sname}_merge -chr-set 38 -allow-extra-chr -maf 0.01 -make-bed -out ${sname}_merge2 
#plink -bfile ${sname}_merge2 -chr-set 38 -allow-extra-chr -indep-pairwise 50 100 0.3 -out ${sname}

plink -bfile ${sname}_merge -chr-set 38 -allow-extra-chr -fst -within redone-clust.txt -out ${sname}_fst
awk '$5>0.375 && $5 != "nan" {print $2}' ${sname}_fst.fst > ${sname}_var.txt
echo $(wc -l ${sname}_var.txt)
wc -l ${sname}_var.txt
cat ${sname}_var.txt /u/project/pellegrini/gkislik/all_geno/final/REF-REDO/bcftools_imputed_filt_bspecvar.txt | sort -u > ${sname}_var2.txt #${sname}.prune.in
plink -bfile ${sname}_merge -chr-set 38 -allow-extra-chr -extract ${sname}_var2.txt -make-bed -out filt-${sname}


echo "merge done, converting to PLINK"
plink -bfile filt-${sname} -keep-allele-order -chr-set 38 -allow-extra-chr -freq -within redone-clust.txt -out ${sname}-freq
head -n 49 ${sname}-freq.frq.strat > hf_${sname}.txt
#echo "Running SCOPE"
$scope_dir/scope -g filt-${sname} -k 48 -seed $seed -nt 64 -o ${sname}_ -freq ${sname}-freq.frq.strat


#generate output
Rscript procQhat.R ${sname}_Qhat.txt ${sname}_merge.fam hf_${sname}.txt redone-clust.txt ${sname}

echo "Done"
