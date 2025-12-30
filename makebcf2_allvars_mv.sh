#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=100G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load plink
module load bwa
module load samtools
module load matlab
module load samtools
module load htslib
module load parallel

inp=$1
pref="${inp##*/}"
echo "$pref"
scope_dir=/u/home/g/gkislik/project-pellegrini/cat-pipeline/SCOPE/build
seed=12345


dir1=/u/home/g/gkislik/project-pellegrini/all_geno/final/
ref=/u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa

$dir1/samtools-1.18/samtools index -c $dir1/REF-REDO/${pref}_test.bam



task () {
	$dir1/samtools-1.18/samtools view -b $dir1/REF-REDO/${pref}_test.bam chr$1 > $dir1/REF-REDO/${pref}_chr$1.bam
	$dir1/bcftools-1.18/bcftools mpileup -f $ref $dir1/REF-REDO/${pref}_chr$1.bam | $dir1/bcftools-1.18/bcftools call -mv -V indels -Ob -o $dir1/ALL_VAR_REF/${pref}_chr$1.bcf.gz
	$dir1/bcftools-1.18/bcftools index $dir1/ALL_VAR_REF/${pref}_chr${1}.bcf.gz
	$dir1/bcftools-1.18/bcftools reheader --samples ${pref}_${i}_rename.txt $dir1/ALL_VAR_REF/${pref}_chr$1.bcf.gz -o $dir1/ALL_VAR_REF/${pref}_chr${1}_corrected.bcf.gz	
	$dir1/bcftools-1.18/bcftools index $dir1/ALL_VAR_REF/${pref}_chr${1}_corrected.bcf.gz
}

rm ${pref}_lc.txt

for i in $(seq 1 38);
	do echo $i
	echo "${pref}" > ${pref}_${i}_rename.txt 
	task ${i} &
	echo "$dir1/ALL_VAR_REF/${pref}_chr${i}_corrected.bcf.gz" >> ${pref}_lc.txt	 
##$dir1//samtools-1.18/samtools view -b $dir1/REF-REDO/${pref}_test.bam chr$i > $dir1/REF-REDO/${pref}_chr$i.bam
##$dir1/bcftools-1.18/bcftools mpileup -f $ref $dir1/REF-REDO/${pref}_test.bam | $dir1/bcftools-1.18/bcftools call -mv -Ob -o $dir1/ALL_VAR_REF/${pref}_test.bcf.gz
done

wait

 

$dir1/bcftools-1.18/bcftools concat -f ${pref}_lc.txt -Ob -o $dir1/ALL_VAR_REF/${pref}_test.bcf.gz
$dir1/bcftools-1.18/bcftools filter -i 'QUAL>=20' $dir1/ALL_VAR_REF/${pref}_test.bcf.gz -Ob -o $dir1/ALL_VAR_REF/${pref}_test_qual20.bcf.gz

$dir1/bcftools-1.18/bcftools norm -f $ref -m -both $dir1/ALL_VAR_REF/${pref}_test_qual20.bcf.gz -Ob \
        -o $dir1/ALL_VAR_REF/${pref}.norm.bcf.gz 



# Give stable, deterministic IDs so matching is unambiguous
bcftools annotate $dir1/ALL_VAR_REF/${pref}.norm.bcf.gz -Ob \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -o $dir1/ALL_VAR_REF/${pref}.norm.id.bcf.gz
$dir1/bcftools-1.18/bcftools index $dir1/ALL_VAR_REF/${pref}.norm.id.bcf.gz

echo $line

rm $dir1/ALL_VAR_REF/${pref}_chr*.bcf.gz
rm ${pref}_*_rename.txt
rm ${pref}_lc.txt

echo "done"



