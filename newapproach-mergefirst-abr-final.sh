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
main_dir=/u/home/g/gkislik/project-pellegrini/all_geno/final
ref=/u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa
variants=$main_dir/ALL_VAR_REF/EXP_REF4_0/sort_NOFST.txt #path to variants to filter for
regions=$main_dir/ALL_VAR_REF/EXP_REF4_0/NOFST_REGIONS.txt # where to call sample file
cluster_file=$main_dir/ALL_VAR_REF/EXP_REF3_FILT_SEP/POST_CORRECTION/BSPEC/EXP-G-corrected-clust.txt #path to 3 column cluster file
seed=12345

refname="${ref##*/}"
sname="${inp##*/}"


##bwa mem -t 64 /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa $main_dir/samples/${inp}_1.fastq* $main_dir/samples/${inp}_2.fastq* | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools sort -l 9  - -o $main_dir/samples/${inp}_test.bam

##/u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools index -c $main_dir/samples/${inp}_test.bam

##$main_dir/bcftools-1.18/bcftools mpileup -f $ref $main_dir/samples/${inp}_test.bam | $main_dir/bcftools-1.18/bcftools call -mA -V indels -Ob -o $main_dir/samples/${inp}_test.bcf.gz
##$main_dir/bcftools-1.18/bcftools filter -i 'QUAL>=20' $main_dir/samples/${inp}_test.bcf.gz -Ob -o $main_dir/samples/${inp}_test_qual20.bcf.gz


task () {
	$main_dir/samtools-1.18/samtools view -b $main_dir/samples/${inp}_test.bam chr$1 > $main_dir/samples/${inp}_chr$1.bam
        $main_dir/samtools-1.18/samtools index $main_dir/samples/${inp}_chr$1.bam
        $main_dir/bcftools-1.18/bcftools mpileup -f $ref -R $regions $main_dir/samples/${inp}_chr$1.bam | $main_dir/bcftools-1.18/bcftools call -m -V indels -Ob -o $main_dir/samples/${inp}_chr$1.bcf.gz
        $main_dir/bcftools-1.18/bcftools index -f $main_dir/samples/${inp}_chr$1.bcf.gz
        $main_dir/bcftools-1.18/bcftools reheader --samples ${sname}_${i}_rename.txt $main_dir/samples/${inp}_chr$1.bcf.gz -o $main_dir/samples/${inp}_chr$1_corrected.bcf.gz
        $main_dir/bcftools-1.18/bcftools index $main_dir/samples/${inp}_chr$1_corrected.bcf.gz

}

##rm ${sname}_lc.txt

##for i in $(seq 1 38);
##       do echo $i
##       echo "${sname}" > ${sname}_${i}_rename.txt
##       task ${i} &
##       echo "$main_dir/samples/${sname}_chr${i}_corrected.bcf.gz" >> ${sname}_lc.txt
##done

##wait


$main_dir/bcftools-1.18/bcftools concat -f ${sname}_lc.txt -Ob -o $main_dir/samples/${inp}_test.bcf.gz

rm ${sname}_lc.txt
rm ${sname}_*_rename.txt

$main_dir/bcftools-1.18/bcftools norm -f $ref -m - $main_dir/samples/${inp}_test.bcf.gz -Ob \
        -o $main_dir/samples/${inp}.norm.bcf.gz

# Give stable IDs
bcftools annotate $main_dir/samples/${inp}.norm.bcf.gz -Ob \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -o $main_dir/samples/${inp}.norm.id.bcf.gz
$main_dir/bcftools-1.18/bcftools index $main_dir/samples/${inp}.norm.id.bcf.gz




echo " beginning merge"
plink -bcf $main_dir/samples/$inp.norm.id.bcf.gz -chr-set 38 -allow-extra-chr -make-bed -out $main_dir/samples/${sname} 

awk '{print $2}' $main_dir/samples/${sname}.bim | sort > ${sname}_snps.txt
sort $variants > sorted_variants.txt
comm -12 ${sname}_snps.txt sorted_variants.txt > ${sname}_common_snps.txt


plink -bfile /u/home/g/gkislik/project-pellegrini/all_geno/final/ALL_VAR_REF/EXP_REF4_0/QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8 -bmerge $main_dir/samples/${sname} -make-bed -chr-set 38 -allow-extra-chr -out $main_dir/samples/merged/${sname}_merge 
awk '{print $1}' $main_dir/samples/merged/${sname}_merge.fam > $main_dir/samples/merged/tempcol.txt
cat $main_dir/samples/merged/tempcol.txt | xargs -I {} basename {} > $main_dir/samples/merged/tempcol2.txt
awk -F "_" '{print $1" "$1" 0 0 0 -9"}' $main_dir/samples/merged/tempcol2.txt > $main_dir/samples/merged/tempcol3.txt
mv $main_dir/samples/merged/${sname}_merge.fam $main_dir/samples/merged/old_${sname}_merge.fam
mv $main_dir/samples/merged/tempcol3.txt $main_dir/samples/merged/${sname}_merge.fam


plink -bfile $main_dir/samples/merged/${sname}_merge -chr-set 38 -allow-extra-chr -extract ${sname}_common_snps.txt -make-bed -out $main_dir/samples/merged/filt-${sname}
plink -bfile $main_dir/samples/merged/filt-${sname} -fst -within $cluster_file -chr-set 38 -allow-extra-chr -out $main_dir/samples/merged/filt-${sname}

awk '$5>0.35 && $5 != "nan" {print $2}' $main_dir/samples/merged/filt-${sname}.fst > $main_dir/samples/merged/filt-${sname}-var.txt

cat $main_dir/samples/merged/filt-${sname}-var.txt $main_dir/ALL_VAR_REF/EXP_REF4_0/BSPEC/10kSNPs.txt | sort -u > $main_dir/samples/merged/filt-${sname}-var2.txt

plink -bfile $main_dir/samples/merged/filt-${sname} -extract $main_dir/samples/merged/filt-${sname}-var2.txt -chr-set 38 -allow-extra-chr -make-bed -out $main_dir/samples/merged/fst-${sname}

plink -bfile $main_dir/samples/merged/fst-${sname} -chr-set 38 -allow-extra-chr -freq -within $cluster_file -out $main_dir/samples/merged/fst-${sname}

echo "merge done, converting to PLINK"
head -n 66 $main_dir/samples/merged/fst-${sname}.frq.strat > $main_dir/samples/outputs/hf_${sname}.txt
#echo "Running SCOPE"
$scope_dir/scope -g $main_dir/samples/merged/fst-${sname} -k 65 -seed $seed -nt 64 -o $main_dir/samples/outputs/${sname}_ -freq $main_dir/samples/merged/fst-${sname}.frq.strat

#generate output
echo $sname
Rscript $main_dir/REF-REDO/procQhat-synth.R $main_dir/samples/outputs/${sname}_Qhat.txt $main_dir/samples/merged/${sname}_merge.fam $main_dir/samples/outputs/hf_${sname}.txt $cluster_file ${sname}

echo "Done"
