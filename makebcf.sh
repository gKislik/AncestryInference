#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=100G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

module load bcftools
module load plink
module load bwa
module load samtools
module load matlab

inp=$1
pref="${inp##*/}"
echo "$pref"
scope_dir=/u/home/g/gkislik/project-pellegrini/cat-pipeline/SCOPE/build
seed=12345


/u/project/pellegrini/gkislik/scripts/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 80 /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa ${inp}_1.fastq.gz ${inp}_2.fastq.gz | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools view -u -  | /u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools sort -l 9  - -o /u/home/g/gkislik/project-pellegrini/all_geno/final/REF-REDO/${pref}_test.bam

qsub makebcf2_allvars_mv.sh ${inp}

#qsub makebcf2.sh ${inp} 



#/u/home/g/gkislik/project-pellegrini/all_geno/final/samtools-1.18/samtools index -c /u/home/g/gkislik/project-pellegrini/all_geno/final/REF-REDO/${pref}_test.bam
#/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools mpileup -f /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa -R /u/project/pellegrini/gkislik/ref-seq/A2P_pos.txt /u/home/g/gkislik/project-pellegrini/all_geno/final/REF-REDO/${pref}_test.bam | /u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools call -mA -Ob -o /u/home/g/gkislik/project-pellegrini/all_geno/final/REF-REDO/${pref}_test.bcf.gz
#plink -bcf /u/home/g/gkislik/project-pellegrini/all_geno/final/REF-REDO/${pref}_test.bcf.gz -set-missing-var-ids @_# -double-id -chr-set 38 -allow-extra-chr -make-bed -snps-only -out /u/home/g/gkislik/project-pellegrini/all_geno/final/LABR-POOD-CKSP/${pref}_test

