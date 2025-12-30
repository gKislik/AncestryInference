#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=4:00:00,h_data=30G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load samtools
module load htslib
module load plink

/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools norm -f /u/project/pellegrini/gkislik/ref-seq/UU_Cfam_GSD_1.0_ROSY.fa -m -both $1.bcf.gz -Ob \
	-o $1.norm.bcf.gz


/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools +missing2ref $1.norm.bcf.gz -Ob -o $1.temp.bcf.gz


# Give stable, deterministic IDs so matching is unambiguous
bcftools annotate $1.temp.bcf.gz -Ob \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -o $1.norm.id.bcf.gz
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools index $1.norm.id.bcf.gz

echo $line
rm $1.tem*
