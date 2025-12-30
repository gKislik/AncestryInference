#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=95G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools

/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools merge -0 -l allids_4.txt -Ob -o QUAL20-0-G.bcf.gz
/u/home/g/gkislik/project-pellegrini/all_geno/final/bcftools-1.18/bcftools index QUAL20-0-G.bcf.gz

