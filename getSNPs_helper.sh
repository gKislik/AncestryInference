!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=70G

source ~/.bashrc
source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load plink
module load R

# the prefix for the reference genome you want to find SNPs with
inp=$1

# the clusters file for plink
clust=$2

# the breed
line=$3

awk -v line=$line '{$3 = ($3 == line ? line : "N")} 1' $clust > $line-clust.txt
plink -bfile $inp -chr-set 38 -allow-extra-chr -fst -within $line-clust.txt -out $line-fst

#try to select top ~2.5k SNPs
awk '$5 != "nan" && $5>0.3 {print}' $line-fst.fst > $line-pre.txt
sort -nk 5 $line-pre.txt | tail -n 2500 > $line-pre2.txt
Rscript getS.R $line-pre2.txt $line
awk '{print $2}' $line-pre2.txt > $line-SNP.txt
        
#generate heatmap to verify clustering
plink -bfile $inp -extract $line-SNP.txt -chr-set 38 -allow-extra-chr -make-bed -out $line-ex
plink -bfile $line-ex -within $clust -freq -chr-set 38 -allow-extra-chr -out $line-freq
Rscript getHM.R $line-freq.frq.strat $line
