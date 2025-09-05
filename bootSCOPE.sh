#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=1:00:00,h_data=10G

source /u/home/g/gkislik/miniconda3/etc/profile.d/conda.sh
. /u/local/Modules/default/init/modules.sh
module load bcftools
module load plink


inp=$1
clust=$2

N=$(cat $clust | wc -l)

echo $N

n_less=$((N-6))
echo $n_less

scope_dir=/u/home/g/gkislik/project-pellegrini/cat-pipeline/SCOPE/build
boot_dir=$PWD/boot_out


for i in {1..100} ; do
	qsub bootSCOPE_helper.sh $inp $clust $i
	#shuf $clust > $boot_dir/shuffled_clust_rep$i.txt
	#head -n $n_less $boot_dir/shuffled_clust_rep$i.txt | sort -k 3 > $boot_dir/trimclust_rep$i.txt
	#nbreeds=$(awk '{print $3}' $boot_dir/trimclust_rep$i.txt | sort -u | wc -l)
	#np1=$((nbreeds+1))
	#plink -bfile $inp -chr-set 38 -allow-extra-chr -keep $boot_dir/trimclust_rep$i.txt -make-bed -out tempinp_rep$i
	#plink -bfile tempinp_rep$i -chr-set 38 -allow-extra-chr -freq -within $boot_dir/trimclust_rep$i.txt -out $boot_dir/tempinpfreq_rep$i
	#head -n $np1 $boot_dir/tempinpfreq_rep$i.frq.strat > $boot_dir/freq_rep$i.txt 
	#$scope_dir/scope -g tempinp_rep$i -k $nbreeds -seed 12345 -freq $boot_dir/tempinpfreq_rep$i.frq.strat -nt 64 -o temp_out
	#cp temp_outQhat.txt $boot_dir/Qhat_rep$i.txt
done 
	 		
