
### SCRIPTS FOR EXPANDING THE REFERENCE POPULATION <br />
loop.sh was used to submit jobs to the UCLA Hoffman cluster <br />
+ it takes a list of forward read fastq files with suffix _1.fastq.gz 
+ makebcf.sh performs bwa alignment 
+ makebcf2_allvars_mv.sh maps and identifies variants, as well as to generate PLINK files for merging <br />

cleannbcf.sh is used to normalize bcf.gz files and set stable variant ids <br />

mergebcfs_with_bcftools_final.sh is used to merge bcf.gz files <br />
+ it takes a list of bcf.gz files to merge in one step (change within script) <br />

### SCRIPTS FOR GETTING BREED SPECIFIC SNPS <br />
getSNPs.sh is used to find the breed-specific SNPs and needs the reference population and a cluster file as arguments. <br />
+ getSNPs_helper.sh, getS.R, and getHM.R are helper files for getSNPs.sh which help to parallelize each breed's process, generate the distribution plots, and heatmaps using breed-specific SNPs, respectively. <br />

### SCRIPT FOR GENERATING SYNTHETIC SAMPLES <br />
synth_helper.sh is used and takes arguments N (desired number of samples to choose), clust (a clusters file) <br />

### SCRIPTS FOR ANALYSIS
synth-FINALPROC.sh and fst-FINALPROC-final.sh are used for submitting jobs, similar to loop.sh. 
They take a text file with samples (full path to forward read of sample fastq - even if you're starting from bcf like for synthetic files, just change .bcf.gz to _1.fastq.gz) on each line. 
+ newapproach-mergefirst-abr-final.sh is used by fst-FINALPROC-final.sh to analyze the non-synthetic samples
+ newapproach-synth-final.sh is used by synth-FINALPROC.sh to analyze synthetic bcf samples made by make-synthetic-bcf.sh
bootSCOPE.sh and bootSCOPE_helper.sh were used to bootstrap standard deviations for the prediction of breeds within the reference population
+ bootSCOPE_helper.sh does most of the heavy lifting, bootSCOPE.sh is used for job submission. 
