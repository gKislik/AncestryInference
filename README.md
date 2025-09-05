
### SCRIPTS FOR EXPANDING THE REFERENCE POPULATION <br />
loop.sh was used to submit jobs to the UCLA Hoffman cluster <br />
+ it takes a list of forward read fastq files with suffix _1.fastq.gz 
+ makebcf.sh performs bwa alignment 
+ makebcf2.sh maps and identifies variants, as well as to generate PLINK files for merging <br />

mergebcfs_with_bcftools_final.sh is used to merge bcf.gz files <br />
+ it takes a list of bcf.gz files to merge in one step (change within script) <br />
+ missing genotypes stay missing (-0 tag not used) <br />

### SCRIPTS FOR GETTING BREED SPECIFIC SNPS <br />
getSNPs.sh is used to find the breed-specific SNPs and needs the reference population and a cluster file as arguments. <br />
+ getSNPs_helper.sh, getS.R, and getHM.R are helper files for getSNPs.sh which help to parallelize each breed's process, generate the distribution plots, and heatmaps using breed-specific SNPs, respectively. <br />

### SCRIPT FOR GENERATING SYNTHETIC SAMPLES <br />
make-synthetic-bcf.sh is used and takes arguments N (desired number of samples to choose), clust (a clusters file), and id (use aaaa, this just serves as a placeholder for analysis after merging) <br />

### SCRIPTS FOR ANALYSIS
+ synthfst-new-FINALPROC.sh and fst-FINALPROC.sh are used for submitting jobs, similar to loop.sh. They take a text file with samples on each line. 
+ proc1-pca-fst.sh is used by fst-FINALPROC.sh to analyze the non-synthetic samples
+ proc1-pca-newsynthvcf-fst.sh is used by synthfst-new-FINALPROC.sh to analyze synthetic PLINK samples made by make-synthetic-plink.sh
