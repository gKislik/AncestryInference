
### SCRIPTS FOR EXPANDING THE REFERENCE POPULATION
loop.sh was used to submit jobs to the UCLA Hoffman cluster
makebcf.sh performs bwa alignment
makebcf2.sh maps and identifies variants, as well as to generate PLINK files for merging
mergeref.sh is used to add the newly made plink files to an existing reference population

### SCRIPTS FOR GETTING BREED SPECIFIC SNPS
getSNPs.sh is used to find the breed-specific SNPs and needs the reference population and a cluster file as arguments. 
getS.R and getHM.R are helper files for getSNPs.sh which help to generate the distribution plots and heatmaps using breed-specific SNPs. 


