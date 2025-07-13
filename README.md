
### SCRIPTS FOR EXPANDING THE REFERENCE POPULATION \n
loop.sh was used to submit jobs to the UCLA Hoffman cluster \n
makebcf.sh performs bwa alignment \n
makebcf2.sh maps and identifies variants, as well as to generate PLINK files for merging \n
mergeref.sh is used to add the newly made plink files to an existing reference population \n

### SCRIPTS FOR GETTING BREED SPECIFIC SNPS \n
getSNPs.sh is used to find the breed-specific SNPs and needs the reference population and a cluster file as arguments. \n
getS.R and getHM.R are helper files for getSNPs.sh which help to generate the distribution plots and heatmaps using breed-specific SNPs. \n


