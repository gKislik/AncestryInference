
### SCRIPTS FOR EXPANDING THE REFERENCE POPULATION <br />
loop.sh was used to submit jobs to the UCLA Hoffman cluster <br />
makebcf.sh performs bwa alignment <br />
makebcf2.sh maps and identifies variants, as well as to generate PLINK files for merging <br />
mergeref.sh is used to add the newly made plink files to an existing reference population <br />
+ it takes a list of PLINK prefixes to merge with the reference population identified within the script <br />

### SCRIPTS FOR GETTING BREED SPECIFIC SNPS <br />
getSNPs.sh is used to find the breed-specific SNPs and needs the reference population and a cluster file as arguments. <br />
+ getS.R and getHM.R are helper files for getSNPs.sh which help to generate the distribution plots and heatmaps using breed-specific SNPs. <br />


