# clear workspace
# rm(list = ls())

# set working directory
setwd("/Users/jaredfrees/Documents/bioinformatics/learning/genomics_bootcamp/")

# run PLINK QC
system("./plink --bfile DataDryad/ADAPTmap_genotypeTOP_20160222_full --cow --nonfounders --allow-no-sex --recode vcf --out goat_plink_output/ADAPTmap_TOP")


"
Missingness per SNP: --geno value
Missingness per individual: --mind value
Minor allele frequency: --maf value
Hardy-Weinberg threshold: --hwe value
"
# Missingness per SNP: 0.1; Missingness per individual: 0.1; Minor allele frequency: 0.05; hwe 0.0000001
system("./plink --bfile DataDryad/ADAPTmap_genotypeTOP_20160222_full --cow --nonfounders --autosome --allow-no-sex \\
       --recode vcf \\
       --geno .1 \\
       --mind .1 \\
       --maf .05 \\
       --hwe .0000001 \\
       --make-bed \\
       --out goat_plink_output/ADAPTmap_TOP")


## PCA ##

## Genetic distances between individuals
system("./plink --cow --allow-no-sex --nonfounders --file goat_plink_output/ADAPTmap_TOP --distance-matrix --out goat_plink_output/data_for_pca")

## Load data
dist_populations <- read.table("goat_plink_output/data_for_pca.mdist", header=F)

### Extract breed names
fam <- data.frame(famids=read.table("goat_plink_output/data_for_pca.mdist.id")[,1])

### Extract individual names 
famInd <- data.frame(IID=read.table("goat_plink_output/data_for_pca.mdist.id")[,2])

## Perform PCA using the cmdscale function 
# Time intensive step - takes a few minutes with the 4.5K animals
mds_populations <- cmdscale(dist_populations, eig=T, 5)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam, famInd, mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig)) * 100, 2)


# Visualize PCA in tidyverse
# Load tidyverse
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

# PCA plot
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`, color = famids), show.legend = FALSE ) + 
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of wordwide goat populations",
       x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) + 
  theme_minimal()

