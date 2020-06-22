remove(list = ls())
graphics.off()

# Load libraries

library(magrittr)

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# List GWAS files

files <- gsub("imputed-","",gsub(".txt.gz","",list.files(path = data_path, pattern = "sbp")))

for (f in files) {
  
  # Read in data
    
  tmp <- data.table::fread(paste0(data_path,"imputed-",f,".txt.gz"))
  
  # Restrict to genome wide significant hits
  
  tmp <- tmp[tmp$P_LINREG < 5e-8,]

  # Add phenotype name
  
  tmp$Phenotype <- f
  
  # Format data for TwoSampleMR package
  
  tmp <- TwoSampleMR::format_data(tmp,
                                 type = "exposure",
                                 snp_col = "SNP", 
                                 phenotype_col = "Phenotype",
                                 effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0",
                                 eaf_col = "A1FREQ",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 pval_col = "P_LINREG",
                                 chr_col = "CHR",
                                 pos_col = "BP")
  
  # Clump data
  
  tmp <- TwoSampleMR::clump_data(tmp,
                                 clump_kb = 10000,
                                 clump_r2 = 0.001)
  
  # Make all effects SBP increasing
  
  tmp$beta.original <- tmp$beta.exposure
  tmp$eaf.original <- tmp$eaf.exposure
  tmp$effect_allele.original <- tmp$effect_allele.exposure
  tmp$other_allele.original <- tmp$other_allele.exposure
  tmp[,c("beta.exposure","eaf.exposure","effect_allele.exposure","other_allele.exposure")] <- NULL
  
  tmp <- tmp %>%
    dplyr::mutate(beta.exposure = ifelse(sign(beta.original)==-1, -1*beta.original, beta.original)) %>%
    dplyr::mutate(effect_allele.exposure = ifelse(sign(beta.original)==-1, other_allele.original, effect_allele.original)) %>%
    dplyr::mutate(other_allele.exposure = ifelse(sign(beta.original)==-1, effect_allele.original, other_allele.original)) %>%
    dplyr::mutate(eaf.exposure = ifelse(beta.original <= 0, 1-eaf.original, eaf.original))
  
  tmp[,c("beta.original","eaf.original","effect_allele.original","other_allele.original")] <- NULL
  
  # Save data
  
  data.table::fwrite(tmp, file = paste0("data/","exposure-",f,".csv"))
  
}