remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# List GWAS files

files <- gsub(".txt","",list.files(path = paste0("raw/"), pattern = "sbp"))

for (f in files) {
  
  # Read in data
    
  tmp <- data.table::fread(input = paste0("raw/",f,".txt"))
  
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
                                 pval_col = "P_BOLT_LMM_INF",
                                 chr_col = "CHR",
                                 pos_col = "BP")
  
  # Clump data
  
  tmp$rsid <- tmp$SNP
  tmp <- ieugwasr::ld_clump(dat = tmp,
                            clump_kb = 10000,
                            clump_r2 = 0.001,
                            clump_p = 1)
  
  # Save data
  
  data.table::fwrite(tmp, file = paste0("data/","exposure-",f,".csv"))
  
}