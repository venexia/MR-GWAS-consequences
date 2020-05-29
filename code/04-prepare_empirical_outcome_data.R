remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Load BMI SNP list

bmi <- TwoSampleMR::extract_instruments(outcomes = c("ieu-a-2"))

# List GWAS files

files <- gsub("imputed-","",gsub(".txt.gz","",list.files(path = data_path, pattern = "imputed-sbp")))

for (f in files) {
  
  # Read in data
  
  tmp <- data.table::fread(input = paste0(data_path,"imputed-",f,".txt.gz"))
  
  # Add phenotype name
  
  tmp$Phenotype <- f
  
  # Format data for TwoSampleMR package
  
  tmp <- TwoSampleMR::format_data(tmp,
                                  type = "outcome",
                                  snps = bmi$SNP,
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
  
  # Save data
  
  data.table::fwrite(tmp, file = paste0("data/","outcome-",f,".csv"))
  
}