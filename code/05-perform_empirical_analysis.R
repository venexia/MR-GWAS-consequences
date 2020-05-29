remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Load BMI exposure data

bmi <- TwoSampleMR::extract_instruments(outcomes = c("ieu-a-2"))

# Conduct MR analyses

df <- NULL

files <- gsub(".txt","",list.files(path = "raw/", pattern = "sbp"))

for (f in files) {
  
  # MR of systolic blood pressure on coronary heart disease
  
  ## Load exposure data
  exp1 <- data.table::fread(paste0("data/exposure-",f,".csv"), data.table = FALSE)
  
  ## Extract outcome SNPs
  out1 <- TwoSampleMR::extract_outcome_data(snps = exp1$SNP,
                                            outcomes = c("ebi-a-GCST003116"))
  
  ## Harmonise exposure and outcome data
  dat1 <- TwoSampleMR::harmonise_data(exp1,out1)
  
  ## Perform Mendelian randomization
  res1 <- TwoSampleMR::mr(dat1, method_list = c("mr_wald_ratio","mr_ivw"))
  
  ## Rename outcome
  res1$outcome <- "CAD"
  
  ## Add to main results data frame
  df <- rbind(df,res1)
  
  # # MR of body mass index on systolic blood pressure
  # 
  # ## Load outcome SNPs
  # out2 <- data.table::fread(paste0("data/outcome-",f,".csv"))
  # 
  # ## Harmonise exposure and outcome data
  # dat2 <- TwoSampleMR::harmonise_data(bmi,out2)
  # 
  # ## Perform Mendelian randomization
  # res2 <- TwoSampleMR::mr(dat2, method_list = c("mr_wald_ratio","mr_ivw"))
  # 
  # ## Rename exposure
  # res2$exposure <- "BMI"
  # 
  # ## Add to main results data frame
  # df <- rbind(df,res2)
  
}

# Add confidence intervals

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Store unit results

df$b_unit <- df$b
df$lci_unit <- df$lci
df$uci_unit <- df$uci

# Scale results

# df$sd <- ifelse(df$outcome=="CAD",19.268,1)
# df$b <- ifelse(df$outcome=="CAD",df$b*sbp_sd,df$b)
# df$lci <- ifelse(df$outcome=="CAD",df$lci*sbp_sd,df$lci)
# df$uci <- ifelse(df$outcome=="CAD",df$uci*sbp_sd,df$uci)

df$sd <- 19.268
df$b <- df$b_unit*df$sd
df$lci <- df$lci_unit*df$sd
df$uci <- df$uci_unit*df$sd

# Save empirical results

data.table::fwrite(df,"output/empirical_results.csv")
