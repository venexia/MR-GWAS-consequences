remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Conduct MR analyses

df <- NULL

files <- gsub("exposure-","",gsub(".csv","",list.files(path = "data/", pattern = "exposure-")))

for (f in files) {
  
  ## Load exposure data
  exp <- data.table::fread(paste0("data/exposure-",f,".csv"), data.table = FALSE)
  
  ## Extract outcome SNPs
  out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP,
                                            outcomes = c("ieu-a-7"))
  
  ## Harmonise exposure and outcome data
  dat <- TwoSampleMR::harmonise_data(exp,out)
  
  ## Perform Mendelian randomization
  res <- TwoSampleMR::mr(dat, method_list = c("mr_wald_ratio","mr_ivw"))
  
  ## Add to main results data frame
  df <- rbind(df,res)
  
}

# Add confidence intervals

df$lci <- df$b - qnorm(0.975)*df$se
df$uci <- df$b + qnorm(0.975)*df$se

# Store unit results

df$b_unit <- df$b
df$lci_unit <- df$lci
df$uci_unit <- df$uci

# Scale results

df$sd <- 19.268
df$b <- df$b_unit*df$sd
df$lci <- df$lci_unit*df$sd
df$uci <- df$uci_unit*df$sd

# Save empirical results

data.table::fwrite(df,"output/empirical.csv")
