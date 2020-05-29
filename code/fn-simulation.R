simulation <- function(name,seed,N,snps,alpha1,alpha2,beta1,beta2,treat,c_mean,c_sd) {
  
  # Example input ------------------------------------------------------------

  # name <- "example"
  # seed <- 1
  # N <- 500000 # Sample size
  # snps <- 200 # Number of SNPs
  # alpha1 <- 10 # Coefficient for Z in exposure
  # alpha2 <- 0.3 # Coefficient for C in exposure
  # beta1 <- log(2) # Coefficient for X in outcome
  # beta2 <- log(1) # Coefficient for C in outcome
  # treat <- 2 # Value added to correct for medication use
  # c_mean <- 10 # Mean value of the covariate
  # c_sd <- 2 # Standard deviation of the covariate

  # Set seed (necessary to recreate analysis) -----------------------------------
  
  set.seed(seed)
  
  # Calculate total sample size -------------------------------------------------
  
  N_total <- 2*N
  
  # Generate genetic variants ---------------------------------------------------
  
  df <- purrr::map_dfc(1:snps, ~ rbinom(N_total,2,0.3))
  colnames(df) <- paste0("z_", 1:ncol(df))
  df$z <- rowSums(df)
  df$z_norm <- (df$z - min(df$z))/(max(df$z)-min(df$z))
  
  # Generate continuous covariate -----------------------------------------------
  
  df$c <- rnorm(N_total,c_mean,c_sd)

  # Generate continuous exposure ------------------------------------------------
  
  df$u <- rnorm(N_total,0,1)     
  df$x <- alpha1*df$z_norm + alpha2*df$c + df$u
  
  # Generatre binary outcome ----------------------------------------------------
  
  df$v <- rnorm(N_total,0,1)
  df$y_cont <- beta1*df$x + beta2*df$c + df$v
  df$p <- exp(df$y_cont) / (1+(exp(df$y_cont)))
  df$y <- rbinom(N_total,1,df$p)
  
  # Generate exposure adjusted for treatment ------------------------------------
  
  df$x_norm <- (df$x - min(df$x))/(max(df$x)-min(df$x))
  df$t <- rbinom(N_total,1,df$x_norm)
  df$xt <- ifelse(df$t==1, df$x+treat, df$x)
  
  # Generate exposure with treated people removed -----------------------------
  
  df$xr <- ifelse(df$t==1, NA, df$x)
  
  # Summary table -------------------------------------------------------------
  
  s <- df[,-c(1:200)]
  
  # Split sample --------------------------------------------------------------
  
  split_sample1 <- sample(1:N_total,N,replace = FALSE)
  sample1 <- df[split_sample1,]
  
  split_sample2 <- setdiff(1:N_total,split_sample1)
  sample2 <- df[split_sample2,]
    
  # Calculate genetic associations ---------------------------------------------
  
  assoc <- data.frame(trait=character(),
                      snp=character(), 
                      est=numeric(),
                      se=numeric(),
                      stringsAsFactors=FALSE) 
  
  df$y <- factor(df$y)
  
  for (v in 1:snps) {
    
    for (u in c("x","xt","xr")) {
      
      f <- paste0(u," ~ z_",v)
      m <- lm(as.formula(f),data = sample1)
      c <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type="HC3")) # type = "HC1" to reproduce Stata
      assoc[nrow(assoc)+1,] <- c(u,paste0("z_",v),c[2,1],c[2,2])
      
      f <- paste0(u," ~ z_",v," + c")
      m <- lm(as.formula(f),data = sample1)
      c <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type="HC3")) # type = "HC1" to reproduce Stata
      assoc[nrow(assoc)+1,] <- c(paste0(u,"c"),paste0("z_",v),c[2,1],c[2,2])
      
    }
    
    f <- paste0("y ~ z_",v)
    m <- glm(as.formula(f),data = sample2, family = "binomial")
    c <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type="HC3")) # type = "HC1" to reproduce Stata
    assoc[nrow(assoc)+1,] <- c("y",paste0("z_",v),c[2,1],c[2,2])
    
  }
  
  assoc$est <- as.numeric(assoc$est)
  assoc$se <- as.numeric(assoc$se)
  
  # Conduct MR ------------------------------------------------------------------
  
  mr <- data.frame(exposure=character(),
                   outcome=character(), 
                   b=numeric(),
                   se=numeric(),
                   lci=numeric(),
                   uci=numeric(),
                   pval=numeric(),
                   analysis=character(),
                   stringsAsFactors=FALSE)
  
  for (p in c("x","xc","xt","xtc","xr","xrc")) {
    
    # MR: Exposure-Outcome ------------------------------------------------------

    mr_input <- MendelianRandomization::mr_input(bx = assoc[assoc$trait==p,]$est,
                                                 bxse = assoc[assoc$trait==p,]$se,
                                                 by = assoc[assoc$trait=="y",]$est,
                                                 byse = assoc[assoc$trait=="y",]$se,
                                                 exposure = p, 
                                                 outcome = "y")
    mr_output <- MendelianRandomization::mr_ivw(mr_input)
    mr[nrow(mr)+1,] <- c(mr_output@Exposure,mr_output@Outcome,
                         mr_output@Estimate,mr_output@StdError,
                         mr_output@CILower,mr_output@CIUpper,
                         mr_output@Pvalue,"MR")
    
  }
  
  mr$b <- as.numeric(mr$b)
  mr$lci <- as.numeric(mr$lci)
  mr$uci <- as.numeric(mr$uci)
  mr$se <- as.numeric(mr$se)
  mr$pval <- as.numeric(mr$pval)
  
  # Record true associations ----------------------------------------------------
  
  for (dat in c("s","sample1","sample2")) {
    
    m <- glm(y ~ x + c,data = get(dat), family = "binomial")
    c <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type="HC3")) # type = "HC1" to reproduce Stata
    mr[nrow(mr)+1,] <- c("x","y",c[2,1],c[2,2],NA,NA,NA,paste0("Truth (",dat,")"))
    
    m <- lm(y_cont ~ x + c,data = get(dat))
    c <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type="HC3")) # type = "HC1" to reproduce Stata
    mr[nrow(mr)+1,] <- c("x","y_cont",c[2,1],c[2,2],NA,NA,NA,paste0("Truth (",dat,")"))
    
  }
  
  # Save ------------------------------------------------------------------------
  
  write.csv(mr,paste0("output/simulation_results-",name,".csv"),row.names = FALSE)
  
  return(mr)
  
}