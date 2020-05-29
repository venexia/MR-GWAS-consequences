remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Load data

df <- NULL

for (c in c("pos","neg")) {
  for (t in c("pos","neg")) {
    tmp <- data.table::fread(paste0("output/simulation_results-c_",c,"_t_",t,".csv"), data.table = FALSE)
    tmp$covariate <- c
    tmp$treatment <- t
    df <- rbind(df,tmp)
  }
}

# Save truth

truth <- df[df$covariate=="pos" & df$treatment=="pos" & df$outcome=="y" & df$analysis=="Truth (s)",]$b
df <- df[df$analysis=="MR",]
  
# Create estimate labels

df$exposure <- ifelse(df$exposure=="xt" & df$treatment=="pos","xt+",df$exposure)
df$exposure <- ifelse(df$exposure=="xt" & df$treatment=="neg","xt-",df$exposure)
df$exposure <- ifelse(df$exposure=="xtc" & df$treatment=="pos","xt+c",df$exposure)
df$exposure <- ifelse(df$exposure=="xtc" & df$treatment=="neg","xt-c",df$exposure)
df$treatment <- NULL

df$covariate <- ifelse(df$covariate=="pos","Positive covariate",df$covariate)
df$covariate <- ifelse(df$covariate=="neg","Negative covariate",df$covariate)

df <- unique(df)

df$estimate <- paste0(sprintf("%.2f",exp(df$b))," (95% CI: ",sprintf("%.2f",exp(df$lci))," to ",sprintf("%.2f",exp(df$uci)),")")

# Create exposure labels

labels <- data.frame(rbind(c("x","No adjustment, correction or selection"),
                           c("xt+","Correction for medication use\n(medication increases exposure)"),
                           c("xt-","Correction for medication use\n(medication decreases exposure)"),
                           c("xr","Selection on medication users"),
                           c("xc","Adjustment for the covariate"),
                           c("xt+c","Adjustment for the covariate and\ncorrection for medication use\n(medication increases exposure)"),
                           c("xt-c","Adjustment for the covariate and\ncorrection for medication use\n(medication decreases exposure)"),
                           c("xrc","Adjustment for the covariate and\nselection on medication users")),
                     stringsAsFactors = FALSE)

colnames(labels) <- c("exposure","exposure_long")

df <- merge(df,labels,by = c("exposure"), all.x = TRUE)

df$exposure_long <- ifelse(!is.na(df$exposure_long),paste0(df$exposure_long,"\n",df$estimate),NA)

df$exposure_long <- factor(df$exposure_long)
df$exposure_long <- factor(df$exposure_long,levels(df$exposure_long)[c(1:2,9:10,3:4,11:12,5:6,15:16,7:8,13:14)])

# Plot results for positive covariate

ggplot(df[df$covariate=="Positive covariate" & df$outcome=="y",], aes(y = exp(b), x = exposure_long)) +
  geom_hline(yintercept = round(exp(truth),2), linetype = "dotted", color = "black") +
  geom_point() + 
  geom_errorbar(aes(ymin = exp(lci), ymax = exp(uci)), width = 0) +
  labs(title = "", x = "Phenotype A genetic associations",
       y = paste0("Odds ratio and 95% confidence interval for the effect of\na unit change in phenotype A on phenotype B\nwhere the dotted line represents the true effect")) +
  scale_y_continuous(trans = "log", limits = c(1.5,2.5), breaks = round(c(exp(truth),seq(1.6,2.4,0.2)),2)) +
  guides(col = guide_legend(ncol=1)) +
  theme_minimal() +
  theme(axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "bottom",
        legend.key.width = unit(0.1,"cm"),
        legend.box.background = element_rect(colour = "dark grey")) +
  coord_flip() 

# Plot results for negative covariate

ggsave("output/simulation_poscovar.jpeg",width = 300, height = 175, unit = "mm", dpi = 600)

ggplot(df[df$covariate=="Negative covariate" & df$outcome=="y",], aes(y = exp(b), x = exposure_long)) +
  geom_hline(yintercept = round(exp(truth),2), linetype = "dotted", color = "black") +
  geom_point() + 
  geom_errorbar(aes(ymin = exp(lci), ymax = exp(uci)), width = 0) +
  labs(title = "", x = "Phenotype A genetic associations",
       y = paste0("Odds ratio and 95% confidence interval for the effect of\na unit change in phenotype A on phenotype B\nwhere the dotted line represents the true effect")) +
  scale_y_continuous(trans = "log", limits = c(1.5,2.5), breaks = round(c(exp(truth),seq(1.6,2.4,0.2)),2)) +
  guides(col = guide_legend(ncol=1)) +
  theme_minimal() +
  theme(axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "bottom",
        legend.key.width = unit(0.1,"cm"),
        legend.box.background = element_rect(colour = "dark grey")) +
  coord_flip() 

ggsave("output/simulation_negcovar.jpeg",width = 300, height = 175, unit = "mm", dpi = 600)