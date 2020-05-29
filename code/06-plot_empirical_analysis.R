remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Define parameters

x <- "sbp"
xt <- "sbp_plus_10"
xr <- "sbp_exc_HTN"
xc <- "sbp_BMI_covariable"
xtc <- "sbp_plus_10_BMI_covariable"
xrc <- "sbp_exc_HTN_BMI_covariable"

# Load data

df <- data.table::fread("output/empirical_results.csv", data.table = FALSE)

df$estimate <- paste0(sprintf("%.2f",exp(df$b))," (95% CI: ",sprintf("%.2f",exp(df$lci))," to ",sprintf("%.2f",exp(df$uci)),")")

# Label exposures 

labels <- data.frame(rbind(c(x,"No adjustment, correction or selection"),
                           c(xt,"Correction for antihypertensive use\n(medication decreases exposure)"),
                           c(xr,"Selection on antihypertensive users"),
                           c(xc,"Adjustment for body mass index"),
                           c(xtc,"Adjustment for body mass index and\ncorrection for antihypertensive use\n(medication decreases exposure)"),
                           c(xrc,"Adjustment for body mass index and\nselection on antihypertensive users")),
                     stringsAsFactors = FALSE)

colnames(labels) <- c("exposure","exposure_long")

df <- merge(df,labels,by = c("exposure"), all.x = TRUE)

df$exposure_long <- ifelse(!is.na(df$exposure_long),paste0(df$exposure_long,"\n",df$estimate),NA)

df$exposure_long <- factor(df$exposure_long)
df$exposure_long <- factor(df$exposure_long,levels(df$exposure_long)[c(1,2,4,6,3,5)])

# Plot x-y

ggplot(df[df$outcome=="CAD",], aes(y = exp(b), x = exposure_long)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  geom_point() + 
  geom_errorbar(aes(ymin = exp(lci), ymax = exp(uci)), width = 0) +
  labs(title = "", x = paste0(toupper(substr("systolic blood pressure", 1, 1)), substr("systolic blood pressure", 2, nchar("systolic blood pressure"))," genetic associations"),
       y = "Odds ratio and 95% confidence interval for the effect of\na standard deviation change in systolic blood pressure on coronary artery disease") +
  scale_y_continuous(trans = "log", breaks = c(1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6)) +
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
        legend.box.background = element_rect(colour = "dark grey"),
        strip.text = element_blank()) +
  coord_flip()

ggsave("output/empirical.jpeg",width = 300, height = 150, unit = "mm", dpi = 600)