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

df <- data.table::fread("output/empirical.csv", data.table = FALSE)

df$estimate <- paste0(sprintf("%.2f",exp(df$b))," (95% CI: ",sprintf("%.2f",exp(df$lci))," to ",sprintf("%.2f",exp(df$uci)),"; #SNPs = ",df$nsnp,")")

# Label exposures 

labels <- data.frame(rbind(c(x,"Model 0\nNo adjustment, correction or selection"),
                           c(xt,"Model 3\nCorrection for antihypertensive use"),
                           c(xr,"Model 2\nSelection on antihypertensive users"),
                           c(xc,"Model 1\nAdjustment for body mass index"),
                           c(xtc,"Model 5\nAdjustment for body mass index and\ncorrection for antihypertensive use"),
                           c(xrc,"Model 4\nAdjustment for body mass index and\nselection on antihypertensive users")),
                     stringsAsFactors = FALSE)

colnames(labels) <- c("exposure","exposure_long")

df <- merge(df,labels,by = c("exposure"), all.x = TRUE)

df$exposure_long <- ifelse(!is.na(df$exposure_long),paste0(df$exposure_long,"\n",df$estimate),NA)

# Plot

ggplot2::ggplot(df, ggplot2::aes(y = exp(b), x = forcats::fct_rev(exposure_long))) +
  ggplot2::geom_hline(yintercept = 1, linetype = "solid", color = "darkgray") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = exp(b)-1e-3, ymax = exp(b)+1e-3), alpha = 1, size = 2, color = "darkgray") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = exp(lci), ymax = exp(uci)), alpha = 0.5, size = 2, color = "darkgray") +
  ggplot2::labs(title = "", 
                x = "",
                y = "Odds ratio and 95% confidence interval for the effect of a standard\ndeviation increase in systolic blood pressure on coronary artery disease") +
  ggplot2::scale_y_continuous(trans = "log", breaks = c(1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6)) +
  ggplot2::guides(col = ggplot2::guide_legend(ncol=1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text = ggplot2::element_text(size=10),
        axis.ticks = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        text = ggplot2::element_text(size=10),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size=8),
        legend.position = "bottom",
        legend.key.width = ggplot2::unit(0.1,"cm"),
        legend.box.background = ggplot2::element_rect(colour = "dark grey"),
        strip.text = ggplot2::element_blank()) +
  ggplot2::coord_flip()

ggplot2::ggsave("output/empirical.jpeg",
                width = 300, height = 150,
                unit = "mm", dpi = 600)