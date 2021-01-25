remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Load data

df <- data.table::fread("output/simulation.csv", 
                        data.table = FALSE,
                        select = c("model","z1_mb","_LCI","_UCI","_WT"),
                        stringsAsFactors = FALSE)

df <- df[df$model!="",]
colnames(df) <- c("model","meta","lci","uci","weight")

# Make percentages

df$meta <- 100*df$meta
df$lci <- 100*df$lci
df$uci <- 100*df$uci

df$estimate <- paste0(sprintf("%.0f",df$meta),"% (95% CI: ",sprintf("%.0f",df$lci)," to ",sprintf("%.0f",df$uci),")")

# Label exposures 

labels <- data.frame(rbind(c("z0_mb","Model 0\nNo adjustment, correction or selection"),
                           c("z1_mb","Model 1\nAdjustment for a covariate"),
                           c("z2_mb","Model 2\nCorrection for medication use"),
                           c("z3_mb","Model 3\nSelection on medication users"),
                           c("z4_mb","Model 4\nAdjustment for a covariate and\ncorrection for medication use"),
                           c("z5_mb","Model 5\nAdjustment for a covariate and\nselection on medication users")),
                     stringsAsFactors = FALSE)

colnames(labels) <- c("model","model_long")

df <- merge(df,labels,by = c("model"), all.x = TRUE)

df$model_long <- ifelse(!is.na(df$model_long),paste0(df$model_long,"\n",df$estimate),NA)



# Plot 

ggplot2::ggplot(df, ggplot2::aes(y = meta, x = forcats::fct_rev(model_long))) +
  ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "darkgray") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = meta-0.1, ymax = meta+0.1), alpha = 1, size = 2, color = "black") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.5, size = 2, color = "darkgray") +
  ggplot2::scale_y_continuous(lim = c(-10,70), breaks = seq(-10,70,10)) +
  ggplot2::labs(title = "", 
                x = "",
                y = "Percentage difference between the estimates and the simulated true\neffect of unmeasured A on B, scaled and meta-analysed across effect sizes") +
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

ggplot2::ggsave("output/Figure3.jpg",
                width = 174, height = 120, 
                unit = "mm", dpi = 320,
                scale = 1.1)
