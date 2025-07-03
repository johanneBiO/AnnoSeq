library(ggplot2)

main_theme <- theme_classic() +
  theme(axis.title.x = element_text(face = "bold", vjust = -1),
        axis.title.y = element_text(face = "bold", vjust = 3),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 11),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8))
