library(ggplot2)

main_theme <- theme_classic() +
  theme(axis.title.x = element_text(face = "bold", vjust = -1),
        axis.title.y = element_text(face = "bold", vjust = 3),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 11),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8))

palette_19 <- c(
  "#2980b9", "#1f618d", "#B7E4C7", "#1B9E77", "#117864",
  "#DDD6B8", "#CBB67C", "#5C4033",
  "#990000", "#FF0000", "#d35400", "#f39c12", "#F0E442", "black", "#999999", "#cccccc",
  "#E7298A", "#CC79A7", "#7d3c98")

palette_26 <- c(
  "#2980b9", "#1f618d", "blue", "#00CED1", 
  "#B7E4C7", "#1B9E77", "#117864", "green", 
  "#DDD6B8", "#CBB67C", "#8B4513", "#5C4033",
  "#990000", "#FF0000", "#d35400", "#f39c12", "#FFD700" , "#F0E442", 
  "black", "#999999", "#cccccc", 
  "pink", "#E7298A", "#CC79A7", "#7d3c98", "#4B0082")