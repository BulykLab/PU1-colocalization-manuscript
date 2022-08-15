library(ggplot2)
library(cowplot)


legend_box = data.frame(x = 0.8, y = seq(0.8, 0.4, -0.1))


r2_panel <- ggplot() +
  geom_rect(data = legend_box,
            aes(xmin = x, xmax = x + 0.1, ymin = y, ymax = y + 0.1),
            color = "black",
            fill = rev(c("navy", "lightskyblue", "green", "orange", "red"))) +
  draw_label("1", x = legend_box$x[1] + 0.1, y = legend_box$y[1]+0.1, hjust = -0.7, size = 14) +
  draw_label("0.8", x = legend_box$x[1] + 0.1, y = legend_box$y[1], hjust = -0.3, size = 14) +
  draw_label("0.6", x = legend_box$x[2] + 0.1, y = legend_box$y[2], hjust = -0.3, size = 14) +
  draw_label("0.4", x = legend_box$x[3] + 0.1, y = legend_box$y[3], hjust = -0.3, size = 14) +
  draw_label("0.2", x = legend_box$x[4] + 0.1, y = legend_box$y[4], hjust = -0.3, size = 14) +
  draw_label("0", x = legend_box$x[5] + 0.1, y = legend_box$y[5], hjust = -0.7, size = 14) +
  draw_label(parse(text = "r^2"), x = legend_box$x[1]+0.06, y = legend_box$y[1]+0.03, vjust = -1.6, size = 14) +
  xlim(c(0.675,1.025)) + ylim(c(0.3,1)) +
  theme_void() + theme(aspect.ratio=2)



ggplot2::ggsave('../figures/r2_panel.pdf',
                plot = r2_panel,
                device='pdf',
                width=70, height=70, units="mm")
