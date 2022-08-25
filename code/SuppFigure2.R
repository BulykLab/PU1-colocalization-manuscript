library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(forcats, quietly = TRUE)
library(patchwork, quietly = TRUE)


#### Supp Figure 3 - JLIM vs Coloc scatterplot ####

# Reading JLIM and Coloc results
jlim_v_coloc <- read.csv('../data/Fig2/JLIM_COLOC.all.stat.txt', header = T, sep ='\t')

## FDR threshold for JLIM statistics
# (0.01172 * 1621) / sum(jlim_v_coloc$jlimp < 0.01172) ~ 0.05
# Therefore, JLIM p value threshold 0.01172 will be used for FDR 5%
jlim_threshold <- 0.01172


## Supp Figure 3
## JLIM X Coloc O
p_supp_2a_right <- ggplot(jlim_v_coloc) + theme_classic() +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_vline(xintercept = -log10(jlim_threshold), linetype = 2) +
  geom_rect(xmin=-0.1, xmax=-log10(jlim_threshold), ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0.01) +
  geom_point(aes(x=-log10(jlim_p), y = coloc_H4, color = status), alpha=1) +
  geom_rect(xmin=-0.1, xmax=-log10(jlim_threshold), ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0, color='black', size=1.1) +
  scale_color_manual("Significance", values = c('red', 'orange', 'skyblue', 'gray'), labels=c("Both", "Coloc only", "JLIM only", "Neither")) +
  xlab(expression(paste(-log[10],"(JLIM ",italic(p),")"))) + ylab("Coloc PP(Colocalization)") +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text = element_text(size=14)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .04))) + scale_x_continuous(expand = expansion(mult = c(.04, .04))) +
  theme(aspect.ratio=1)

## JLIM O Coloc X
p_supp_2b_right <- ggplot(jlim_v_coloc) + theme_classic() +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_vline(xintercept = -log10(jlim_threshold), linetype = 2) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.15, ymin=-0.035, ymax=0.5, fill="#FFF4C0", alpha=0.01) +
  geom_point(aes(x=-log10(jlim_p), y = coloc_H4, color = status), alpha=1) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.15, ymin=-0.035, ymax=0.5, fill="#FFF4C0", alpha=0, color='black', size=1.1) +
  scale_color_manual("Significance", values = c('red', 'orange', 'skyblue', 'gray'), labels=c("Both", "Coloc only", "JLIM only", "Neither")) +
  xlab(expression(paste(-log[10],"(JLIM ",italic(p),")"))) + ylab("Coloc PP(Colocalization)") +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text = element_text(size=14)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .04))) + scale_x_continuous(expand = expansion(mult = c(.04, .04))) +
  theme(aspect.ratio=1)


## Association plots
## 1
pu1_stat <- read.csv('../data/SuppFig2/PU1_17402.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_17402.mono_count.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1 <- pu1_stat %>% filter((start_var > start - 500000) & (start_var < end + 500000)) %>%
        ggplot(aes(x=start_var, y=-log10(pval))) +
          geom_vline(xintercept = start, color = 'red') +
          geom_point() +
          theme_classic() +
          labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
          scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
          scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 500000, end + 500000)) +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
                axis.text.x = element_blank(), axis.text.y = element_text(size=14),
                plot.title = element_text(
                  hjust = 1e-2,
                  margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Mono count GWAS"
p_gwas <- gwas_stat %>% filter((Position > start - 500000) & (Position < end + 500000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05)), breaks = c(0,3,6,9)) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

p_pu1_17402 <- p_pu1 + p_gwas + plot_layout(nrow = 2)

## 2
pu1_stat <- read.csv('../data/SuppFig2/PU1_27131.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_27131.MCV.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1 <- pu1_stat %>% filter((start_var > start - 500000) & (start_var < end + 500000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05)), breaks = c(0,3,6,9,12)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "MCV GWAS"
p_gwas <- gwas_stat %>% filter((Position > start - 500000) & (Position < end + 500000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

p_pu1_27131 <- p_pu1 + p_gwas + plot_layout(nrow = 2)

## 3
pu1_stat <- read.csv('../data/SuppFig2/PU1_54831.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_54831.lym_count.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1 <- pu1_stat %>% filter((start_var > start - 150000) & (start_var < end + 150000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05)), breaks = c(0,3,6,9)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 150000, end + 150000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Lym count GWAS"
p_gwas <- gwas_stat %>% filter((Position > start - 150000) & (Position < end + 150000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 150000, end + 150000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

p_pu1_54831 <- p_pu1 + p_gwas + plot_layout(nrow = 2)

## 4
pu1_stat <- read.csv('../data/SuppFig2/PU1_64115.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_64115.mono_count.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1 <- pu1_stat %>% filter((start_var > start - 500000) & (start_var < end + 500000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05)), breaks = c(0,3,6,9)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Mono count GWAS"
p_gwas <- gwas_stat %>% filter((Position > start - 500000) & (Position < end + 500000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

p_pu1_64115 <- p_pu1 + p_gwas + plot_layout(nrow = 2)






p_supp_2a_left <- plot_grid(p_pu1_17402, p_pu1_27131,p_pu1_54831, p_pu1_64115,ncol = 2)
p_supp_2a <- plot_grid(p_supp_2a_left, NULL, p_supp_2a_right, ncol = 3, rel_widths = c(2,0.1,1))

ggsave('../figures/Supp2_a.pdf',
       plot = p_supp_2a,
       device='pdf',
       width=400, height=140, units="mm")



##
## 1
pu1_stat <- read.csv('../data/SuppFig2/PU1_35142.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_35142.eosino_count.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1_35142 <- pu1_stat %>% filter((start_var > start - 350000) & (start_var < end + 350000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 350000, end + 350000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Eosino count GWAS"
p_gwas_35142 <- gwas_stat %>% filter((Position > start - 350000) & (Position < end + 350000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 350000, end + 350000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))


## 2
pu1_stat <- read.csv('../data/SuppFig2/PU1_47192.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_47192.mono_percentage.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1_47192 <- pu1_stat %>% filter((start_var > start - 500000) & (start_var < end + 500000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Mono % GWAS"
p_gwas_47192 <- gwas_stat %>% filter((Position > start - 500000) & (Position < end + 500000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))


## 3
pu1_stat <- read.csv('../data/SuppFig2/PU1_58189.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_58189.eosino_percentage.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1_58189 <- pu1_stat %>% filter((start_var > start - 500000) & (start_var < end + 500000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Eosino % GWAS"
p_gwas_58189 <- gwas_stat %>% filter((Position > start - 500000) & (Position < end + 500000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))


## 4
pu1_stat <- read.csv('../data/SuppFig2/PU1_65827.qtl.txt', header = T, sep =' ')
gwas_stat <- read.csv('../data/SuppFig2/PU1_65827.eosino_percentage.sumstat', header = T, sep =' ')

chr <- pu1_stat[1,2]
start <- pu1_stat[1,3]
end <- pu1_stat[1,4]

title <- "PU.1 bQTL"
p_pu1_65827 <- pu1_stat %>% filter((start_var > start - 500000) & (start_var < end + 500000)) %>%
  ggplot(aes(x=start_var, y=-log10(pval))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))

title <- "Eosino % GWAS"
p_gwas_65827 <- gwas_stat %>% filter((Position > start - 500000) & (Position < end + 500000)) %>%
  ggplot(aes(x=Position, y=-log10(PV))) +
  geom_vline(xintercept = start, color = 'red') +
  geom_point() +
  theme_classic() +
  labs(x = paste0(chr,' (Mb)'), y = expression(paste(-log[10],"(",italic(p),")")), title = title)  +
  scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)),
                     lim = c(start - 500000, end + 500000)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))


##

## Plotting together
p_supp_2b_left <- p_pu1_35142  + p_pu1_47192 + p_gwas_35142 + p_gwas_47192 +
  p_pu1_58189  + p_pu1_65827 + p_gwas_58189 + p_gwas_65827 + plot_layout(nrow=4, ncol=2)

p_supp_2b <- plot_grid(p_supp_2b_left, NULL, p_supp_2b_right, ncol = 3, rel_widths = c(2,0.1,1))



ggsave('../figures/Supp2_b.pdf',
       plot = p_supp_2b,
       device='pdf',
       width=400, height=140, units="mm")
