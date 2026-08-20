# Drought stress network establishment, as bipartite model
# Lin Chung-wen
# 07.02.2023 13:53

# environment
library(igraph)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(gridExtra)

data_dir <- "/Volumes/INET/Daten/AG Falter-Braun/Ram/CW_Network_Analysis"

species <- c("Ath", "Aly", "Esa")

# load data
raw_PDI <- read.xlsx(file.path(data_dir, "DGRN_2.xlsx"), sheet = 1)

# network graph establishment and degree calculation
net <- list()
net_adj <- list()

for (n in species) {
    # network graph
    net[[n]] <- graph_from_data_frame(raw_PDI[raw_PDI$Species == n, c("Promoter", "TF_Common.name")], directed = FALSE)
    V(net[[n]])$type <- V(net[[n]])$name %in% raw_PDI[raw_PDI$Species == n, "TF_Common.name"]

    # degree detection
    tf <- unique(raw_PDI[raw_PDI$Species == n, "TF_Common.name"]) # as column
    promoter <- unique(raw_PDI[raw_PDI$Species == n, "Promoter"]) # as raw
    adj_tmp <- as_adjacency_matrix(net[[n]])
    net_adj[[n]] <- as.matrix(adj_tmp[promoter, tf])
}

write.xlsx(net_adj, file = file.path(data_dir, "network_as_adjacency_matrix.xlsx"), 
    sheetName = species, rowNames = TRUE, overwrite = TRUE)

# degree distribution calculation
## promoter level
promoter_deg_c <- c()
tf_deg_c <- c()

for (n in species) {
    pro_deg_tmp <- as.numeric(rowSums(net_adj[[n]]))
    promoter_deg_c <- c(promoter_deg, pro_deg_tmp)
    tf_deg_tmp <- as.numeric(colSums(net_adj[[n]]))
    tf_deg <- c(tf_deg, tf_deg_tmp)
}

promoter_deg <- data.frame(
    degree = as.numeric(promoter_deg),
    species = rep(species, sapply(net_adj, nrow))
    )
tf_deg <- data.frame(
    degree = tf_deg,
    species = rep(species, sapply(net_adj, ncol))
    )

# statistical analysis
## check degree is normal distribution or not

pdf(file.path(data_dir, "degree_distribution.pdf"), height = 8)
par(mfrow = c(3, 2))
for (n in species) {
    stest <- shapiro.test(promoter_deg[promoter_deg$species == n, "degree"])
    message("Based on Shapiro test, degree of promoter of species ", n, " with p value = \n", stest$p.value, "\n")
    hist(promoter_deg[promoter_deg$species == n, "degree"], main = paste0("Promoter degree of ", n), xlab = "degree")
    mtext(side = 3, line = -0.3, cex = 0.5, text = paste0("Shapiro test P = ", formatC(stest$p.value, format = "e", digits = 2), "\n", ifelse(stest$p.value < 0.05, "match non-normal distribution", "match normal distribution")))

    stest <- shapiro.test(tf_deg[tf_deg$species == n, "degree"])
    message("Based on Shapiro test, degree of TF of species ", n, " with p value = \n", stest$p.value, "\n")
    hist(tf_deg[tf_deg$species == n, "degree"], main = paste0("TF degree of ", n), xlab = "degree")
    mtext(side = 3, line = -0.3, cex = 0.5, text = paste0("Shapiro test P = ", formatC(stest$p.value, format = "e", digits = 2), "\n", ifelse(stest$p.value < 0.05, "match non-normal distribution", "match normal distribution")))
}
dev.off()

# > some degree distributions are matched normal distribution, some are not!
# > stay with Wilcox test (U test), which does not make an assumption about the distribution of a dataset.

# data visulization
theme_gp <- theme_pubr() +
    theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        # axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))

jitter_width <- 0.1
jitter_alpha <- 0.2
promoter_gp <- ggplot(promoter_deg, aes(x = species, y = degree, fill = species)) + 
    # geom_violin() +
    geom_boxplot(width = 0.5) + 
    geom_jitter(width = jitter_width, color = "red", alpha = jitter_alpha) +
    labs(x = "", y = "Degree", title = "Promoter degree distribution") +
    theme_gp
tf_gp <- ggplot(tf_deg, aes(x = species, y = degree, fill = species)) + 
    # geom_violin() +
    geom_boxplot(width = 0.5) + 
    geom_jitter(width = jitter_width, color = "red", alpha = jitter_alpha) +
    labs(x = "", y = "Degree", title = "TF degree distribution") +
    theme_gp

my_comparisons <- list(c("Ath", "Aly"), c("Ath", "Esa"), c("Aly", "Esa"))

# add wilcoxon test, by default of 'stat_compare_means'
to_plot <- list(
    promoter = promoter_gp + stat_compare_means(aes(label = ..p.format..), comparisons = my_comparisons),
    tf = tf_gp + stat_compare_means(aes(label = ..p.format..), comparisons = my_comparisons) 
)

ggsave(file.path(data_dir, "degree_with_statistics.pdf"),
    marrangeGrob(grobs = to_plot, nrow = 2, ncol = 1, top = NULL),
    device = "pdf", width = 3, height = 8
)