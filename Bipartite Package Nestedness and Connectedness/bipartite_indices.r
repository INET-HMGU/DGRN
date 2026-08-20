# bipartite analysis of Y1H under plant drought stress stimulation
# 21.02.2023
# Lin Chung-wen

# environment setup
data_dir <- "~/Documents/INET-work/Y1H_promoter/data"
result_dir <- "~/Documents/INET-work/Y1H_promoter/result"
target_dir <- "/Volumes/INET/Daten/AG\ Falter-Braun/Ram/CW_Network_Analysis"

library(openxlsx)
library(bipartite)
library(igraph)

# load data
edge <- read.xlsx(file.path(data_dir, "DGRN_2.xlsx"))

# network and adjacency matrix establishment
species <- edge$Species

networks <- list()
networks_matrix <- list()

for (s in unique(species)) {
    df <- edge[edge$Species == s, c(5, 3)]
    net <- graph_from_data_frame(df, directed = T)
    V(net)$type <- V(net)$name %in% df$Orthogroup

    net_matrix <- as.matrix(as.matrix(net))
    rows <- unique(df$TF_Common.name)
    cols <- unique(df$Orthogroup)
    net_matrix <- net_matrix[rows, colnames(net_matrix) %in% cols]
    networks_matrix[[s]] <- net_matrix
    networks[[s]] <- net
}

# indices calculation
network_indices <- lapply(networks_matrix, networklevel)
group_indices <- lapply(networks_matrix, grouplevel)
node_indices <- lapply(networks_matrix, specieslevel)

write.csv(t(do.call(rbind, network_indices)), file.path(target_dir, "networkLevel_indices.csv"), row.names = TRUE)

wb <- createWorkbook()

addWorksheet(wb, "Ath_higherlevel")
writeData(wb, sheet = 1, node_indices[["Ath"]]$"higher level", rowNames = TRUE)
addWorksheet(wb, "Ath_lowerlevel")
writeData(wb, sheet = 2, node_indices[["Ath"]]$"lower level", rowNames = TRUE)

addWorksheet(wb, "Aly_higherlevel")
writeData(wb, sheet = 3, node_indices[["Aly"]]$"higher level", rowNames = TRUE)
addWorksheet(wb, "Aly_lowerlevel")
writeData(wb, sheet = 4, node_indices[["Aly"]]$"lower level", rowNames = TRUE)

addWorksheet(wb, "Esa_higherlevel")
writeData(wb, sheet = 5, node_indices[["Esa"]]$"higher level", rowNames = TRUE)
addWorksheet(wb, "Esa_lowerlevel")
writeData(wb, sheet = 6, node_indices[["Esa"]]$"lower level", rowNames = TRUE)

saveWorkbook(wb, file.path(target_dir, "nodeLevel_indices.xlsx"), overwrite = TRUE)

# module detection
modules <- lapply(networks_matrix, computeModules)

pdf(file.path(target_dir, "bipartite_module.pdf"), width = 5, height =  8)
for (s in unique(species)) {
    plotModuleWeb(modules[[s]], labsize = 0.5)
    title(main = s, line = -1)
    mtext("TF", side = 2, line = 2.5, cex = 0.7)
    mtext("Orthologue", side = 1, line = 0.5, cex = 0.7)
}
dev.off()

# save data
save.image(file = file.path(result_dir, "bipartite.RData"))