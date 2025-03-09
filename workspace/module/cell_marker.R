.libPaths("/data2/platform/cell_type_workspace/RLibrary")
source("/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/module/util.R")

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE, repos = "https://mirror-hk.koddos.net/CRAN/")
  }
  library(package, character.only = TRUE)
}
# Seurat has a dependency Matrix which requires R >=4.4.0, manually install here
# install_if_missing("spatstat.utils")
# install_if_missing("Seurat")
# install_if_missing("SeuratObject")
# install_if_missing("devtools")
# install_if_missing("anndata")
# devtools::install_github('immunogenomics/presto')

library(Seurat)
library(readxl)
library(dplyr)
library(jsonlite)
library(tibble)

args              <- commandArgs(trailingOnly = TRUE)
outputdir         <- args[1]
species           <- args[2]
tissue_class      <- args[3]
override_cluster_json <- ifelse(length(args) > 3, args[4], NULL)

# outputdir <- "/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/user_data/1734203196_5900/result"
# species <- "Mouse"
# tissue_class <- "Brain"
# override_cluster_json <- "/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/user_data/1734203196_5900/result/override_cluster.json"

important_genes <- c("CD4", "CD45", "CD44", "CD3e", "CD20", "CD68", "CD11c", "CD31", "CD8", "CD45R", "B220", "Pan-Cytokeratin",
                     "Caveolin", "Ter19", "COL1A", "Vimentin", "S100A9", "IBA1", "CD206", "MRC1", "F4/80", "Ly6g", "CD36", "FcRγ",
                     "FOXP3", "Ki67", "CD45RO", "CD8", "HLA-A", "HLA-DR", "CD14", "CD56", "IFNG", "IDO1", "PD-1", "ICOS", "PD-L1",
                     "VISTA", "LAG3", "E-cadherin", "SMA", "CD34", "Beta-actin", "Podoplanin", "Collagen IV", "β-Catenin", "CD107",
                     "CD57", "CD21", "CD38", "Granzyme B", "TOX", "TCF-1", "CD79a", "CD39", "iNOS", "CD66", "MPO", "CD163", "CD11b",
                     "CD209", "Histone H3 Phospho (Ser28)", "PCNA", "Keratin-14", "Keratin-8", "Keratin-18", "Keratin-5", "SOX2",
                     "PanCK", "ER", "Bcl-2", "EpCAM", "CP100", "TP63")

obj <- readRDS(file.path(outputdir, "aKNNO_clustered.rds"))
# head(obj@meta.data)

if (!is.null(override_cluster_json)) {
  override_cluster_result <- jsonlite::fromJSON(override_cluster_json)
  for (barcode in names(override_cluster_result)) {
    obj$aKNN_O_res.0.8[barcode] <- override_cluster_result[barcode]
    obj$seurat_clusters[barcode] <- override_cluster_result[barcode]
  }
  # If the override cluster json is provided, save the updated cluster result to update gexf
  names(color_aknno) <- sort(as.integer(levels(obj$aKNN_O_res.0.8)))
  cluster_result <- data.frame(
    Cell_id = colnames(obj),
    Cluster = obj$aKNN_O_res.0.8,
    Color = color_aknno[as.character(obj$aKNN_O_res.0.8)])
  write.csv(cluster_result, file.path(outputdir, "aKNNO_cluster_result_overridden.csv"), row.names = FALSE)
}

# Find markers
Idents(obj) <- "aKNN_O_res.0.8"
min_pct_threshold <- 0.75
all_markers <- Seurat::FindAllMarkers(object = obj, logfc.threshold = 1, min.pct=min_pct_threshold, only.pos = TRUE, return.thresh=0.01)
all_markers <- all_markers %>% filter(pct.1 >= min_pct_threshold)
write.csv(all_markers, file.path(outputdir, "marker_genes.csv"), row.names = FALSE)
# all_markers <- read.csv(file.path(outputdir, "marker_genes.csv"))

# Read predefined markers and group by cell name, count the number of markers in each cell type
predefined_markers <- readxl::read_excel("/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/module/Cell_marker_All.xlsx")
predefined_markers <- predefined_markers[predefined_markers$species == species & predefined_markers$tissue_class == tissue_class, ]
predefined_marker_cell_count <- predefined_markers %>%
  group_by(cell_name) %>%
  summarise(marker_count = n())

# For each cluster, calculate the echarts option data
cluster_list <- unique(obj$aKNN_O_res.0.8)
dir.create(file.path(outputdir, "echarts_data_heatmap"), showWarnings = FALSE)
dir.create(file.path(outputdir, "marker_expression"), showWarnings = FALSE)
for(i in seq_along(cluster_list)) {
  cluster <- cluster_list[i]
  cluster_markers <- all_markers %>% filter(cluster == !!cluster)

  # Filter the predefined markers by the markers found in the clustering
  filtered_predefined_markers <- predefined_markers %>% filter(tolower(marker) %in% tolower(cluster_markers$gene))
  filtered_predefined_markers_cell_count <- filtered_predefined_markers %>%
    group_by(cell_name) %>%
    summarise(marker_count = n())

  # Calculate the proportion of the target markers in the total markers
  target_marker_genes_proportion <- filtered_predefined_markers_cell_count %>%
    left_join(predefined_marker_cell_count, by = "cell_name", suffix = c("_target", "_total")) %>%
    mutate(
      proportion = marker_count_target / marker_count_total,
      is_important = ifelse(marker %in% important_genes, 1, 0)
    ) %>%
    arrange(desc(is_important), desc(marker_count_target)) %>%
    filter(!is.na(proportion))

  # Calculate Echarts data
  x_headers <- target_marker_genes_proportion$cell_name
  y_headers <- filtered_predefined_markers %>%
    group_by(marker) %>%
    summarise(cell_name_count = n()) %>%
    arrange(desc(cell_name_count)) %>%
    pull(marker)

  # get expression data for each marker and write to separate files
  for (marker in y_headers) {
    # if marker output file exists, skip
    json_path <- file.path(outputdir, "marker_expression", paste0(marker, ".json"))
    if (file.exists(json_path) && file.size(json_path) > 0) {
      next
    }
    marker_expression <- as.list(Seurat::GetAssayData(obj, slot = "data")[marker, ])
    jsonlite::write_json(marker_expression, json_path, pretty = FALSE, auto_unbox = TRUE)
  }

  # get data
  data <- lapply(seq_along(target_marker_genes_proportion$proportion), function(i) c(i-1, 0, target_marker_genes_proportion$proportion[i]))
  for (i in seq_along(y_headers)) {
    for (j in seq_along(x_headers)) {
      if (any(filtered_predefined_markers$cell_name == x_headers[j] & filtered_predefined_markers$marker == y_headers[i])) {
        data <- append(data, list(c(j-1, i, 1)))
      } else {
        data <- append(data, list(c(j-1, i, "-")))
      }
    }
  }

  json_data <- list(
    x_headers = x_headers,
    y_headers = c("ALL", y_headers),
    data = data
  )
  jsonlite::write_json(json_data, file.path(outputdir, "echarts_data_heatmap", paste0("cluster_", cluster, "_marker_genes.json")), pretty = TRUE)
}
