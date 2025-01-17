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
library(RColorBrewer)
library(readxl)
library(dplyr)
library(jsonlite)

args              <- commandArgs(trailingOnly = TRUE)
outputdir         <- args[1]
species           <- args[2]
tissue_class      <- args[3]
override_cluster_json <- ifelse(length(args) > 3, args[4], NULL)

# outputdir <- "/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/user_data/1731495626_7121/result"
# species <- "Mouse"
# tissue_class <- "Brain"

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
  write.csv(cluster_result, file.path(outputdir, "aKNNO_cluster_result.csv"), row.names = FALSE)
}

# Find markers
Idents(obj) <- "aKNN_O_res.0.8"
all_markers <- Seurat::FindAllMarkers(object = obj, logfc.threshold = 1, min.pct=0.25, only.pos = TRUE)
all_markers <- all_markers %>% filter(p_val_adj < 0.01)
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
for(i in seq_along(cluster_list)) {
  cluster <- cluster_list[i]
  cluster_markers <- all_markers %>% filter(cluster == !!cluster)

  # Filter the predefined markers by the markers found in the clustering
  filtered_predefined_markers <- predefined_markers %>% filter(marker %in% cluster_markers$gene)
  filtered_predefined_markers_cell_count <- filtered_predefined_markers %>%
    group_by(cell_name) %>%
    summarise(marker_count = n())

  # Calculate the proportion of the target markers in the total markers
  target_marker_genes_proportion <- filtered_predefined_markers_cell_count %>%
    left_join(predefined_marker_cell_count, by = "cell_name", suffix = c("_target", "_total")) %>%
    mutate(proportion = marker_count_target / marker_count_total) %>%
    arrange(desc(proportion)) %>%
    filter(!is.na(proportion))

  # Calculate Echarts data
  x_headers <- target_marker_genes_proportion$cell_name
  y_headers <- filtered_predefined_markers %>%
    group_by(marker) %>%
    summarise(cell_name_count = n()) %>%
    arrange(desc(cell_name_count)) %>%
    pull(marker)

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
