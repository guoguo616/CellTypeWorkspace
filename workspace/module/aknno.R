.libPaths("/data2/platform/cell_type_workspace/RLibrary")

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

if (!require("aKNNO", character.only = TRUE)) {
  devtools::install_github("liuqivandy/aKNNO")
}
# library(Seurat)
library(aKNNO)
library(RColorBrewer)
library(readxl)
library(dplyr)

args         <- commandArgs(trailingOnly = TRUE)
inputfile    <- args[1]
outputdir    <- args[2]
species      <- args[3]
tissue_class <- args[4]

# outputdir <- "/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/user_data/1731495626_7121/result"
# inputfile <- "/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/user_data/1731495626_7121/result/visium_data"
# species <- "Mouse"
# tissue_class <- "Brain"
# if the input file is a zip, unzip it, and load using Load10X_Spatial
# if (tools::file_ext(inputfile) == "zip") {
#   tempdir <- tempdir()
#   unzip(inputfile, exdir = tempdir)
#   inputfile <- file.path(tempdir, tools::file_path_sans_ext(basename(inputfile)))
#   obj <- Seurat::Load10X_Spatial(data.dir = inputfile)
# } else if (tools::file_ext(inputfile) == "rds") {
#   obj <- readRDS(inputfile)
# } else if (tools::file_ext(inputfile) == "h5seurat") {
#   obj <- SeuratDisk::LoadH5Seurat(inputfile, meta.data = FALSE, misc = FALSE)
# } else if (tools::file_ext(inputfile) == "h5") {
#   obj <- Seurat::Read10X_h5(inputfile)
# } else {
#   # raise an error if the file extension is not supported
#     stop(paste("Unsupported file extension", tools::file_ext(inputfile)))
# }

# input_filename <- basename(inputfile)
# sample_name <- tools::file_path_sans_ext(input_filename)
obj <- Seurat::Load10X_Spatial(data.dir = inputfile)

# clustering on the optimized adaptive k-nearest neighbor graph
obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, verbose = FALSE)
obj <- aKNNO::FindNeighbors_aKNNO(obj, verbose = FALSE)
obj <- Seurat::FindClusters(obj, graph.name = "aKNN_O", verbose = FALSE)
saveRDS(obj, file.path(outputdir, "aKNNO_clustered.rds"))
# obj <- readRDS(file.path(outputdir, "aKNNO_clustered.rds"))
# head(obj@meta.data)

# define colors
color_aknno <- c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"), brewer.pal(8, "Pastel2"),
                 brewer.pal(8, "Accent"))[-c(6, 8, 10, 17, 20)]
color_aknno[c(22, 26, 27, 31, 28, 23)] <- c("#2171b5", "#FF7F00", "#238b45", "#E41A1C", "#dd3497", "#984EA3")
color_aknno[c(3, 12, 11, 24, 2, 30)] <- c(color_aknno[c(11, 3)], "#6a51a3", "#cc4c02", "#E6AB02", "#525252")
names(color_aknno) <- sort(as.integer(levels(obj$aKNN_O_res.0.8)))

p_aKNNO <- SpatialDimPlot(obj, group.by = "aKNN_O_res.0.8", cols = color_aknno) +
  ggtitle("aKNNO") +
  NoLegend()
ggsave(file.path(outputdir, "aKNNO.png"), plot = p_aKNNO)

# Output the cluster result and its color to a CSV file
cluster_result <- data.frame(Cell_id = colnames(obj), Cluster = obj$aKNN_O_res.0.8, Color = color_aknno[as.character(obj$aKNN_O_res.0.8)])
write.csv(cluster_result, file.path(outputdir, "aKNNO_cluster_result.csv"), row.names = FALSE)
# cluster_result <- read.csv(file.path(outputdir, "aKNNO_cluster_result.csv"))

# Find markers
Idents(obj) <- "aKNN_O_res.0.8"
all_markers <- FindAllMarkers(object = obj, logfc.threshold = 1, min.pct=0.25, only.pos = TRUE)
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
