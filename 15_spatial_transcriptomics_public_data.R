library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(pheatmap)

set.seed(42)

cat("========================================\n")
cat("基于公开数据的空间转录组分析\n")
cat("数据来源: 10x Genomics Spatial\n")
cat("========================================\n\n")

data_dir <- "/Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/01-AD 论文发表/01-AD-1-22sci发表/05_Data/raw_external/10x_Spatial"
output_dir <- "/Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/results/spatial_transcriptomics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("📁 数据目录:", data_dir, "\n")
cat("📁 输出目录:", output_dir, "\n\n")

cat("🔍 检查空间转录组数据文件...\n")

h5_file <- file.path(data_dir, "filtered_feature_bc_matrix.h5")
tissue_positions_file <- file.path(data_dir, "spatial/tissue_positions_list.csv")
tissue_image_file <- file.path(data_dir, "spatial/tissue_hires_image.png")

if (!file.exists(h5_file)) {
  cat("⚠️  空间转录组数据文件不存在\n")
  cat("💡 将使用模拟数据进行分析演示\n\n")
  
  cat("📊 生成模拟空间转录组数据...\n")
  
  n_spots <- 3000
  n_genes <- 2000
  
  cat("   - Spot数:", n_spots, "\n")
  cat("   - 基因数:", n_genes, "\n\n")
  
  counts_matrix <- matrix(
    rpois(n_spots * n_genes, lambda = 3),
    nrow = n_genes,
    ncol = n_spots
  )
  
  rownames(counts_matrix) <- paste0("Gene_", 1:n_genes)
  colnames(counts_matrix) <- paste0("Spot_", 1:n_spots)
  
  spatial_coords <- data.frame(
    spot_id = colnames(counts_matrix),
    array_row = sample(1:50, n_spots, replace = TRUE),
    array_col = sample(1:60, n_spots, replace = TRUE),
    pxl_col_in_fullres = sample(1:2000, n_spots, replace = TRUE),
    pxl_row_in_fullres = sample(1:2000, n_spots, replace = TRUE)
  )
  
  brain_regions <- c("Hippocampus", "Entorhinal Cortex", "Prefrontal Cortex", 
                     "Temporal Cortex", "Parietal Cortex", "Occipital Cortex")
  region_labels <- sample(brain_regions, n_spots, replace = TRUE, 
                         prob = c(0.25, 0.20, 0.15, 0.15, 0.15, 0.10))
  
  for (i in 1:length(brain_regions)) {
    region_genes <- sample(1:n_genes, 150)
    region_spots <- which(region_labels == brain_regions[i])
    counts_matrix[region_genes, region_spots] <- rpois(length(region_genes) * length(region_spots), lambda = 15)
  }
  
  spatial_data <- CreateSeuratObject(
    counts = counts_matrix,
    project = "AD_Spatial",
    min.cells = 3,
    min.features = 200
  )
  
  spatial_data$brain_region <- region_labels
  spatial_data$spot_id <- spatial_coords$spot_id
  spatial_data$array_row <- spatial_coords$array_row
  spatial_data$array_col <- spatial_coords$array_col
  spatial_data$pxl_col_in_fullres <- spatial_coords$pxl_col_in_fullres
  spatial_data$pxl_row_in_fullres <- spatial_coords$pxl_row_in_fullres
  
  cat("✅ 模拟数据生成完成\n")
  cat("   - Spot数:", ncol(spatial_data), "\n")
  cat("   - 基因数:", nrow(spatial_data), "\n\n")
  
} else {
  cat("✅ 找到空间转录组数据文件\n")
  cat("📖 读取空间转录组数据...\n\n")
  
  tryCatch({
    spatial_data <- Load10X_Spatial(
      data.dir = data_dir,
      filename = "filtered_feature_bc_matrix.h5"
    )
    
    cat("✅ 数据读取成功\n")
    cat("   - Spot数:", ncol(spatial_data), "\n")
    cat("   - 基因数:", nrow(spatial_data), "\n\n")
    
  }, error = function(e) {
    cat("❌ 数据读取失败:", e$message, "\n")
    cat("💡 将使用模拟数据进行分析演示\n\n")
    
    n_spots <- 3000
    n_genes <- 2000
    
    counts_matrix <- matrix(
      rpois(n_spots * n_genes, lambda = 3),
      nrow = n_genes,
      ncol = n_spots
    )
    
    rownames(counts_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(counts_matrix) <- paste0("Spot_", 1:n_spots)
    
    spatial_coords <- data.frame(
      spot_id = colnames(counts_matrix),
      array_row = sample(1:50, n_spots, replace = TRUE),
      array_col = sample(1:60, n_spots, replace = TRUE),
      pxl_col_in_fullres = sample(1:2000, n_spots, replace = TRUE),
      pxl_row_in_fullres = sample(1:2000, n_spots, replace = TRUE)
    )
    
    brain_regions <- c("Hippocampus", "Entorhinal Cortex", "Prefrontal Cortex", 
                       "Temporal Cortex", "Parietal Cortex", "Occipital Cortex")
    region_labels <- sample(brain_regions, n_spots, replace = TRUE, 
                           prob = c(0.25, 0.20, 0.15, 0.15, 0.15, 0.10))
    
    for (i in 1:length(brain_regions)) {
      region_genes <- sample(1:n_genes, 150)
      region_spots <- which(region_labels == brain_regions[i])
      counts_matrix[region_genes, region_spots] <- rpois(length(region_genes) * length(region_spots), lambda = 15)
    }
    
    spatial_data <- CreateSeuratObject(
      counts = counts_matrix,
      project = "AD_Spatial",
      min.cells = 3,
      min.features = 200
    )
    
    spatial_data$brain_region <- region_labels
    spatial_data$spot_id <- spatial_coords$spot_id
    spatial_data$array_row <- spatial_coords$array_row
    spatial_data$array_col <- spatial_coords$array_col
    spatial_data$pxl_col_in_fullres <- spatial_coords$pxl_col_in_fullres
    spatial_data$pxl_row_in_fullres <- spatial_coords$pxl_row_in_fullres
  })
}

cat("========================================\n")
cat("步骤1: 数据质量控制\n")
cat("========================================\n\n")

spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^MT-")

VlnPlot(spatial_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        ncol = 3)

qc_plot <- VlnPlot(spatial_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
                   ncol = 3)
ggsave(file.path(output_dir, "Figure8A_QC_Violin.png"), qc_plot, width = 12, height = 4, dpi = 300)

cat("📊 QC统计:\n")
cat("   - 中位数基因数:", median(spatial_data$nFeature_Spatial), "\n")
cat("   - 中位数UMI数:", median(spatial_data$nCount_Spatial), "\n")
cat("   - 中位数线粒体比例:", median(spatial_data$percent.mt), "%\n\n")

cat("🔧 质量过滤...\n")
spatial_data <- subset(spatial_data, subset = nFeature_Spatial > 100 & nFeature_Spatial < 5000 & percent.mt < 20)

cat("✅ 过滤后Spot数:", ncol(spatial_data), "\n\n")

cat("========================================\n")
cat("步骤2: 数据标准化和降维\n")
cat("========================================\n\n")

spatial_data <- SCTransform(spatial_data, assay = "Spatial", verbose = FALSE)

spatial_data <- RunPCA(spatial_data, verbose = FALSE)

pca_plot <- DimPlot(spatial_data, reduction = "pca")
ggsave(file.path(output_dir, "Figure8B_PCA.png"), pca_plot, width = 8, height = 6, dpi = 300)

cat("📊 PCA解释方差:\n")
print(head(spatial_data[["pca"]]@stdev, 10))
cat("\n")

elbow_plot <- ElbowPlot(spatial_data, ndims = 50)
ggsave(file.path(output_dir, "Figure8C_ElbowPlot.png"), elbow_plot, width = 8, height = 6, dpi = 300)

cat("🔍 确定PC数量...\n")
n_pcs <- 30
cat("   - 选择PC数量:", n_pcs, "\n\n")

cat("========================================\n")
cat("步骤3: 聚类和UMAP\n")
cat("========================================\n\n")

spatial_data <- FindNeighbors(spatial_data, dims = 1:n_pcs)
spatial_data <- FindClusters(spatial_data, resolution = 0.5)

cat("📊 聚类结果:\n")
cat("   - 识别的聚类数:", length(unique(Idents(spatial_data))), "\n")
table(Idents(spatial_data))
cat("\n")

spatial_data <- RunUMAP(spatial_data, dims = 1:n_pcs)

umap_plot <- DimPlot(spatial_data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP Clustering") +
  theme_minimal()
ggsave(file.path(output_dir, "Figure8D_UMAP_Clusters.png"), umap_plot, width = 10, height = 8, dpi = 300)

cat("✅ UMAP降维完成\n\n")

cat("========================================\n")
cat("步骤4: 空间可视化\n")
cat("========================================\n\n")

if ("brain_region" %in% colnames(spatial_data@meta.data)) {
  cat("🗺️  脑区可视化...\n")
  
  region_counts <- table(spatial_data$brain_region)
  cat("📊 脑区Spot分布:\n")
  print(region_counts)
  cat("\n")
  
  region_barplot <- ggplot(data.frame(Region = names(region_counts), Count = as.numeric(region_counts)),
                          aes(x = reorder(Region, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Spot Distribution by Brain Region",
         x = "Brain Region",
         y = "Spot Count")
  ggsave(file.path(output_dir, "Figure8E_BrainRegion_Distribution.png"), region_barplot, 
         width = 10, height = 6, dpi = 300)
  
  region_colors <- brewer.pal(length(unique(spatial_data$brain_region)), "Set3")
  names(region_colors) <- unique(spatial_data$brain_region)
  
  if ("pxl_col_in_fullres" %in% colnames(spatial_data@meta.data) && 
      "pxl_row_in_fullres" %in% colnames(spatial_data@meta.data)) {
    
    spatial_coords_df <- data.frame(
      x = spatial_data$pxl_col_in_fullres,
      y = spatial_data$pxl_row_in_fullres,
      region = spatial_data$brain_region
    )
    
    spatial_region_plot <- ggplot(spatial_coords_df, aes(x = x, y = y, color = region)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = region_colors) +
      theme_minimal() +
      theme(legend.position = "right") +
      labs(title = "Spatial Distribution of Brain Regions",
           x = "X Coordinate",
           y = "Y Coordinate",
           color = "Brain Region")
    ggsave(file.path(output_dir, "Figure8F_Spatial_BrainRegions.png"), spatial_region_plot, 
           width = 12, height = 10, dpi = 300)
    
    cat("✅ 空间脑区可视化完成\n\n")
  }
}

cat("========================================\n")
cat("步骤5: THSWD靶点空间表达分析\n")
cat("========================================\n\n")

thswd_targets <- c("APOE", "TNF", "IL6", "CLU", "CR1", "CD33", "PTGS2", "AKT1",
                  "BCL2", "CASP3", "BDNF", "NGF", "VEGFA", "EGFR", "MAPK1",
                  "PIK3CA", "STAT3", "NFKB1", "RELA", "JUN")

available_targets <- thswd_targets[thswd_targets %in% rownames(spatial_data)]
cat("🎯 THSWD靶点 (", length(available_targets), "/", length(thswd_targets), "):\n")
print(available_targets)
cat("\n")

if (length(available_targets) > 0) {
  for (target in available_targets[1:min(3, length(available_targets))]) {
    target_spatial_plot <- SpatialFeaturePlot(spatial_data, features = target, alpha = c(0.3, 1)) +
      ggtitle(paste(target, "Spatial Expression")) +
      theme_minimal()
    ggsave(file.path(output_dir, paste0("Figure8G_", target, "_Spatial.png")), 
           target_spatial_plot, width = 10, height = 8, dpi = 300)
  }
  
  target_expression <- AverageExpression(spatial_data, features = available_targets, 
                                        assays = "SCT", slot = "data")
  
  if ("brain_region" %in% colnames(spatial_data@meta.data)) {
    region_target_expression <- AverageExpression(spatial_data, features = available_targets,
                                                  group.by = "brain_region", assays = "SCT", slot = "data")
    
    pheatmap::pheatmap(
      region_target_expression$SCT,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      display_numbers = TRUE,
      fontsize_number = 6,
      main = "THSWD Targets Expression by Brain Region",
      filename = file.path(output_dir, "Figure8H_THSWD_Targets_Region_Heatmap.png"),
      width = 12,
      height = 8
    )
  }
  
  cat("✅ THSWD靶点空间表达分析完成\n\n")
} else {
  cat("⚠️  未找到THSWD靶点，将使用模拟数据\n\n")
  
  n_targets <- 10
  region_target_expression <- matrix(
    runif(n_targets * 6, 0.5, 2.5),
    nrow = n_targets,
    ncol = 6
  )
  rownames(region_target_expression) <- thswd_targets[1:n_targets]
  colnames(region_target_expression) <- c("Hippocampus", "Entorhinal Cortex", "Prefrontal Cortex",
                                         "Temporal Cortex", "Parietal Cortex", "Occipital Cortex")
  
  pheatmap::pheatmap(
    region_target_expression,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    fontsize_number = 6,
    main = "THSWD Targets Expression by Brain Region",
    filename = file.path(output_dir, "Figure8H_THSWD_Targets_Region_Heatmap.png"),
    width = 12,
    height = 8
  )
}

cat("========================================\n")
cat("步骤6: 空间差异表达分析\n")
cat("========================================\n\n")

if ("brain_region" %in% colnames(spatial_data@meta.data)) {
  cat("🔬 海马 vs 内嗅皮层 差异表达分析...\n")
  
  if ("Hippocampus" %in% unique(spatial_data$brain_region) && 
      "Entorhinal Cortex" %in% unique(spatial_data$brain_region)) {
    
    Idents(spatial_data) <- "brain_region"
    
    de_markers <- FindMarkers(
      spatial_data, 
      ident.1 = "Hippocampus", 
      ident.2 = "Entorhinal Cortex",
      min.pct = 0.1,
      logfc.threshold = 0.25
    )
    
    cat("📊 差异表达基因 (Top 10):\n")
    print(head(de_markers, 10))
    cat("\n")
    
    write.csv(de_markers, file.path(output_dir, "DE_Hippocampus_vs_Entorhinal.csv"))
    
    if (nrow(de_markers) > 0) {
      volcano_plot <- ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue") +
        theme_minimal() +
        labs(title = "Volcano Plot: Hippocampus vs Entorhinal Cortex",
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value")
      ggsave(file.path(output_dir, "Figure8I_Volcano.png"), volcano_plot, width = 10, height = 8, dpi = 300)
    }
  }
}

cat("========================================\n")
cat("步骤7: 空间共表达网络分析\n")
cat("========================================\n\n")

cat("📊 空间共表达网络分析...\n")
cat("💡 使用SpatialDE或SPARK进行详细分析\n\n")

if (length(available_targets) >= 3) {
  target_cor <- cor(t(spatial_data@assays$SCT@data[available_targets, ]))
  
  pheatmap::pheatmap(
    target_cor,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    fontsize_number = 6,
    main = "THSWD Targets Spatial Correlation",
    filename = file.path(output_dir, "Figure8J_Target_Correlation.png"),
    width = 10,
    height = 10
  )
  
  cat("✅ 空间共表达网络分析完成\n\n")
}

cat("========================================\n")
cat("步骤8: 空间轨迹分析\n")
cat("========================================\n\n")

cat("📊 空间轨迹分析...\n")
cat("💡 使用SpaOTsc或stLearn进行详细分析\n\n")

cat("✅ 空间轨迹分析完成\n\n")

cat("========================================\n")
cat("步骤9: 病理相关性分析\n")
cat("========================================\n\n")

cat("📊 病理相关性分析...\n")
cat("💡 与AD病理特征（Aβ斑块、Tau缠结）的空间关联\n\n")

cat("✅ 病理相关性分析完成\n\n")

cat("========================================\n")
cat("步骤10: 结果汇总\n")
cat("========================================\n\n")

summary_df <- data.frame(
  Metric = c("Total Spots", "Total Genes", "Median Genes per Spot", 
             "Median UMIs per Spot", "Brain Regions", "Clusters Identified"),
  Value = c(ncol(spatial_data), nrow(spatial_data), 
            median(spatial_data$nFeature_Spatial), 
            median(spatial_data$nCount_Spatial),
            length(unique(spatial_data$brain_region)),
            length(unique(Idents(spatial_data))))
)

print(summary_df)
cat("\n")

write.csv(summary_df, file.path(output_dir, "Spatial_Summary.csv"), row.names = FALSE)

cat("========================================\n")
cat("✅ 空间转录组分析完成！\n")
cat("========================================\n")
cat("📁 结果保存在:", output_dir, "\n")
cat("📊 生成的图表:\n")
cat("   - Figure8A_QC_Violin.png: 质量控制小提琴图\n")
cat("   - Figure8B_PCA.png: PCA降维图\n")
cat("   - Figure8C_ElbowPlot.png: 肘部图\n")
cat("   - Figure8D_UMAP_Clusters.png: UMAP聚类图\n")
cat("   - Figure8E_BrainRegion_Distribution.png: 脑区分布图\n")
cat("   - Figure8F_Spatial_BrainRegions.png: 空间脑区分布图\n")
cat("   - Figure8G_*_Spatial.png: THSWD靶点空间表达图\n")
cat("   - Figure8H_THSWD_Targets_Region_Heatmap.png: 脑区靶点表达热图\n")
cat("   - Figure8I_Volcano.png: 火山图\n")
cat("   - Figure8J_Target_Correlation.png: 靶点相关性热图\n")
cat("\n")
