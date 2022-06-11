pkgs <- c(
    "Seurat", "SeuratWrappers", "ggplot2", "batchelor",
    "dplyr", "optparse", "reshape2", "data.table", "magrittr",
    "patchwork", "scales", "GSVA", "RColorBrewer", "ggridges",
    "clusterProfiler", "survminer", "survminer", "monocle",
    "psych", "ggrepel", "pheatmap", "escape", "multcomp", "agricolae"
)

lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) # nolint

# lymphocyte--9
yr <- c(
    "#fafac0", "#f5eca6", "#fee391", "#fec44f",
    "#fe9929", "#ec7014", "#cc4c02", "#8c2d04",
    "#611f03"
)
# endothelium--5
gr <- c(
    "#b9e9e0", "#7dc9b7", "#59b898", "#41ae76",
    "#16d355", "#238b45", "#116d37", "#025826",
    "#003516"
)
# fibroblast--9
bl <- c(
    "#82cbf5", "#7ba7e0", "#5199cc", "#488dbe",
    "#3690c0", "#0570b0", "#0d71aa", "#045a8d",
    "#023858"
)
# myeloid--7
pur <- c(
    "#8c97c6", "#a28abd", "#997abd", "#9362cc",
    "#88419d", "#810f7c", "#4d004b"
)
# plasma--6
bro <- c(
    "#8c510a", "#995401", "#be7816", "#be9430",
    "#ad8d36", "#a07540"
)

color1 <- c(bl[1:8], yr, pur, gr[c(3, 5, 6, 7, 9)], bro, "#A9A9A9")
color2 <- c(
    "#E0D39B", "#D05146", "#748EAE", "#567161",
    "#574F84", "#967447"
)

data <- readRDS("final_input1.Rds")


## Tissue and sample number analysis
tissue_sample_number <- read.table("number.txt", header = T, sep = "\t")
tissue_sample_number <- melt(tissue_sample_number)
tissue_sample_number$variable <- factor(
    tissue_sample_number$variable,
    levels = c("Normal", "Adjacent", "Tumor_1", "Tumor_2")
)

pdf("F2B.pdf", width = 6, height = 3)
ggplot(tissue_sample_number, aes(x = tissue, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#9362cc", "#5199cc", "#fe9929", "#fcbf47")) +
    theme_bw() +
    labs(x = "", y = "Number of samples") +
    theme(
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10), # nolint
        axis.text.y = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)
    )
dev.off()


## The Umap of clusters and celltype
mat <- data.frame(data@reductions$umap@cell.embeddings, group = data$celltype)
mt <- mat[order(mat$group != "Epithelium"), ]


pdf("umap_subtype.pdf", width = 10, height = 9)
ggplot(mt, aes(x = UMAP_1, y = UMAP_2, color = group)) +
    geom_point(size = 1e-5) +
    theme_classic() +
    scale_color_manual(values = color2) +
    theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size = 8)))
dev.off()

pdf("umap_subcluster.pdf", width = 10, height = 10)
DimPlot(data, label = F) +
    NoLegend() +
    scale_color_manual(values = color1) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )
dev.off()


## Proportion chart of cells from different sources in each cluster
percent_result <- read.table("percent_result.txt", header = T, sep = "\t")
percent_result <- melt(percent_result)
percent_result$variable <- factor(
    percent_result$variable,
    levels = c("Normal", "Adjacent", "Tumor")
)

pdf("Percent_types.pdf", width = 4, height = 3)
ggplot(percent_result, aes(x = tissue, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(y = "Proportion (%)") +
    scale_fill_manual(values = c("#9362cc", "#5199cc", "#fe9929")) +
    coord_flip() +
    theme(
        axis.line = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.6, "cm"),
        axis.text.x = element_text(
            angle = 45, hjust = 1,
            color = "black", size = 7
        ),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = percent)
dev.off()

percent_subtypes <- read.table("percent_subtype.txt", header = T, sep = "\t")
percent_subtypes <- melt(percent_subtypes)
percent_subtypes$variable <- factor(
    percent_subtypes$variable,
    levels = c("Normal", "Adjacent", "Tumor")
)
percent_subtypes$cluster <- factor(
    percent_subtypes$cluster,
    levels = c(paste0("c", 34:1), "total")
)

plot_list <- foreach::foreach(gp = unique(percent_subtypes$group)) %do% {
    plot_dat <- percent_subtypes[percent_subtypes$group == gp, ]
    ggplot(plot_dat, aes(x = cluster, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = "fill") +
        labs(y = "Proportion (%)", title = i) +
        scale_fill_manual(values = c("#9362cc", "#5199cc", "#fe9929")) +
        coord_flip() +
        theme(
            axis.line = element_blank(),
            legend.title = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = 12),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_blank(),
            legend.key.size = unit(0.6, "cm"),
            axis.text.x = element_text(
                angle = 45, hjust = 1,
                color = "black", size = 7
            ),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_continuous(expand = c(0, 0.01), labels = percent)
}
final <- wrap_plots(plot_list, ncol = length(plot_list), guides = "collect")
ggsave("percent_subtype.pdf", final, width = 10, height = 3)


## FeaturePlot or each tissue
subtype <- c(
    "Bladder", "Breast", "Colorectal", "Gastric", "Intrahepatic duct",
    "Lung", "Ovarian", "Pancreas", "Prostate", "Thyroid"
)

plot_tissue <- foreach::foreach(ts = unique(data$tissue)) %do% {
    data$co <- ifelse(
        data$tissue != ts | data$clusters == 34, "Others", subtype[k]
    )
    mat <- data.frame(data@reductions$umap@cell.embeddings, group = data$co)
    mt <- rbind(mat[which(mat$group == "Others"), ], mat[which(mat$group == subtype[k]), ]) # nolint
    mt$group <- factor(mt$group, levels = c("Others", subtype[k]))

    ggplot(mt, aes(x = cluster, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = "fill") +
        labs(y = "Proportion (%)", title = i) +
        scale_fill_manual(values = c("#9362cc", "#5199cc", "#fe9929")) +
        coord_flip() +
        theme(
            axis.line = element_blank(),
            legend.title = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = 12),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_blank(),
            legend.key.size = unit(0.6, "cm"),
            axis.text.x = element_text(
                angle = 45, hjust = 1,
                color = "black", size = 7
            ),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_continuous(expand = c(0, 0.01), labels = percent)
}

final_plot_tissue <- wrap_plots(plot_tissue, ncol = length(plot_tissue) / 2, guides = "collect") # nolint
ggsave("featureplot_each_cluster.pdf", final_plot_tissue, width = 30, height = 13) # nolint


## Subcluster analysis for DC
dc_clusters <- subset(data, ident = 13)
dc_clusters %<>% NormalizeData(object = ., normalization.method = "LogNormalize") %>% # nolint
    FindVariableFeatures(selection.method = "vst")

all.genes <- rownames(dc_clusters)
dc_clusters %<>% ScaleData(object = ., features = all.genes) %>%
    ScaleData(object = ., vars.to.regress = "percent.mt")
dc_clusters %<>% RunPCA(object = ., features = VariableFeatures(object = .)) %>%
    RunFastMNN(object = ., object.list = SplitObject(., split.by = "SampleID"))

dc_clusters %<>% FindNeighbors(reduction = "mnn", dims = 1:30) %>%
    FindClusters(resolution = 0.7) %<>%
    RunTSNE(reduction = "mnn", dims = 1:30) %>%
    RunUMAP(reduction = "mnn", dims = 1:30)

dc_clusters <- RenameIdents(dc_clusters,
    "0" = "1", "1" = "1", "2" = "1", "3" = "1",
    "8" = "1", "7" = "5", "5" = "2", "4" = "4",
    "6" = "3", "10" = "3", "9" = "6"
)

dc_clusters$seurat_clusters <- dc_clusters@active.ident
dc_clusters$seurat_clusters <- factor(dc_clusters$seurat_clusters, levels = 1:6)
saveRDS(dc_clusters, "DC_sub-clusters.Rds")

pdf("DC_sub-clusters_VlnPlot.pdf", width = 3, height = 6)
VlnPlot(dc_clusters,
    features = c("CD3D", "CD3E"),
    pt.size = 0, cols = pal_npg()(6), ncol = 1
) &
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.title.x = element_blank()
    )
dev.off()

dot_plot_genes <- c(
    "CLEC10A", "FCGR2B", "FCER1G", "FCGR2A", "S100B", "LTB",
    "CD1A", "CD1E", "STMN1", "TUBB", "TYMS", "TOP2A", "CLEC9A",
    "LGALS2", "CPVL", "XCR1", "CCL22", "BIRC3", "CCL19", "CCR7",
    "TRAC", "CD3D", "CD3E", "CD2"
)

dc_clusters@active.ident <- factor(dc_clusters@active.ident, levels = 6:1)
pdf("DC_sub-clusters_DotPlot.pdf", width = 8, height = 3.5)
DotPlot(dc_clusters,
    features = dot_plot_genes,
    cols = c(
        "#FAFCFA", "#4D9157"
    )
) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank()
    )
dev.off()

## Subcluster analysis for Endothelium
endothelium_clusters <- subset(data, ident = c(25, 26, 28))
endothelium_clusters$div <- ifelse(endothelium_clusters@active.ident == 25, "TEC", "NEC") # nolint
endothelium_clusters@active.ident <- factor(endothelium_clusters$div, levels = c("NEC", "TEC")) # nolint

pdf("endothelium_cluster_vln_plot.pdf", width = 4.5, height = 7)
VlnPlot(endothelium_clusters, features = c(
    "IGFBP2", "IGFBP4", "INSR",
    "SPRY1", "CD320", "IGHG4"
), pt.size = 0, cols = pal_npg()(6), ncol = 2) &
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.title.x = element_blank()
    )
dev.off()

TEC_NEC_DEG_Marker <- FindAllMarkers(endothelium_clusters, min.pct = 0.25, only.pos = T) # nolint

TEC_NEC_DEG_Marker$fc <- ifelse(
    TEC_NEC_DEG_Marker$cluster == "NEC",
    -1 * TEC_NEC_DEG_Marker$avg_log2FC, TEC_NEC_DEG_Marker$avg_log2FC
) # nolint
TEC_NEC_DEG_Marker$q_value <- -log(TEC_NEC_DEG_Marker$p_val_adj, 10)
TEC_NEC_DEG_Marker$group <- ifelse(
    TEC_NEC_DEG_Marker$fc > 0.5 & TEC_NEC_DEG_Marker$p_val_adj < 0.05,
    "sig up",
    ifelse(TEC_NEC_DEG_Marker$fc < -0.5 & TEC_NEC_DEG_Marker$p_val_adj < 0.05,
        "sig down", "not sig"
    )
)
select_genes <- c(
    TEC_NEC_DEG_Marker[which(TEC_NEC_DEG_Marker$group == "sig up"), "gene"][1:10], # nolint
    TEC_NEC_DEG_Marker[which(TEC_NEC_DEG_Marker$group == "sig down"), "gene"][1:10] # nolint
)
TEC_NEC_DEG_Marker$label <- ifelse(
    TEC_NEC_DEG_Marker$gene %in% select_genes,
    TEC_NEC_DEG_Marker$gene, NA
)

p <- ggplot(TEC_NEC_DEG_Marker, aes(x = fc, y = q_value)) +
    geom_point(aes(color = group), size = 1) +
    scale_color_manual(values = c("gray", "blue", "red")) +
    labs(x = "log2 fold change", y = "-log10 padj", title = "CAMR vs stable") +
    geom_text_repel(aes(label = label), size = 2) +
    geom_hline(yintercept = 1.30103, linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("TEC_NEC_DEG_Volcano_Plot.pdf", p, width = 6.5, height = 6)


## Cluster analysis of tissue origin specificity
Percent_Calculated <- read.table(
    "Percent_Calculated.txt",
    header = T, sep = "\t"
)
Percent_Calculated <- melt(Percent_Calculated)
Percent_Calculated$variable <- factor(
    Percent_Calculated$variable,
    levels = c("c21", "c17", "c6", "c5", "c4")
)
npg_color <- pal_npg(alpha = 0.8)(10)
color_for_use <- rev(
    c(
        "#9269A2", "#7dc9b7", colo[3], "#5199cc", "#6479cc",
        npg_color[c(10, 1, 8, 2, 6, 9)],
        "#EEC877", npg_color[5]
    )
)

p <- ggplot(Percent_Calculated, aes(x = variable, y = value, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(y = "Proportion (%)") +
    scale_fill_manual(values = color_for_use) + # nolint
    coord_flip() +
    theme(
        axis.line = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.6, "cm"),
        axis.text.x = element_text(
            angle = 45, hjust = 1,
            color = "black", size = 10
        ),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = percent)
ggsave(
    "Tissue_Origin_Specificity_Percent_Calculated.pdf", p,
    width = 5, height = 3.8
)


## Analysis of the distribution of cells of tumor, adjacent and normal tissue origin, respectively
plot_list <- foreach::foreach(gp = c("Normal", "Adjacent", "Tumor")) %do% {
    select_subsets <- data
    select_subsets$group_for_color <- ifelse(
        select_subsets$clusters != 35 & select_subsets$group == gp, gp, "Others"
    )
    plot_data <- data.frame(
        select_subsets@reductions$umap@cell.embeddings,
        group = select_subsets$group_for_color
    )
    plot_data$group <- factor(plot_data$group, levels = c("Others", gp))
    plot_data <- plot_data[order(plot_data$group), ]

    ggplot(tt, aes(x = UMAP_1, y = UMAP_2, color = group)) +
        geom_point(size = 0.1) +
        scale_color_manual(values = c("#D9D9D9", "#9362cc")) +
        labs(title = gp) +
        theme_classic() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 25, color = "black"),
            legend.position = "none"
        )
}
final <- wrap_plots(plot_list, ncol = length(plot_list), guides = "collect")
ggsave("Subtype_Umap_Plot", final, width = 24, height = 8)


## Subcluster analysis for C24
C24_cluster <- subset(data, ident = 24)
C24_cluster_clusters %<>% NormalizeData(object = ., normalization.method = "LogNormalize") %>% # nolint
    FindVariableFeatures(selection.method = "vst")

all.genes <- rownames(C24_cluster)
C24_cluster %<>% ScaleData(object = ., features = all.genes) %>%
    ScaleData(object = ., vars.to.regress = "percent.mt")
C24_cluster %<>% RunPCA(object = ., features = VariableFeatures(object = .)) %>%
    RunFastMNN(object = ., object.list = SplitObject(., split.by = "SampleID"))

C24_cluster %<>% FindNeighbors(reduction = "mnn", dims = 1:30) %>%
    FindClusters(resolution = 0.7) %<>%
    RunTSNE(reduction = "mnn", dims = 1:30) %>%
    RunUMAP(reduction = "mnn", dims = 1:30)

pdf("C24_Sub_Cluster_Umap.pdf", width = 4, height = 3.5)
DimPlot(C24_cluster, label = F) +
    scale_color_manual(values = pal_npg()(6))
dev.off()

pdf("C24_Sub_Cluster_Vln_Plot.pdf", width = 3, height = 4.5)
VlnPlot(C24_cluster,
    features = "RGS13",
    pt.size = 0, cols = pal_npg()(6), ncol = 1
) &
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.title.x = element_blank()
    )
dev.off()

## Analysis of FABP4+ macrophage (c20) proportion
Validation_Results <- readRDS("Validation_Results/input.Rds") # nolint
Validation_Matrix <- table(Validation_Results@meta.data[, c("sampleid", "seurat_clusters")]) # nolint
Validation_Matrix_Sum <- 100 * rowSums(Validation_Matrix[, c(4, 32)]) / as.vector(table(Validation_Results$sampleid)) # nolint
Validation_Matrix_Value <- data.frame(name = names(Validation_Matrix_Sum), value = as.numeric(Validation_Matrix_Sum)) # nolint
Reference <- unique(Validation_Results@meta.data[, c("sampleid", "group")]) # nolint
Validation_Matrix_Draw <- merge(Validation_Matrix_Value, Reference, by.x = "name", by.y = "sampleid") # nolint

Validation_Matrix_Draw <- Validation_Matrix_Draw[
    which(Validation_Matrix_Draw$group %in% c("Lung_N", "Lung_T", "tLB", "mBrain")), # nolint
]
Validation_Matrix_Draw$group <- factor(
    Validation_Matrix_Draw$group,
    levels = c("Lung_N", "Lung_T", "tLB", "mBrain")
)

model <- aov(value ~ group, data = tt)
rht <- glht(model, linfct = mcp(group = "Dunnett"), alternative = "two.side")
summary(rht)

p <- ggplot(tt, aes(x = group, y = value)) +
    geom_boxplot(
        size = 0.8, fill = "white", outlier.fill = NA,
        outlier.color = NA, outlier.size = 0
    ) +
    geom_point(aes(fill = group, color = group), shape = 21, size = 3) +
    scale_color_manual(values = pal_npg()(4)) +
    scale_fill_manual(values = pal_npg()(4)) +
    theme_classic() +
    labs(y = "Proportion (%)") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(
            size = 10, color = "black", angle = 45, hjust = 1
        ),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
ggsave("Validation_FABP4_Macrophage_Proportion.pdf", width = 2.3, height = 3)
