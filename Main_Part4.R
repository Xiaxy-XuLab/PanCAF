pkgs <- c(
    "Seurat", "SeuratWrappers", "ggplot2", "batchelor", "circlize",
    "dplyr", "optparse", "reshape2", "data.table", "magrittr",
    "patchwork", "scales", "GSVA", "RColorBrewer", "ggridges", "ggridges",
    "clusterProfiler", "survminer", "survminer", "monocle", "tidyverse",
    "psych", "ggrepel", "pheatmap", "escape", "multcomp", "agricolae"
)

color_for_use <- c(
    "#8F2A47", "#B83A4D", "#D25C52", "#E27E56",
    "#ECA86B", "#F4CB85", "#F8E8A2", "#FAF8C7", "#EBF0AF",
    "#CEE2A2", "#ABD3A6", "#82C3A5", "#609EB0", "#4C78B1",
    "#5C519B"
)

lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) # nolint
source("Dependent_SCENIC.R")
source("Dependent_Monocle.R")
source("Dependent_Escape.R")

data <- readRDS("final_input1.Rds")

FIB <- subset(data, ident = c(1, 2, 4, 5, 3, 6, 7, 8))
pdf("FIB_Markers1.pdf", width = 8, height = 3)
DotPlot(FIB,
    features = c(
        "ACTA2", "DCN", "CD86", "FCGR1A",
        "MRC1", "CD163", "CD68", "CD14"
    ), cols = c("#FAFCFA", "#4D9157")
) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank()
    )

dev.off()

pdf("FIB_Markers_Vln1.pdf", width = 8, height = 3)
VlnPlot(FIB,
    features = c("ACTA2", "HLA-DRA", "CD74"),
    pt.size = 0, ncol = 3
) &
    scale_fill_manual(values = color1) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
dev.off()

pdf("FIB_Markers_Vln1.pdf", width = 8, height = 3)
VlnPlot(FIB,
    features = c("ACTA2", "HLA-DRA", "CD74"),
    pt.size = 0, ncol = 3
) &
    scale_fill_manual(values = pal_npg()(8)) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
dev.off()

pdf("FIB_Markers_Vln2.pdf", width = 8, height = 3)
VlnPlot(FIB,
    features = c("MPZ", "S100B", "PLP1", "LGI4"),
    pt.size = 0, ncol = 2
) &
    scale_fill_manual(values = pal_npg()(8)) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
dev.off()

## Trajactory analysis of CAFmyo, TAM and CAFapi
subset1 <- subset(data, ident = 8)

Select_Cells <- c()
for (i in c(20, 1)) {
    object_select <- subset(data, ident = i)
    set.seed(1)
    Select_Cells <- append(
        Select_Cells, rownames(object_select@meta.data[
            sample(1:nrow(object_select@meta.data), 500),
        ])
    )
}

monotj <- subset(
    data,
    cells = c(Select_Cells, rownames(object_select@meta.data))
)
monocds <- RunMonocle(monotj)
monocds <- orderCells(monocds, root_state = 5)

pdf("CAF_TAM_Trajactory_State.pdf", width = 3.5, height = 4)
plot_cell_trajectory(
    monocds,
    cell_size = 1, color_by = "State", show_tree = FALSE,
    show_branch_points = FALSE
) +
    scale_color_manual(values = pal_aaas(alpha = 0.8)(9)[c(1, 2, 4, 3, 9)]) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
    )
dev.off()

monocds@phenoData@data$clusters <- paste0("c", monocds@phenoData@data$Clusters)
pdf("CAF_TAM_Trajactory_Cluster.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
    monocds,
    cell_size = 1, color_by = "clusters", show_tree = FALSE,
    show_branch_points = FALSE
) +
    scale_color_manual(values = pal_npg()(9)[c(4, 2, 5)]) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
    )
dev.off()

pdf("CAF_TAM_Trajactory_Pseudotime1.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
    monocds,
    cell_size = 1, color_by = "Pseudotime",
    show_tree = FALSE, show_branch_points = FALSE
) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
    )
dev.off()

my_pseudotime_de <- differentialGeneTest(
    monocds,
    fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 30
)
sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 10^-20))

pdf("CAF_TAM_Trajactory_Heatmap.pdf", width = 6, height = 8)
plot_pseudotime_heatmap(
    monocds[sig_gene_names, ],
    num_clusters = 4, cores = 20,
    use_gene_short_name = TRUE, show_rownames = TRUE
)
dev.off()

Gene_label <- plot_pseudotime_heatmap(
    monocds[sig_gene_names, ],
    num_clusters = 4, cores = 60,
    use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE
)
clusters <- cutree(Gene_label$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[, 1] <- as.character(clustering[, 1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.table(clustering, "CAF_TAM_Trajactory_label.txt", quote = F, sep = "\t")
Gene_label <- data.frame(cluster = clustering, gene = rownames(clustering))

plotdf <- pData(monocds)
plotdf$clusters <- factor(plotdf$clusters, levels = c("c8", "c20", "c1"))

p <- ggplot(plotdf, aes(x = Pseudotime, y = clusters, fill = clusters)) +
    geom_density_ridges(scale = 1) +
    scale_fill_manual(values = pal_npg()(9)[c(5, 2, 4)]) +
    geom_vline(xintercept = c(0, 5, 10, 15), linetype = 2) +
    scale_y_discrete("") +
    theme_minimal() +
    theme(panel.grid = element_blank())
ggsave("CAF_TAM_Trajactory_Density.pdf", width = 4, height = 1.5)

## SCENIC analysis
Select_Counts <- as.matrix(monotj@assays$RNA@counts)
RunSCENIC(Select_Counts)

scenic <- read.table("moreSCENIC.txt", header = T, sep = "\t")
name <- gsub("-", "_", rownames(monotj@meta.data))
colnames(scenic) <- name

mat <- data.frame(
    cell = name,
    celltype = paste0("c", monotj@meta.data[, "Clusters"]),
    state = paste0("s", monocds@phenoData@data$State)
)

mat$celltype <- factor(mat$celltype, levels = c("c1", "c20", "c8"))
mat$state <- factor(mat$state, levels = paste0("s", c(1, 5, 2, 3, 4)))
mat <- mat[order(mat$celltype, mat$state), ]
rownames(mat) <- mat$cell
mat <- as.data.frame(mat[, c(2, 3)])

col1 <- pal_npg()(9)[c(2, 5, 4)]
names(col1) <- levels(mat$celltype)
col2 <- pal_aaas(alpha = 0.8)(9)[c(1, 9, 2, 4, 3)]
names(col2) <- paste0("s", c(1, 5, 2, 3, 4))

ann_colors <- list(
    celltype = col1,
    state = col2
)
scen <- scenic[
    c(
        "MYLK (60g)",
        "MAFB_extended (258g)", "SPI1 (1010g)"
    ),
]
bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("CAF_TAM_SCENIC.pdf", width = 8, height = 6)
pheatmap(scen[, rownames(mat)],
    show_colnames = FALSE, annotation_col = mat,
    annotation_colors = ann_colors,
    show_rownames = T, scale = "row", cluster_cols = FALSE,
    cluster_rows = T,
    color = c(
        colorRampPalette(
            colors = brewer.pal(11, "RdYlBu")[11:6]
        )(length(bk) / 2),
        colorRampPalette(
            colors = brewer.pal(11, "RdYlBu")[6:1]
        )(length(bk) / 2)
    ),
    breaks = bk, annotation_legend = TRUE,
    legend_breaks = c(-1, 1), legend_labels = c("Low", "High"),
    fontsize = 12, annotation_names_col = TRUE, border_color = "gray"
)
dev.off()

## Similarity analysis
col1 <- c(
    "#8F2A47", "#B83A4D", "#D25C52", "#E27E56",
    "#ECA86B", "#F4CB85", "#F8E8A2", "#FAF8C7", "#EBF0AF",
    "#CEE2A2", "#ABD3A6", "#82C3A5", "#609EB0", "#4C78B1",
    "#5C519B"
)
simi <- function(x) {
    Object_Select <- subset(data, ident = x)
    Object_Select <- FindVariableFeatures(
        Object_Select,
        selection.method = "vst", nfeatures = 5000
    )
    Expr_Matrix <- as.data.frame(t(
        as.data.frame(
            Object_Select@assays$RNA@data[VariableFeatures(Object_Select), ]
        )
    ))
    Expr_Matrix$cluster <- as.character(Object_Select@active.ident)
    Expr_Matrix_Mean <- aggregate(
        Expr_Matrix[, 1:(ncol(Expr_Matrix) - 1)],
        by = list(Expr_Matrix$Cluster), FUN = mean
    )
    rownames(Expr_Matrix_Mean) <- paste0("c", Expr_Matrix_Mean$Group.1)
    Expr_Matrix_Mean <- Expr_Matrix_Mean[, 2:ncol(Expr_Matrix_Mean)]
    Expr_Matrix_Mean <- as.data.frame(t(Expr_Matrix_Mean))
    Corr_Result <- cor(Expr_Matrix_Mean)
    pdf("CAF_TAM_Cor.pdf", width = 7, height = 7)
    corrplot(Corr_Result,
        method = "color",
        col = rev(brewer.pal(11, "RdYlBu")),
        type = "upper",
        tl.col = "black",
        order = "hclust", is.corr = FALSE, addCoef.col = "grey"
    )
    dev.off()
}
x <- c(1:9, 19:21, 23)
simi(x)
x <- c(1:8, 25:28)
simi(x)

## Proportion analysis
CAF_Select <- subset(data, ident = 6)
Tissue_Origin <- melt(
    table(
        CAF_Select@meta.data[, c("tissue", "group")]
    ) / as.vector(table(CAF_Select$tissue))
)
Tissue_Origin$group <- factor(
    Tissue_Origin$group,
    levels = c("Normal", "Adjacent", "Tumor")
)

p <- ggplot(mat, aes(x = name, y = percent, fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = pal_aaas()(3)[c(1, 3, 2)]) +
    theme_bw() +
    labs(x = "", y = "Proportion (%)") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()
    )
ggsave("Tissue_Origin_Proportion.pdf", p, width = 5, height = 3)

## VlnPlot for ESM1
subsets <- subset(data, ident = c(1:8, 26:29))
VlnPlot(subsets, features = "ESM1", pt.size = 0) +
    labs(x = "Cluster", y = "Expression level")
