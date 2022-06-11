pkgs <- c(
    "Seurat", "SeuratWrappers", "ggplot2", "batchelor", "circlize",
    "dplyr", "optparse", "reshape2", "data.table", "magrittr",
    "patchwork", "scales", "GSVA", "RColorBrewer", "ggridges",
    "clusterProfiler", "survminer", "survminer", "monocle",
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

## The interaction of TME in normal, adjacent and tumor.
Count_Network <- read.table("/work/xiaxy/work/Pancancer/GO/go4/r15/out/count_network.txt", header = T, sep = "\t") # nolint

screen <- function(x) {
    Count_Network_Screen <- Count_Network[
        intersect(
            grep(x, Count_Network$SOURCE), grep(x, Count_Network$TARGET)
        ),
    ]
    Count_Network_Screen <- dcast(Count_Network_Screen, SOURCE ~ TARGET)
    rownames(Count_Network_Screen) <- Count_Network_Screen$SOURCE
    Count_Network_Screen <- Count_Network_Screen[, -1]
    return(Count_Network_Screen)
}

bk <- c(seq(0, 5, by = 0.01))
pdf("TME_Interaction_Heatmap_Tumor.pdf", width = 5.5, height = 5)
p <- pheatmap(log1p(screen("Tumor")),
    show_colnames = TRUE,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = TRUE, scale = "none",
    cluster_rows = T, color = c(
        colorRampPalette(colors = col1[15:11])(length(bk) * 2 / 3),
        colorRampPalette(colors = col1[11:1])(length(bk) / 3)
    ),
    breaks = bk, annotation_legend = TRUE,
    treeheight_row = 0, treeheight_col = 0,
    legend_breaks = c(0, 5), legend_labels = c("Low", "High"),
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

Adjacent_Order <- c(
    "Fibroblast_Adjacent", "Myeloid_Adjacent", "Endothelium_Adjacent",
    "Epithelium_Adjacent", "Lymphocyte_Adjacent", "Plasma_Adjacent"
)
pdf("TME_Interaction_Heatmap_Adjacent.pdf", width = 5.5, height = 5)
p <- pheatmap(log1p(screen("Adjacent")[Adjacent_Order, Adjacent_Order]),
    show_colnames = TRUE,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE, scale = "none",
    cluster_rows = FALSE, color = c(
        colorRampPalette(colors = col1[15:11])(length(bk) * 2 / 3),
        colorRampPalette(colors = col1[11:1])(length(bk) / 3)
    ),
    breaks = bk, annotation_legend = TRUE,
    treeheight_row = 0, treeheight_col = 0,
    legend_breaks = c(0, 5), legend_labels = c("Low", "High"),
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

Normal_Order <- c(
    "Fibroblast_Normal", "Myeloid_Normal", "Endothelium_Normal",
    "Epithelium_Normal", "Lymphocyte_Normal", "Plasma_Normal"
)
pdf("TME_Interaction_Heatmap_Normal.pdf", width = 5.5, height = 5)
p <- pheatmap(log1p(screen("Normal")[Normal_Order, Normal_Order]),
    show_colnames = TRUE,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE, scale = "none",
    cluster_rows = FALSE, color = c(
        colorRampPalette(colors = col1[15:11])(length(bk) * 2 / 3),
        colorRampPalette(colors = col1[11:1])(length(bk) / 3)
    ),
    breaks = bk, annotation_legend = TRUE,
    treeheight_row = 0, treeheight_col = 0,
    legend_breaks = c(0, 5), legend_labels = c("Low", "High"),
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

## Comparison Dendrograms for sub_clusters of fibroblast
FIBROBLAST <- subset(data, ident = 1:34)
FIBROBLAST$Label <- paste("c", FIBROBLAST$clusters, " ", FIBROBLAST$tissue, sep = "") # nolint
FIBROBLAST <- FindVariableFeatures(FIBROBLAST,
    selection.method = "vst",
    nfeatures = 2000
)
FIBROBLAST_EXP <- FIBROBLAST@assays$RNA@data[VariableFeatures(FIBROBLAST), ]
FIBROBLAST_EXP <- t(as.matrix(FIBROBLAST_EXP))
FIBROBLAST_EXP_Mean <- aggregate(FIBROBLAST_EXP, list(FIBROBLAST$Label), mean)
rownames(FIBROBLAST_EXP_Mean) <- FIBROBLAST_EXP_Mean$Group.1
FIBROBLAST_EXP_Mean <- FIBROBLAST_EXP_Mean[, -1]

FIBROBLAST_EXP_R <- corr.test(t(FIBROBLAST_EXP_Mean), method = "pearson")
Distance <- hclust(as.dist((1 - FIBROBLAST_EXP_R$r) / 2), method = "ward.D2")
pdf("Comparison_Dendrograms_Fibroblast.pdf", width = 9, height = 9)
fviz_dend(Distance,
    k = 7,
    cex = 0.8,
    k_colors = pal_npg()(7),
    color_labels_by_k = FALSE,
    rect_border = pal_npg()(5),
    rect = TRUE,
    type = "circular",
    rect_fill = TRUE,
    horiz = TRUE,
    main = "", ylab = ""
)
dev.off()

## circlize analysis
for (i in c("Normal", "Adjacent", "Normal")) {
    Select_Object <- subset(data, cells = rownames(data@meta.data[which(data$group == i), ])) # nolint
    Select_Cells <- c()
    for (i in c(1:34)) {
        Select_Subset <- subset(Select_Object, ident = i)
        set.seed(1)
        if (nrow(Select_Subset@meta.data) > 2000) {
            Select_Cells <- append(
                Select_Cells,
                rownames(Select_Subset@meta.data[sample(1:nrow(Select_Subset@meta.data), 2000), ]) # nolint
            )
        } else {
            Select_Cells <- append(Select_Cells, rownames(Select_Subset@meta.data)) # nolint
        }
    }

    Object_Final <- subset(data, cells = Select_Cells)
    write.table(
        Select_Object@assays$RNA@counts,
        paste0(i, "/count.txt"),
        quote = F, sep = "\t"
    )
    meta <- data.frame(
        Cell = rownames(Select_Object@meta.data),
        Cell_type = paste0(Select_Object$celltype)
    )

    write.table(
        meta,
        paste0(i, "/meta.txt"),
        quote = F, sep = "\t", row.names = F
    )
}
Tumor <- read.table("Tumor/out/count_network.txt", header = T, sep = "\t") # nolint
Normal <- read.table("Normal/out/count_network.txt", header = T, sep = "\t") # nolint
Adjacent <- read.table("Adjacent/out/count_network.txt", header = T, sep = "\t") # nolint

gr <- function(x, y) {
    x <- x[which(x$SOURCE != x$TARGET), ]
    x <- dcast(x, SOURCE ~ TARGET)
    x$SOURCE <- paste(x$SOURCE, y, sep = "_")
    rownames(x) <- x$SOURCE
    x <- x[, -1]
    colnames(x) <- paste(colnames(x), y, sep = "_")
    return(x)
}
Interaction_Counts_1 <- as.matrix(gr(Tumor, "tumor"))
Interaction_Counts_2 <- as.matrix(gr(Normal, "normal"))
Interaction_Counts_3 <- as.matrix(gr(Adjacent, "adjacent"))
for (i in 1:5) {
    Interaction_Counts_1[i, i] <- 0
}
for (i in 1:5) {
    Interaction_Counts_2[i, i] <- 0
}
for (i in 1:5) {
    Interaction_Counts_3[i, i] <- 0
}
Interaction_Counts <- matrix(0, nrow = 15, ncol = 15)
rownames(Interaction_Counts) <- c(rownames(Interaction_Counts_1), rownames(Interaction_Counts_2), rownames(Interaction_Counts_3)) # nolint
colnames(Interaction_Counts) <- c(colnames(Interaction_Counts_1), colnames(Interaction_Counts_2), colnames(Interaction_Counts_3)) # nolint
Interaction_Counts[rownames(Interaction_Counts_1), colnames(Interaction_Counts_1)] <- Interaction_Counts_1 # nolint
Interaction_Counts[rownames(Interaction_Counts_2), colnames(Interaction_Counts_2)] <- Interaction_Counts_2 # nolint
Interaction_Counts[rownames(Interaction_Counts_3), colnames(Interaction_Counts_3)] <- Interaction_Counts_3 # nolint
Name_Order <- sort(rownames(Interaction_Counts))[c(2, 1, 3, 5, 4, 6, 8, 7, 9, 11, 10, 12, 14, 13, 15)] # nolint
Interaction_Counts <- Interaction_Counts[Name_Order, Name_Order]
Label_Name <- unique(unlist(dimnames(Interaction_Counts)))
group <- structure(rep(1:5, each = 3), names = Label_Name)

grid.col <- structure(
    rep(pal_aaas()(3)[c(1, 3, 2)], time = 5),
    names = names(group)
)

pdf("cir.pdf", width = 4, height = 4)
# 最外层添加一个空白轨迹
chordDiagram(
    Interaction_Counts,
    group = group, grid.col = grid.col,
    annotationTrack = c("grid", "axis"),
    preAllocateTracks = list(
        track.height = mm_h(4),
        track.margin = c(mm_h(4), 0)
    )
)
circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
        sector.index <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(
            mean(xlim), mean(ylim),
            sector.index,
            cex = 0.6,
            niceFacing = TRUE
        )
    },
    bg.border = NA
)
highlight.sector(
    names(group)[1:3],
    track.index = 1, col = "#fb8072",
    text = "Endothelium", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    names(group)[4:6],
    track.index = 1, col = "#80b1d3",
    text = "Fibroblast", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    names(group)[7:9],
    track.index = 1, col = "#fdb462",
    text = "Lymphocyte", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    names(group)[10:12],
    track.index = 1, col = "#73428A",
    text = "Myeloid", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    names(group)[13:15],
    track.index = 1, col = "#6EA056",
    text = "Plasma", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
dev.off()

## SCENIC analysis for sub_clusters of fibroblast
Select_Cells <- c()
for (i in 1:8) {
    select_cluster <- subset(data, ident = i)
    set.seed(1)
    if (nrow(select_cluster@meta.data) > 500) {
        Select_Cells <- append(
            Select_Cells,
            rownames(
                select_cluster@meta.data[sample(1:nrow(select_cluster@meta.data), 500), ] # nolint
            )
        )
    } else {
        Select_Cells <- append(Select_Cells, rownames(select_cluster@meta.data))
    }
}
Select_Object <- subset(data, cells = Select_Cells)
Select_Counts <- as.matrix(Select_Object@assays$RNA@counts)
RunSCENIC(Select_Counts)

scenic <- read.table("moreSCENIC.txt", header = T, sep = "\t")
name <- gsub("-", "_", rownames(Select_Object@meta.data))
colnames(scenic) <- name

annotation <- data.frame(
    celltype = paste0("c", Select_Object$Clusters)
)
rownames(annotation) <- name
colnames(annotation) <- "celltype"

color1 <- pal_npg()(9)[c(2, 5, 6, 3, 4, 1, 7, 8)]
names(color1) <- paste0("c", 1:8)
ann_colors <- list(
    celltype = color1
)
Cellname <- c()
for (i in paste0("c", Select_Object$Clusters)) {
    Cellname <- append(Cellname, rownames(annotation)[which(annotation$celltype == i)]) # nolint
}

SCENIC_Select <- scenic[
    c(
        "EGR1 (14g)", "CREB3L1 (84g)", "TCF21 (13g)",
        "TWIST2_extended (17g)",
        "SOX2_extended (101g)", "SOX10 (123g)", "SOX17_extended (13g)",
        "ETS1 (305g)", "ETS2 (60g)", "SOX18 (59g)", "ERG (90g)",
        "MAFB_extended (51g)", "SPI1 (338g)", "IKZF1 (64g)",
        "MEF2C_extended (50g)", "TBX2_extended (24g)"
    ), Cellname
]

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))

pdf("FIBROBLAST_SCENIC_Heatmap.pdf", width = 8, height = 3.5)
pheatmap(SCENIC_Select,
    show_colnames = FALSE, annotation_col = annotation,
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

## SCENIC analysis for sub_clusters of fibroblast
Select_cells <- c()
for (i in 1:8) {
    Fibroblast_Select <- subset(data, ident = i)
    set.seed(1)
    if (nrow(Fibroblast_Select@meta.data) > 5000) {
        Select_cells <- append(
            Select_cells,
            rownames(Fibroblast_Select@meta.data[sample(1:nrow(Fibroblast_Select@meta.data), 5000), ]) # nolint
        )
    } else {
        Select_cells <- append(Select_cells, rownames(Fibroblast_Select@meta.data)) # nolint
    }
}
sce <- subset(data, cells = Select_cells)
RunEscape(sce)

esc <- read.table("Hall_Marker.txt", header = T, sep = "\t")
mat <- aggregate(esc, list(sce$Clusters), mean)
rownames(mat) <- mat$Group.1
mat <- mat[c(1, 2, 4, 6:8), -1]
colnames(mat) <- gsub("HALLMARK_", "", colnames(mat))
mat <- as.data.frame(t(mat))
annotation_col <- data.frame(
    celltype = colnames(mat)
)
rownames(annotation_col) <- colnames(mat)
col1 <- pal_npg(alpha = 0.9)(9)[c(1, 2, 4, 6, 7, 8)]
names(col1) <- colnames(mat)

ann_colors <- list(
    celltype = col1
)

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("pheatmap.pdf", width = 4.2, height = 4)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, scale = "row", cluster_cols = FALSE,
    cluster_rows = T, color = c(
        colorRampPalette(colors = brewer.pal(11, "RdYlBu")[11:6])(length(bk) / 2), # nolint
        colorRampPalette(colors = brewer.pal(11, "RdYlBu")[6:1])(length(bk) / 2) # nolint
    ),
    breaks = bk, annotation_legend = TRUE,
    legend_breaks = c(-1, 1),
    legend_labels = c("Low", "High"),
    fontsize = 6, annotation_names_col = FALSE,
    border_color = "gray", treeheight_row = 10
)
dev.off()

## Interaction comparison in tumor, adjacent and normal tissue
nom <- read.table("Normal/out/count_network.txt", header = T, sep = "\t") # nolint
adj <- read.table("Adjacent/out/count_network.txt", header = T, sep = "\t") # nolint
tum <- read.table("Tumor/out/count_network.txt", header = T, sep = "\t") # nolint

deal <- function(x) {
    x <- x[which(x$SOURCE != x$TARGET), ]
    interaction_table <- matrix(NA, nrow = length(unique(x$SOURCE)), ncol = 2) # nolint
    k <- 1
    for (i in unique(x$SOURCE)) {
        table_select <- x[which(x$SOURCE == i), ]
        interaction_table[k, 1] <- i
        interaction_table[k, 2] <- sum(table_select$count)
        k <- k + 1
    }
    colnames(interaction_table) <- c("iterm", "count")
    interaction_table <- as.data.frame(interaction_table)
    interaction_table$type <- sapply(
        strsplit(interaction_table$iterm, "_", fixed = T), "[", 1
    )
    interaction_table$tissue <- sapply(
        strsplit(interaction_table$iterm, "_", fixed = T), "[", 2
    )
    return(interaction_table)
}

npg_color <- pal_npg(alpha = 0.8)(10)
color_new <- rev(c(
    "#9269A2", colo[3], "#6479cc",
    npg_color[c(10, 8, 2, 6, 9)], "#EEC877", npg_color[5]
))

p1 <- ggplot(deal(tum), aes(x = type, y = as.numeric(count))) +
    geom_boxplot() +
    geom_jitter(
        aes(color = tissue, fill = tissue),
        shape = 21, width = 0.2, size = 3
    ) +
    scale_color_manual(values = color_new) +
    scale_fill_manual(values = color_new) +
    theme_bw() +
    labs(x = "", y = "The counts of interaction") +
    theme(
        axis.text.x = element_text(
            size = 10, angle = 45,
            hjust = 1, color = "black"
        ),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank()
    )
ggsave("Interaction_Tumor.pdf", p1, width = 4.2, height = 3.2)

p1 <- ggplot(deal(adj), aes(x = type, y = as.numeric(count))) +
    geom_boxplot() +
    geom_jitter(
        aes(color = tissue, fill = tissue),
        shape = 21, width = 0.2, size = 3
    ) +
    scale_color_manual(values = color_new[c(1, 3:8, 10)]) +
    scale_fill_manual(values = color_new[c(1, 3:8, 10)]) +
    theme_bw() +
    labs(x = "", y = "The counts of interaction") +
    theme(
        axis.text.x = element_text(
            size = 10, angle = 45,
            hjust = 1, color = "black"
        ),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank()
    )
ggsave("Interaction_Adjacent.pdf", p1, width = 4.2, height = 3.2)


p1 <- ggplot(deal(nom), aes(x = type, y = as.numeric(count))) +
    geom_boxplot() +
    geom_jitter(
        aes(color = tissue, fill = tissue),
        shape = 21, width = 0.2, size = 3
    ) +
    scale_color_manual(values = color_new[c(1, 3:4, 6, 8, 9, 10)]) +
    scale_fill_manual(values = color_new[c(1, 3:4, 6, 8, 9, 10)]) +
    theme_bw() +
    labs(x = "", y = "The counts of interaction") +
    theme(
        axis.text.x = element_text(
            size = 10, angle = 45,
            hjust = 1, color = "black"
        ),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank()
    )
ggsave("Interaction_Normal.pdf", p1, width = 4.2, height = 3.2)

## FeaturePlot for sub-clusters of fibroblast
Fibroblast_Subset <- subset(data, ident = 1:5)
Fibroblast_Subset@active.ident <- factor(
    paste0("c", Fibroblast_Subset@active.ident),
    levels = c("c3", "c5", "c4", "c2", "c1")
)
names(Fibroblast_Subset@active.ident) <- rownames(Fibroblast_Subset@meta.data)

Embeddings <- as.data.frame(Fibroblast_Subset@reductions$umap@cell.embeddings)
Select_Cells <- rownames(Embeddings[which(Embeddings$UMAP_1 < 2 & Embeddings$UMAP_2 > -2 & Embeddings$UMAP_2 < 5), ]) # nolint
Select_Object <- subset(Fibroblast_Subset, cells = Select_Cells)

pdf("newout/sf2/c_fea.pdf", width = 7, height = 6.5)
FeaturePlot(Select_Object,
    features = c("ACTA2", "CFD", "FAP", "TGFB1"),
    cols = rev(color_for_use), ncol = 2, order = T, pt.size = 0.5
)
dev.off()

## DotPlot analysis for sub-clusters of fibroblast
CAF_Clusters <- subset(data, ident = 1:8) # nolint
CAF_Clusters_DEG_Markers <- FindAllMarkers( # nolint
    CAF_Clusters,
    only.pos = T, min.pct = 0.25
)

CAF_Clusters_DEG_Genes <- do.call( # nolint
    rbind, lapply(
        split(data.table(CAF_Clusters_DEG_Markers), by = "cluster"),
        function(x) x[1:5, "gene"]
    )
)

CAF_Clusters@active.ident <- factor(CAF_Clusters@active.ident, levels = 8:1)
pdf("CAF_Clusters_DEG_DotPlot.pdf", width = 8, height = 3)
DotPlot(CAF_Clusters,
    features = unique(CAF_Clusters_DEG_Genes$gene),
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


## Marker analysis of CAF and NIF
CAF_NIF <- subset(data, ident = 1:5)
CAF_NIF$group <- ifelse(CAF_NIF@active.ident %in% c(1, 2, 4), "CAF", "NIF")
CAF_NIF$label <- ifelse(CAF_NIF@active.ident %in% c(3, 5), "NIF",
    ifelse(CAF_NIF@active.ident == 1, "CAFmyo",
        ifelse(CAF_NIF@active.ident == 2, "CAFinfla", "CAFadi")
    )
)
CAF_NIF@active.ident <- factor(sub$group)
names(CAF_NIF@active.ident) <- rownames(CAF_NIF@meta.data)

pdf("CAF_NIF_VlnPlot.pdf", width = 9, height = 4.5)
VlnPlot(sub,
    features = c("PDGFRA", "PDGFRB", "NOTCH3", "HES4", "FAP", "THY1"),
    pt.size = 0, ncol = 4
) &
    scale_fill_manual(values = pal_npg()(10)) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
dev.off()

## Trajactory analysis of CAF and NIF
CAF_NIF <- subset(data, ident = 1:5)
CAF_NIF@active.ident <- factor(
    paste0("c", CAF_NIF@active.ident),
    levels = c("c1", "c2", "c4", "c3", "c5")
)
names(CAF_NIF@active.ident) <- rownames(CAF_NIF@meta.data)

Randomly_Selected_Cells <- c()
for (i in c("c1", "c2", "c3", "c4", "c5")) {
    CAF_NIF <- subset(CAF_NIF, ident = i)
    set.seed(1)
    Randomly_Selected_Cells <- append(
        Randomly_Selected_Cells,
        rownames(CAF_NIF@meta.data[sample(1:nrow(CAF_NIF@meta.data), 1500), ]) # nolint
    )
}

CAF_NIF_Randomly <- subset(data, cells = Randomly_Selected_Cells)
CAF_NIF_Randomly$group <- ifelse(
    CAF_NIF_Randomly@active.ident %in% c(1, 2, 4),
    "CAF", "NIF"
)

monocds <- RunMonocle(CAF_NIF_Randomly)
monocds <- orderCells(monocds, root_state = 2)

pdf("CAF_NIF_Trajactory_State.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
    monocds,
    cell_size = 1, color_by = "State",
    show_tree = FALSE, show_branch_points = FALSE
) +
    scale_color_manual(
        values = c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
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

pdf("CAF_NIF_Trajactory_Label.pdf", width = 4, height = 4)
plot_cell_trajectory(monocds,
    cell_size = 1,
    color_by = "label", show_tree = FALSE, show_branch_points = FALSE
) +
    scale_color_manual(values = colo[c(1, 4, 3, 9)])
dev.off()

## The percentage of clusters in each state
mat <- melt(table(monocds@phenoData@data[, c("State", "label")]))

p <- ggplot(mat, aes(x = State, y = value, fill = label)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(y = "Proportion (%)") +
    scale_fill_manual(values = pal_npg()(4)) +
    coord_flip() +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.6, "cm"),
        axis.text.x = element_text(
            color = "black", size = 10
        ),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = percent)
ggsave("CAF_NIF_Trajactory_Proportion.pdf", p, height = 2.5, width = 5)

beam_res <- BEAM(monocds, branch_point = 1, cores = 40)
beam_res <- beam_res[order(beam_res$qval), ]
beam_res <- beam_res[, c("gene_short_name", "pval", "qval")]

pdf("CAF_NIF_Trajactory_Heatmap.pdf", height = 8, width = 6)
plot_genes_branched_heatmap(monocds[row.names(subset(
    beam_res,
    qval < 1e-4
)), ],
branch_point = 1,
num_clusters = 4,
cores = 30,
use_gene_short_name = T,
show_rownames = F
)
dev.off()

p <- plot_genes_branched_heatmap(monocds[row.names(subset(
    beam_res,
    qval < 1e-4
)), ],
branch_point = 1,
num_clusters = 4,
cores = 30,
use_gene_short_name = T,
show_rownames = T,
return_heatmap = T
)
write.table(
    p$annotation_row, "CAF_NIF_Trajactory_Heatmap_Label.txt",
    quote = F, sep = "\t"
)

hallmark <- getGeneSets(library = "H")
data_set <- list(enriched = hallmark[[14]]@geneIds)
EMT_Scores <- gsva(
    as.matrix(CAF_NIF_Randomly@assays$RNA@data), data_set,
    min.sz = 5, kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 20
)

monocds@phenoData@data$EMT <- scale(as.vector(t(EMT_Scores)))
monocds@phenoData@data$EMT <- ifelse(
    monocds@phenoData@data$EMT > 1.5, 1.5, monocds@phenoData@data$EMT
)

pdf("CAF_NIF_Trajactory_EMT.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
    monocds,
    cell_size = 0.5, color_by = "EMT",
    show_tree = FALSE, show_branch_points = FALSE
) +
    scale_color_gradientn(colors = rev(col1)) +
    theme(
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
dev.off()

EMT_Value <- data.frame(
    State = monocds@phenoData@data$State,
    EMT = as.vector(t(EMT_Scores)),
    subtype = monocds@phenoData@data$tissue
)

p <- ggplot(EMT_Value, aes(x = State, y = EMT)) +
    geom_jitter(aes(fill = State, color = State),
        width = 0.2, shape = 21, size = 0.5
    ) +
    scale_color_manual(
        values = c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
    ) +
    scale_fill_manual(
        values = c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
    ) +
    geom_boxplot(
        size = 0.6, fill = "white",
        outlier.fill = NA, outlier.color = NA, outlier.size = 0
    ) +
    theme_classic() +
    labs(title = "EMT scores") +
    theme(
        axis.title.y = element_blank(),
        axis.text = element_text(size = 11, color = "black"),
        axis.ticks.length = unit(0.4, "lines"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)
    )
ggsave("CAF_NIF_Trajactory_EMT_BoxPlot.pdf", p, width = 2.8, height = 3.3)

CREB3L1_Matrix <- data.frame(
    t(monocds@reducedDimS), CAF_NIF_Randomly@assays$RNA@data["CREB3L1", ]
)
colnames(CREB3L1_Matrix) <- c("Component1", "Component2", "CREB3L1")
CREB3L1_Matrix$group <- monocds@phenoData@data$State
CREB3L1_Matrix <- CREB3L1_Matrix[order(mat$CREB3L1), ]
pdf("creb3l1.pdf", width = 5, height = 3)
ggplot(CREB3L1_Matrix, aes(x = Component1, y = Component2, color = CREB3L1)) +
    geom_point(size = 0.5) +
    scale_color_gradientn(colors = rev(color_for_use)) +
    theme_classic() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
    )
dev.off()

CREB3L1_Matrix$group <- paste0("State", CREB3L1_Matrix$group)
pdf("CAF_NIF_Trajactory_CREB3L1_Box.pdf", width = 2, height = 3)
ggplot(CREB3L1_Matrix, aes(x = group, y = CREB3L1, fill = group)) +
    geom_boxplot(outlier.color = "white", outlier.size = 0) +
    theme_classic() +
    scale_fill_manual(
        values =
            c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
    ) +
    labs(x = "", y = "Expression levels of CREB3L1") +
    theme(legend.position = "none")
dev.off()

CREB3L1_Matrix$tissue <- CAF_NIF_Randomly$tissue

p <- ggplot(CREB3L1_Matrix, aes(x = tissue, y = CREB3L1, fill = group)) +
    scale_fill_manual(values = c(pal_aaas()(10), rev(pal_aaas()(5)))) +
    geom_boxplot(
        size = 0.6,
        outlier.fill = NA, outlier.color = NA, outlier.size = 0
    ) +
    theme_bw() +
    labs(title = "Expression levels of CREB3L1") +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(
            size = 11, angle = 45, hjust = 1, color = "black"
        ),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
    )
ggsave("e_state_creb3l1.pdf", p, width = 5.5, height = 3.3)

CREB3L1_Target_Genes <- read.table(
    "CREB3L1_Target_Genes.txt",
    header = T, sep = "\t"
)
CREB3L1_Target_Genes_Enrichment <- gsva(
    as.matrix(CAF_NIF_Randomly@assays$RNA@data),
    list(CREB3L1_Target_Genes$gene),
    min.sz = 5, kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 20
)
monocds@phenoData@data$target <- scale(
    as.vector(t(CREB3L1_Target_Genes_Enrichment))
)
monocds@phenoData@data$target <- ifelse(
    monocds@phenoData@data$target > 1.5, 1.5, monocds@phenoData@data$target
)

pdf("CREB3L1_Target_Genes_Enrichment.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
    monocds,
    cell_size = 0.5, color_by = "target",
    show_tree = FALSE, show_branch_points = FALSE
) +
    scale_color_gradientn(colors = rev(col1)) +
    theme(
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
dev.off()

CREB3L1_Target_Value <- data.frame( # nolint
    target = scale(as.vector(t(CREB3L1_Target_Genes_Enrichment))),
    t(monocds@reducedDimS),
    group = paste0("State", monocds@phenoData@data$State)
)
colnames(CREB3L1_Target_Value)[2:3] <- c("Component1", "Component2")
pdf("CREB3L1_Target_BoxPlot.pdf", width = 2, height = 2.6)
ggplot(CREB3L1_Target_Value, aes(x = group, y = target, fill = group)) +
    geom_boxplot(outlier.color = "white", outlier.size = 0) +
    theme_classic() +
    scale_fill_manual(
        values = c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
    ) +
    labs(x = "", y = "Regulon activity of CREB3L1") +
    theme(
        legend.position = "none",
        axis.title.y = element_text(size = 7, color = "black")
    )
dev.off()
