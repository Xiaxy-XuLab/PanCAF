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

## Subcluster analysis for Lymphocyte
lymphocyte_clusters <- subset(data, ident = 9:17)
lymphocyte_clusters %<>% NormalizeData(object = ., normalization.method = "LogNormalize") %>% # nolint
    FindVariableFeatures(selection.method = "vst")

all.genes <- rownames(lymphocyte__clusters)
lymphocyte_clusters %<>% ScaleData(object = ., features = all.genes) %>%
    ScaleData(object = ., vars.to.regress = "percent.mt")
lymphocyte_clusters %<>% RunPCA(
    object = .,
    features = VariableFeatures(object = .)
) %>%
    RunFastMNN(object = ., object.list = SplitObject(., split.by = "SampleID"))

lymphocyte_clusters %<>% FindNeighbors(reduction = "mnn", dims = 1:30) %>%
    FindClusters(resolution = 1.2) %<>%
    RunTSNE(reduction = "mnn", dims = 1:30) %>%
    RunUMAP(reduction = "mnn", dims = 1:30)

saveRDS(lymphocyte_clusters, "lymphocyte_clusters.Rds")

lymphocyte_clusters_anno <- subset(lymphocyte_clusters, ident = c(0:12, 14, 16))
lymphocyte_clusters_anno <- RenameIdents(lymphocyte_clusters_anno,
    "0" = "CD8+ TEM", "5" = "CD8+ TEM", "8" = "CD8+ TEM",
    "10" = "CD8+ TEM", "1" = "CD4+ TCM", "7" = "CD4+ TCM",
    "4" = "Treg", "14" = "Cycling CD4+ T", "2" = "DN T",
    "9" = "CD8+ CTL", "3" = "PDCD1+ CD8+ T",
    "12" = "Cycling CD8+ T", "16" = "Cycling NKC",
    "6" = "TIM3+ NKC", "11" = "CD160+ NKC"
)

pdf("Lymphocyte_Umap.pdf", width = 4, height = 4)
DimPlot(lymphocyte_clusters_anno, label = T) + NoLegend() +
    scale_color_manual(values = c(pal_npg()(9), pal_jama()(5))) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )
dev.off()

Lymphocyte_Marker <- FindAllMarkers(
    ss,
    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)
write.table(Lymphocyte_Marker, "Lymphocyte_Marker.txt", quote = F, sep = "\t")

Select_Genes <- c()
for (i in unique(Lymphocyte_Marker$cluster)) {
    Select_Genes <- append(
        Select_Genes,
        marker[which(Lymphocyte_Marker$cluster == i), "gene"][1:5]
    )
}

lymphocyte_clusters_anno@active.ident <- factor(
    lymphocyte_clusters_anno@active.ident,
    levels = rev(levels(lymphocyte_clusters_anno))
)
pdf("Lymphocyte_.pdf", width = 11, height = 4.5)
DotPlot(lymphocyte_clusters_anno,
    features = unique(Select_Genes),
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

lymphocyte_clusters_anno$seurat_clusters <- lymphocyte_clusters_anno@active.ident # nolint
lymphocyte_table <- table(
    lymphocyte_clusters_anno@meta.data[, c("group", "seurat_clusters")]
)
lymphocyte_table <- melt(lymphocyte_table / as.vector(table(data$group)))
lymphocyte_table$group <- factor(
    lymphocyte_table$group,
    levels = c("Normal", "Adjacent", "Tumor")
)
p <- ggplot(
    lymphocyte_table,
    aes(x = seurat_clusters, y = value, fill = group)
) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic() +
    coord_flip() +
    scale_fill_manual(
        values = c(
            pal_aaas()(9)[c(1, 3, 2)]
        )
    ) +
    labs(y = "Proportion", x = "") +
    theme(
        axis.line = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        legend.key.size = unit(0.6, "cm")
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = percent)

ggsave("Lymphocyte_Clusters_Percent.pdf", p, width = 5, height = 3.5)

## CAF vs DC sub-clusters
CAF_DC <- read.table("GO/go4/r11/out/count_network.txt", header = T, sep = "\t")
set.seed(15)
pdf("CAF_DC.pdf", width = 3.5, height = 3.5)
chordDiagram(CAF_DC, transparency = 0.5)
dev.off()

mypvals <- read.delim("GO/go4/r11/out/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("GO/go4/r11/out/means.txt", check.names = FALSE)

costimulatory1 <- read.table("lire/DC-CAF.txt", header = F, sep = "\t")
costimulatory2 <- read.table("lire/CAF-DC.txt", header = F, sep = "\t")

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory2$V1) %>%
    dplyr::select("interacting_pair", starts_with("state")) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("[|]state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}

mat <- mat[-nrow(mat), ]
annotation_col <- data.frame(
    celltype = rep(paste0("state2", "state3"), each = 6),
    group = rep(paste0("dc", 0:5), time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- pal_aaas()(6)
names(col2) <- paste0("dc", 0:5)
ann_colors <- list(
    group = col2,
    celltype = col1
)

bk <- c(seq(0, 4, by = 0.01))
pdf("CAF-DC_Ligand_Receptor1.pdf", width = 6.5, height = 3)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory1$V1) %>%
    dplyr::select("interacting_pair", grep("state[1-3]$", names(mypvals))) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("^state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}
ord <- colnames(mat)[c(
    grep("state1", colnames(mat)),
    grep("state2", colnames(mat)),
    grep("state3", colnames(mat))
)]
mat <- mat[, ord]
annotation_col <- data.frame(
    celltype = rep(paste0("state", 1:3), each = 6),
    group = rep(paste0("dc", 0:5), time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- pal_aaas()(6)
names(col2) <- paste0("dc", 0:5)
ann_colors <- list(
    group = col2,
    celltype = col1
)

bk <- c(seq(0, 4, by = 0.01))
pdf("DC-CAF_Ligand_Receptor2.pdf", width = 6.5, height = 3)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

## CAF vs B sub-clusters
CAF_B <- read.table("GO/go4/r10/out/count_network.txt", header = T, sep = "\t")
set.seed(15)
pdf("caf_B.pdf", width = 3.5, height = 3.5)
chordDiagram(CAF_B, transparency = 0.5)
dev.off()

mypvals <- read.delim("GO/go4/r10/out/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("GO/go4/r10/out/means.txt", check.names = FALSE)

costimulatory1 <- read.table("lire/CAF-B.txt", header = F, sep = "\t")

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory1$V1) %>%
    dplyr::select("interacting_pair", starts_with("state")) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("[|]state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}

annotation_col <- data.frame(
    celltype = rep(paste0("state", 1:3), each = 8),
    group = rep(paste0("B_sc", 0:7), time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- pal_aaas()(8)
names(col2) <- paste0("B_sc", 0:7)
ann_colors <- list(
    group = col2,
    celltype = col1
)

bk <- c(seq(0, 4, by = 0.01))
pdf("CAF_B_Ligand_Receptor1.pdf", width = 6.5, height = 3)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory1$V1) %>%
    dplyr::select("interacting_pair", grep("state[1-3]$", names(mypvals))) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("^state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}
ord <- colnames(mat)[c(
    grep("state1", colnames(mat)),
    grep("state2", colnames(mat)),
    grep("state3", colnames(mat))
)]
mat <- mat[, ord]

annotation_col <- data.frame(
    celltype = rep(paste0("state", 1:3), each = 8),
    group = rep(paste0("B_sc", 0:7), time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- pal_aaas()(8)
names(col2) <- paste0("B_sc", 0:7)
ann_colors <- list(
    group = col2,
    celltype = col1
)

bk <- c(seq(0, 4, by = 0.01))
pdf("B_CAF_Ligand_Receptor2.pdf", width = 6.5, height = 3)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

## CAF vs TEC/NEC
CAF_TECNEC <- read.table("GO/go4/r2/out/count_network.txt", header = T, sep = "\t") # nolint
set.seed(15)
pdf("CAF_TECNEC.pdf", width = 3.5, height = 3.5)
chordDiagram(CAF_TECNEC, transparency = 0.5)
dev.off()

mypvals <- read.delim("GO/go4/r2/out/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("GO/go4/r2/out/means.txt", check.names = FALSE)

costimulatory1 <- read.table("lire/TEC-CAF.txt", header = F, sep = "\t")
costimulatory2 <- read.table("lire/CAF-TEC.txt", header = F, sep = "\t")

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory2$V1) %>%
    dplyr::select("interacting_pair", starts_with("state")) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("[|]state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}

annotation_col <- data.frame(
    celltype = rep(paste0("state", 1:3), each = 2),
    group = rep(c("NEC", "TEC"), time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- pal_aaas()(2)
names(col2) <- c("NEC", "TEC")
ann_colors <- list(
    group = col2,
    celltype = col1
)

bk <- c(seq(0, 4, by = 0.01))
pdf("CAF_TECNEC_Ligand_Receptor1.pdf", width = 4.5, height = 3)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory1$V1) %>%
    dplyr::select("interacting_pair", grep("state[1-3]$", names(mypvals))) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("^state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}
ord <- colnames(mat)[c(
    grep("state1", colnames(mat)),
    grep("state2", colnames(mat)),
    grep("state3", colnames(mat))
)]
mat <- mat[, ord]
annotation_col <- data.frame(
    celltype = rep(paste0("state", 1:3), each = 2),
    group = rep(c("NEC", "TEC"), time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- pal_aaas()(2)
names(col2) <- c("NEC", "TEC")
ann_colors <- list(
    group = col2,
    celltype = col1
)

bk <- c(seq(0, 4, by = 0.01))
pdf("TECNEC_CAF_Ligand_Receptor2.pdf", width = 4.5, height = 3.5)
pheatmap(mat,
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

## CAF vs NKT sub-clusters
CAF_NKT <- read.table("figure3/3C/out/count_network.txt", header = T, sep = "\t") # nolint
set.seed(15)
pdf("CAF_NKT.pdf", width = 3.5, height = 3.5)
chordDiagram(CAF_NKT, transparency = 0.5)
dev.off()

mypvals <- read.delim("figure3/3C/out/pvalues.txt", check.names = FALSE)
mymeans <- read.delim("figure3/3C/out/means.txt", check.names = FALSE)

costimulatory1 <- read.table("lire/NKT-CAF.txt", header = F, sep = "\t")
costimulatory2 <- read.table("lire/CAF-NKT.txt", header = F, sep = "\t")

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory2$V1) %>%
    dplyr::select("interacting_pair", starts_with("state")) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("[|]state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}

annotation_col <- data.frame(
    celltype = rep(c("state1", "state2", "state3"), each = 11),
    group = rep(sapply(strsplit(colnames(mat), "|", fixed = T), "[", 2)[1:11], time = 3)
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- c("state1", "state2", "state3")
col2 <- c(pal_aaas()(9), "#C7A162", "black")
names(col2) <- sapply(strsplit(colnames(mat), "|", fixed = T), "[", 2)[1:11]
ann_colors <- list(
    group = col2,
    celltype = col1
)

annotation_col$celltype <- factor(
    annotation_col$celltype,
    levels = paste0("state", 1:3)
)
annotation_col$group <- factor(
    annotation_col$group,
    levels = annotation_col$group[1:11]
)
annotation_col <- annotation_col[
    order(
        annotation_col$group,
        annotation_col$celltype
    ),
]

bk <- c(seq(0, 4, by = 0.01))
pdf("CAF-NKT_Ligand_Receptor1.pdf", width = 8, height = 4)
pheatmap(mat[, rownames(annotation_col)],
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

mypvals %>%
    dplyr::filter(interacting_pair %in% costimulatory1$V1) %>%
    dplyr::select("interacting_pair", grep("state[1-3]$", names(mypvals))) %>%
    reshape2::melt() -> pvalsdf
pvalsdf <- pvalsdf[-grep("^state", pvalsdf$variable), ]
colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
mat <- dcast(pvalsdf, CC ~ interacting_pair)
rownames(mat) <- mat$CC
mat <- t(mat[, -1])
mat <- -log10(mat + 0.0001)
for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
        if (mat[i, j] < 0) {
            mat[i, j] <- 0
        }
    }
}
ord <- colnames(mat)[
    c(
        grep("state1", colnames(mat)),
        grep("state2", colnames(mat)),
        grep("state3", colnames(mat))
    )
]
mat <- mat[, ord]

annotation_col <- data.frame(
    celltype = rep(paste0("state", 1:3), each = 11),
    group = rep(
        unique(sapply(strsplit(colnames(mat), "|", fixed = T), "[", 1)),
        time = 3
    )
)
rownames(annotation_col) <- colnames(mat)

col1 <- pal_npg("nrc")(3)
names(col1) <- paste0("state", 1:3)
col2 <- c(pal_aaas()(9), "#C7A162", "black")
names(col2) <- unique(sapply(strsplit(colnames(mat), "|", fixed = T), "[", 1))
ann_colors <- list(
    group = col2,
    celltype = col1
)
annotation_col$celltype <- factor(
    annotation_col$celltype,
    levels = paste0("state", 1:3)
)
annotation_col$group <- factor(
    annotation_col$group,
    levels = annotation_col$group[1:11]
)
annotation_col <- annotation_col[
    order(annotation_col$group, annotation_col$celltype),
]

bk <- c(seq(0, 4, by = 0.01))
pdf("NKT-CAF_Ligand_Receptor2.pdf", width = 8, height = 4)
pheatmap(mat[, rownames(annotation_col)],
    show_colnames = FALSE, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, cluster_cols = FALSE,
    cluster_rows = T, color = colorRampPalette(
        colors = rev(c(brewer.pal(11, "RdGy")[1:6], "white"))
    )(length(bk)),
    annotation_legend = TRUE, treeheight_row = 5,
    fontsize = 12, annotation_names_col = FALSE, border_color = "gray",
)
dev.off()

## survival analysis
## UCC_CAF state3
data <- read.table("UCC.txt", sep = "\t", header = T)

data <- data[order(data$Fibroblast_3), ]
mt <- rbind(
    data.frame(data[1:270, ],
        group = rep("LOW", 270)
    ),
    data.frame(data[(nrow(data) - 77):nrow(data), ],
        group = rep("HIGH", 78)
    )
)

mt$event <- ifelse(mt$event == 0, 1, 2)
fit <- survfit(Surv(time, event) ~ group, data = mt)
p <- ggsurvplot(fit,
    data = mt, pval = TRUE, conf.int = TRUE,
    risk.table = TRUE, risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    palette = c("red", "blue")
)

pdf("UCC_CAF_State3.pdf")
print(p, newpage = FALSE)
dev.off()

### UUS_CAF state3
data <- read.table("UUS.txt", sep = "\t", header = T)

data <- data[order(data$Fibroblast_3), ]
mt <- rbind(
    data.frame(data[1:16, ],
        group = rep("LOW", 16)
    ),
    data.frame(data[(nrow(data) - 33):nrow(data), ],
        group = rep("HIGH", 34)
    )
)
mt$event <- ifelse(mt$event == 0, 1, 2)
fit <- survfit(Surv(time, event) ~ group, data = mt)
p <- ggsurvplot(fit,
    data = mt, pval = TRUE, conf.int = TRUE,
    risk.table = TRUE, risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    palette = c("red", "blue")
)
pdf("UUS_CAF_State3.pdf")
print(p, newpage = FALSE)
dev.off()

####### Melanoma_CAF state3
data <- read.table("melanoma_liu.txt", sep = "\t", header = T)

data <- data[order(data$Fibroblast_3), ]
mt <- rbind(
    data.frame(data[1:41, ],
        group = rep("LOW", 41)
    ),
    data.frame(data[(nrow(data) - 79):nrow(data), ],
        group = rep("HIGH", 80)
    )
)

mt$event <- ifelse(mt$event == 0, 1, 2)
fit <- survfit(Surv(time, event) ~ group, data = mt)
p <- ggsurvplot(fit,
    data = mt, pval = TRUE, conf.int = TRUE,
    risk.table = TRUE, risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    palette = c("red", "blue")
)
pdf("Melanoma_CAF_State3.pdf")
print(p, newpage = FALSE)
dev.off()
