pkgs <- c(
    "Seurat", "SeuratWrappers", "ggplot2", "batchelor", "circlize",
    "dplyr", "optparse", "reshape2", "data.table", "magrittr",
    "patchwork", "scales", "GSVA", "RColorBrewer", "ggridges",
    "clusterProfiler", "survminer", "survminer", "monocle", "nichenet",
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

## FeaturePlot
subset <- subset(data, ident = 6)
pdf("FeaturePlot_Markers.pdf", width = 10, height = 10)
FeaturePlot(subset,
    features =
        c("RGS5", "ACTA2", "PLVAP", "VWF"), order = T, ncol = 2
) &
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
dev.off()

pdf("VlnPlot_Markers.pdf", width = 6, height = 6)
VlnPlot(subset,
    features = c("RGS5", "PLVAP"),
    pt.size = 0, ncol = 2
) &
    scale_fill_manual(values = pal_npg()(10)) &
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
dev.off()

## TME interaction in normal, adjacent and tumor
subsets <- subset(data, ident = c(1:34))
inter <- function(x) {
    subset_select <- subset(
        subsets,
        cells = rownames(subsets@meta.data[which(subsets$group == x), ])
    )

    Select_Cells <- c()
    for (i in levels(subset_select)) {
        subset <- subset(subset_select, ident = i)
        set.seed(1)
        if (nrow(subset@meta.data) > 1000) {
            Select_Cells <- append(
                Select_Cells,
                rownames(subset@meta.data[sample(1:nrow(subset@meta.data), 1000), ]) # nolint
            )
        } else {
            Select_Cells <- append(Select_Cells, rownames(subset@meta.data))
        }
    }

    Object_Select <- subset(subset_select, cells = Select_Cells)
    Object_Select$seurat_clusters <- Object_Select@active.ident

    write.table(
        Object_Select@assays$RNA@counts,
        paste0("figure5/1A/", x, "/count.txt"),
        quote = F, sep = "\t"
    )

    meta <- data.frame(
        name = cellname,
        cluster = paste0("c", Object_Select@meta.data[cellname, "Clusters"])
    )
    write.table(
        meta,
        paste0("figure5/1A/", x, "/meta.txt"),
        quote = F, sep = "\t",
        row.names = F
    )
}
inter("Tumor")
inter("Adjacent")
inter("Normal")

Normal <- read.table(
    "figure5/1A/Normal/out/count_network.txt",
    header = T, sep = "\t"
)
Adjacent <- read.table(
    "figure5/1A/Adjacent/out/count_network.txt",
    header = T, sep = "\t"
)
Tumor <- read.table(
    "figure5/1A/Tumor/out/count_network.txt",
    header = T, sep = "\t"
)

redu <- function(x) {
    return(x[which(x$SOURCE == "c6"), ])
}

mat <- cbind(redu(nature), redu(normal), redu(tumor))
mat <- mat[, c(1, 2, 3, 6, 9)]
colnames(mat) <- c("source", "target", "Normal", "Adjacent", "Tumor")
mt <- melt(mat)

model <- aov(value ~ variable, data = mt)
rht <- glht(model, linfct = mcp(variable = "Dunnett"), alternative = "two.side")
summary(rht)

p <- ggplot(mt, aes(variable, value)) +
    geom_boxplot(
        size = 0.8, fill = "white",
        outlier.fill = NA, outlier.color = NA, outlier.size = 0
    ) +
    geom_point(aes(color = variable), size = 6) +
    geom_line(aes(group = target), size = 0.6, colour = "#9C9C9C") +
    theme_classic() +
    scale_color_manual(values = c("#9362cc", "#5199cc", "#fe9929")) +
    labs(y = "Interaction counts") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(
            size = 14, color = "black",
            angle = 45, hjust = 1
        ),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(0.4, "lines"),
        # legend.position = "none",
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    )
ggsave("TME_Interaction_Boxplot.pdf", p, width = 3.5, height = 4)

Tumor <- dcast(Tumor, SOURCE ~ TARGET)
bk <- c(seq(0, 5, by = 0.01))
pdf("TME_Interaction_Heatmap_Tumor.pdf", width = 8, height = 8)
p <- pheatmap(log1p(Tumor),
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

Cluster_Order <- rownames(Tumor)[p$tree_row$order]
pdf("TME_Interaction_Heatmap_Adjacent.pdf", width = 8, height = 8)
p <- pheatmap(log1p(Adjacent[Cluster_Order, Cluster_Order]),
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

pdf("TME_Interaction_Heatmap_Normal.pdf", width = 5.5, height = 5)
p <- pheatmap(log1p(Normal[Cluster_Order, Cluster_Order]),
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

## Trajactory analysis
Select_Cells <- c()
for (i in c(1, 26)) {
    subsets <- subset(data, ident = i)
    set.seed(1)
    Select_Cells <- append(
        Select_Cells,
        rownames(subsets@meta.data[sample(1:nrow(subsets@meta.data), 3000), ])
    )
}
cluster6 <- subset(data, ident = 6)
cluster6_cells <- rownames(cluster6@meta.data)
monotj <- subset(data, cells = c(Select_Cells, cluster6_cells))
RunMonocle(monotj)

monocds <- orderCells(monocds, root_state = 3)
pdf("Trajactory_State.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "State",
    show_tree = FALSE, show_branch_points = FALSE
) +
    scale_color_manual(values = pal_aaas()(9)) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
    )
dev.off()

pdf("Trajactory_Cluster.pdf", width = 3.5, height = 4)
plot_cell_trajectory(monocds,
    cell_size = 0.2,
    color_by = "groo", show_tree = FALSE, show_branch_points = FALSE
) +
    scale_color_manual(values = colo[c(1, 2, 3, 4)]) +
    theme(
        legend.title = element_blank(),
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
    fullModelFormulaStr = "~sm.ns(Pseudotime)",
    cores = 40
)
sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 10^-20))

pdf("Trajactory_Pheatmap.pdf", width = 6, height = 8)
plot_pseudotime_heatmap(
    monocds[sig_gene_names, ],
    num_clusters = 3, cores = 40,
    use_gene_short_name = TRUE,
    show_rownames = FALSE
)
dev.off()

cds_subset <- monocds[c("ITGA9", "ITGB1"), ]
cds_subset@phenoData@data$clusters <- paste0("c", cds_subset@phenoData@data$Clusters) # nolint
pdf("Trajactory_Genes.pdf", width = 5, height = 3)
plot_genes_in_pseudotime(cds_subset, color_by = "Clusters") +
    scale_color_manual(values = pal_npg()(3))
dev.off()

gshallmark <- getGeneSets(library = "H")

data_set <- list(enriched = gshallmark@.Data[[4]]@geneIds)
result <- gsva(
    as.matrix(monotj@assays$RNA@data), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
dat <- data.frame(
    group = paste0("s", monocds@phenoData@data$State), value = as.vector(result)
)
dat$groupp <- ifelse(dat$group %in% c("s1", "s7"), "state1",
    ifelse(dat$group %in% c("s4", "s5"), "state2", "state3")
)
p <- ggplot(dat, aes(x = groupp, y = value)) +
    geom_jitter(aes(fill = groupp, color = groupp),
        width = 0.2, shape = 21, size = 0.5
    ) +
    scale_color_npg(alpha = 0.9) +
    scale_fill_npg(alpha = 0.9) +
    geom_boxplot(
        size = 0.6, fill = "white", outlier.fill = NA,
        outlier.color = NA, outlier.size = 0
    ) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        legend.position = "none"
    )

ggsave("Trajactory_State_Angiogenesis.pdf", p, width = 2.5, height = 3)

monocds@phenoData@data$ANGIOGENESIS <- scale(as.vector(result))
monocds@phenoData@data$ANGIOGENESIS <- ifelse(
    monocds@phenoData@data$ANGIOGENESIS > 2, 2,
    ifelse(monocds@phenoData@data$ANGIOGENESIS < -2, -2,
        monocds@phenoData@data$ANGIOGENESIS
    )
)

pdf("Trajactory_Angiogenesis.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(monocds,
    cell_size = 1,
    color_by = "ANGIOGENESIS", show_tree = FALSE,
    show_branch_points = FALSE
) +
    scale_color_gradientn(colors = rev(col1)) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
dev.off()

## Nichenet analysis
dd <- subset(data, ident = c(25, 27, 6))
dd$gro <- ifelse(dd$seurat_clusters == 6, "CAFendomt", "Endothelium")

dd@active.ident <- factor(dd$gro)
marker <- FindAllMarkers(dd, min.pct = 0.25, only.pos = T)
gg <- marker[which(marker$cluster == "CAFendomt"), ]
colo <- c(
    pal_npg(alpha = 0.9)(9),
    pal_aaas(alpha = 0.7)(9), pal_lancet(alpha = 0.5)(9)
)

gene <- read.table("CAF_Markers.txt", header = T, sep = "\t")
gene <- gg[which(gg$p_val_adj == 0), "gene"]


gene <- read.table("CAF_Markers.txt", header = T, sep = "\t")
ligand_target_matrix <- readRDS("/work/xiaxy/learn/7.nechenet/ligand_target_matrix.rds")
lr_network <- readRDS("/work/xiaxy/learn/7.nechenet/lr_network.rds")
weighted_networks <- readRDS("/work/xiaxy/learn/7.nechenet/weighted_networks.rds")
weighted_networks_lr <- weighted_networks$lr_sig %>%
    inner_join(lr_network %>% distinct(from, to), by = c("from", "to"))

receiver <- "6"
receobject <- subset(data, ident = receiver)
expressed_genes_receiver <- get_expressed_genes(
    receiver, receobject,
    pct = 0.10
)
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] # nolint

sender_celltypes <- 19
sendobject <- subset(data, ident = sender_celltypes)
list_expressed_genes_sender <- get_expressed_genes(sender_celltypes, sendobject, pct = 0.10) # nolint
expressed_genes_sender <- list_expressed_genes_sender %>%
    unlist() %>%
    unique()

geneset_oi <- gene$gene
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands <- lr_network %>%
    pull(from) %>%
    unique()
receptors <- lr_network %>%
    pull(to) %>%
    unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()

ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
)
ligand_activities <- ligand_activities %>%
    arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)))

best_upstream_ligands <- ligand_activities %>%
    top_n(30, pearson) %>%
    arrange(-pearson) %>%
    pull(test_ligand) %>%
    unique()

pdf("Ligands_Expr.pdf", width = 9, height = 5)
DotPlot(sendobject,
    features = best_upstream_ligands %>% rev(),
    cols = "RdYlBu"
) +
    RotatedAxis() + coord_flip() +
    theme(axis.title = element_blank())
dev.off()

active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
        geneset = geneset_oi,
        ligand_target_matrix = ligand_target_matrix, n = 200
    ) %>%
    bind_rows() %>%
    drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix, cutoff = 0.33
)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% # nolint
    rev() %>%
    make.names()
order_targets <- active_ligand_target_links_df$target %>%
    unique() %>%
    intersect(rownames(active_ligand_target_links)) %>%
    make.names()
rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>%
    make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>%
    make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t() # nolint

p_ligand_target_network <- vis_ligand_target %>%
    make_heatmap_ggplot("Prioritized ligands", "Predicted target genes",
        color = "purple", legend_position = "top",
        x_axis_position = "top", legend_title = "Regulatory potential"
    ) +
    theme(axis.text.x = element_text(face = "italic")) +
    scale_fill_gradient2(
        low = "whitesmoke",
        high = "purple", breaks = c(0, 0.0045, 0.0090)
    )
ggsave(
    "Ligand_Target_Network.pdf", p_ligand_target_network,
    width = 10, height = 5
)

vis_ligand_target <- vis_ligand_target[, c(
    "CD44", "COL1A1", "TGFBR2", "COL3A1",
    "ZEB2", "PALLD", "EMP1", "C1S", "BGN", "MFGE8"
)]
bk <- c(seq(0, 0.012, by = 0.0001))
colo <- colorRampPalette(c(
    "#324C74", "#377D8C",
    "#6CB39D", "#B2DDA6", "#F9FEB7"
))(50)
p2 <- pheatmap(vis_ligand_target,
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo[50:48])(length(bk) / 2),
        colorRampPalette(colors = colo[47:1])(length(bk) / 2)
    ))
)
ggsave("Lt-Net.pdf", p2, width = 4.5, height = 4)
gge <- rownames(vis_ligand_target)
gge <- gsub("[.]", "-", gge)
sub <- subset(data, ident = 4)
mat <- sub@assays$RNA@data[gge, ]
mt <- rowMeans(mat)
tt <- as.data.frame(mt)

bk <- c(seq(0, 3, by = 0.01))
colo2 <- colorRampPalette(c("white", "#18196A"))(50)
p1 <- pheatmap(tt,
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo2[50:45])(length(bk) / 2),
        colorRampPalette(colors = colo2[44:1])(length(bk) / 2)
    ))
)
ggsave("Lt_Expr.pdf", p1, width = 1.5, height = 4)

aab <- colnames(vis_ligand_target)
aab <- gsub("[.]", "-", aab)
sb <- subset(data, ident = 56)
mat <- sb@assays$RNA@data[aab, ]
mt <- rowMeans(mat)
tt <- as.data.frame(t(mt))

bk <- c(seq(0, 1.4, by = 0.01))
colo2 <- colorRampPalette(c("white", "#F8DB85"))(50)
p1 <- pheatmap(tt,
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo2[50:25])(length(bk) / 2),
        colorRampPalette(colors = colo2[25:0])(length(bk) / 2)
    ))
)
ggsave("Lt_Expression.pdf", p1, width = 4, height = 1)

layout <- "AABBBBBB"
final <- as.ggplot(p1) + as.ggplot(p2) +
    plot_layout(design = layout) &
    theme(legend.position = "bottom")
ggsave("Nichenet_Merge.pdf", final, width = 6, height = 3.5)

lr_network_top <- lr_network %>%
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
    distinct(from, to)
best_upstream_receptors <- lr_network_top %>%
    pull(to) %>%
    unique()

lr_network_top_df_large <- weighted_networks_lr %>%
    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df <- lr_network_top_df_large %>%
    spread("from", "weight", fill = 0)
lr_network_top_matrix <- lr_network_top_df %>%
    select(-to) %>%
    as.matrix() %>%
    magrittr::set_rownames(lr_network_top_df$to)

dist_receptors <- dist(lr_network_top_matrix, method = "binary")
hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
order_receptors <- hclust_receptors$labels[hclust_receptors$order]

dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix)) # nolint
order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix)) # nolint

vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor] # nolint
rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()
aa <- vis_ligand_receptor_network %>% t()
bk <- c(seq(0, 1.25, by = 0.01))
colo3 <- colorRampPalette(c("white", "#7E294E"))(50)
p1 <- pheatmap(aa[gge, which(colSums(aa[gge, ]) > 0.1)],
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo3[50:45])(length(bk) / 2),
        colorRampPalette(colors = colo3[44:1])(length(bk) / 2)
    ))
)
ggsave("LT_LP.pdf", p1, width = 8, height = 5)
p_ligand_receptor_network <- vis_ligand_receptor_network %>%
    t() %>%
    make_heatmap_ggplot("Ligands", "Receptors",
        color = "#7E294E",
        x_axis_position = "top",
        legend_title = "Prior interaction potential"
    )
ggsave("Network.pdf", p_ligand_receptor_network, width = 6, height = 5)

ligand_pearson_matrix <- ligand_activities %>%
    select(pearson) %>%
    as.matrix() %>%
    magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) <- rownames(ligand_pearson_matrix) %>%
    make.names()
colnames(ligand_pearson_matrix) <- colnames(ligand_pearson_matrix) %>%
    make.names()

vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>%
    as.matrix(ncol = 1) %>%
    magrittr::set_colnames("Pearson")
vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>%
    as.matrix(ncol = 1) %>%
    magrittr::set_colnames("Pearson")

# ligand expression Seurat dotplot
monotj$celltype <- monotj@active.ident
order_ligands_adapted <- gsub("[.]", "-", order_ligands_adapted)
order_ligands_adapted <- order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] <- "H2-M3"
order_ligands_adapted[order_ligands_adapted == "H2.T23"] <- "H2-T23"
rotated_dotplot <- DotPlot(
    monotj %>% subset(celltype %in% sender_celltypes),
    features = order_ligands_adapted, cols = "RdYlBu"
) +
    coord_flip() +
    theme(
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave("Ligand_Expr", rotated_dotplot)
layout <- "ABBB
           CCCC"
final <- p + p_ligand_target_network + p_ligand_receptor_network +
    plot_layout(design = layout) &
    theme(legend.position = "bottom")
ggsave("Merge.pdf", final, width = 15, height = 13)
