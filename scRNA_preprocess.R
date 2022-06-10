pkgs <- c(
    "Seurat", "SeuratWrappers", "ggplot2", "batchelor",
    "dplyr", "optparse", "reshape2", "data.table", "magrittr"
)

lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))

data_combined <- readRDS("deal.Rds")

## calculate the percentage of MT
data_combined[["percent.mt"]] <- PercentageFeatureSet(
    data_combined,
    pattern = "^MT-"
)

## Remove low quality cells, normalization and select top 2,000 variably expressed genes
data_combined %<>% subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20) %>%
    NormalizeData(normalization.method = "LogNormalize") %>%
    FindVariableFeatures(selection.method = "vst")

## Scale the data
all.genes <- rownames(data_combined)
data_combined %<>% ScaleData(features = all.genes) %>%
    ScaleData(vars.to.regress = "percent.mt")

## RunPCA for reducing the dimensionality
data_combined %<>% RunPCA(object = ., features = VariableFeatures(object = .)) %>%
    RunFastMNN(object = ., object.list = SplitObject(., split.by = "Patient"))

## batch effect correction and findclusters
data_combined %<>% FindNeighbors(reduction = "mnn", dims = 1:30) %>%
    FindClusters(resolution = 1.4)

data_combined %<>% RunTSNE(reduction = "mnn", dims = 1:30) %>%
    RunUMAP(reduction = "mnn", dims = 1:30)

## save RDS file for subsequent analysis
saveRDS(data_combined, "final_input.Rds")

pdf("umap.pdf", width = 9, height = 9)
DimPlot(data_combined, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

pdf("tsne.pdf", width = 9, height = 9)
DimPlot(data_combined, reduction = "tsne", label = TRUE) + NoLegend()
dev.off()

## Calculate the markers for each cluster
markers <- FindAllMarkers(data_combined, only.pos = TRUE, min.pct = 0.25)
write.table(markers, "umap_allmarker.txt", quote = F, sep = "\t")
