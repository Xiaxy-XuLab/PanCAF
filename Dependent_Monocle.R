# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Trajactory analysis
#'
#' This function is used to construct the possible evolutionary
#' trajectory of the incorporated cells based on their
#' @param object the object
#' 

RunMonocle <- function(object) {
    monomatrix <- as(
        as.matrix(GetAssayData(object, slot = "counts")), "sparseMatrix"
    )
    feature_ann <- data.frame(
        gene_id = rownames(monomatrix), gene_short_name = rownames(monomatrix)
    )
    rownames(feature_ann) <- rownames(monomatrix)
    monofd <- new("AnnotatedDataFrame", data = feature_ann)
    sample_ann <- object@meta.data
    monopd <- new("AnnotatedDataFrame", data = sample_ann)

    monocds <- newCellDataSet(monomatrix,
        phenoData = monopd,
        featureData = monofd,
        lowerDetectionLimit = 0.1,
        expressionFamily = negbinomial.size()
    )

    head(pData(monocds))
    head(fData(monocds))

    monocds <- estimateSizeFactors(monocds)
    monocds <- estimateDispersions(monocds)

    monocds <- detectGenes(monocds, min_expr = 0.1)
    print(head(fData(monocds)))
    expressed_genes <- row.names(subset(fData(monocds), num_cells_expressed >= 50)) # nolint
    monocds <- monocds[expressed_genes, ]

    disp_table <- dispersionTable(monocds)
    unsup_clustering_genes <- subset(
        disp_table, mean_expression >= 0.05 &
            dispersion_empirical >= 2 * dispersion_fit
    ) #
    monocds <- setOrderingFilter(monocds, unsup_clustering_genes$gene_id)

    monocds <- reduceDimension(
        monocds,
        max_components = 2,
        method = "DDRTree"
    )
    monocds <- orderCells(monocds)
    return(monocds)
}
