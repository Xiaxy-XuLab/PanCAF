# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Read files for Seurat
#'
#' This function will read in the expression file and integrate
#' the seurat object. In addition, it can optionally detect
#' and reject doublet using the appropriate software.
#' @param inputdir A path containing the cellranger output files or matrix files
#' @param outputdir Output path of the merged RDS file
#' @param type Data format type (cellranger output files or matrix files)
#' @param mincells Include features detected in at least this many cells.
#' @param deledoublet Whether to test and reject doublet
#' @param doubletpercent Percentage of excluded doublet
#' @param doubletpc Number of PC quadrants included
#' @param doubletresolution resoluttion value for doublet
#' @param wid width of picture
#' @param hei height of picture
#'

read_data <- function(inputdir,
                            outputdir,
                            type,
                            mincells = 3,
                            deledoublet = TRUE,
                            doubletpercent = 0.075,
                            doubletpc = seq(30),
                            doubletresolution = 1,
                            wid = 6,
                            hei = 6) {
    if (is.null(inputdir)) stop("inputdir is null")
    setwd(inputdir)
    files <- list.files()
    info <- data.table(filepath = files, sampleid = files)
    samplelist <- list()

    mat <- matrix(NA, nrow = nrow(info), ncol = 3)

    # read files
    for (i in seq(nrow(info))) {
        dir_of_10x <- info$filepath[i]
        sampleid <- info$sampleid[i]

        print(paste("Starting processing", sampleid, "at", Sys.time()))

        if (type == "cellranger") {
            sc_mat <- Read10X(dir_of_10x)
            sc_tumor <- CreateSeuratObject(
                counts = sc_mat, project = sampleid, min.cells = mincells
            )
        } else {
            sc_mat <- read.table(
                info$filepath[i],
                header = T, sep = "\t", row.names = 1
            )
            sc_mat <- data.frame(sc_mat)
            sc_mat <- as(as.matrix(sc_mat), "dgCMatrix")
            sc_tumor <- CreateSeuratObject(
                counts = sc_mat, project = sampleid, min.cells = mincells
            )
        }

        mat[i, 1] <- sampleid
        mat[i, 2] <- nrow(sc_tumor@meta.data)

        # rename cellnames
        sc_tumor <- RenameCells(sc_tumor, add.cell.id = sampleid)
        sc_tumor <- RenameCells(
            sc_tumor,
            new.names = gsub("-", "_", rownames(sc_tumor@meta.data))
        )

        if (deledoublet) {
            sc_tumor <- NormalizeData(
                sc_tumor,
                normalization.method = "LogNormalize",
                scale.factor = 10000
            )
            sc_tumor <- FindVariableFeatures(sc_tumor,
                selection.method = "vst", nfeatures = 2000
            )
            sc_tumor <- ScaleData(sc_tumor)
            sc_tumor <- RunPCA(sc_tumor)
            sc_tumor <- FindNeighbors(sc_tumor, dims = doubletpc)
            sc_tumor <- FindClusters(sc_tumor, resolution = doubletresolution)
            sc_tumor$seurat_clusters <- sc_tumor@active.ident

            sweep.res <- paramSweep_v3(sc_tumor, PCs = doubletpc, sct = T)
            sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
            bcmvn <- find.pK(sweep.stats)
            mpk <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

            annotations <- sc_tumor@meta.data$seurat_clusters
            homotypic <- modelHomotypic(annotations)
            nexp_poi <- round(doubletpercent * ncol(sc_tumor@assays$RNA@data))
            nexp_poiadj <- round(nexp_poi * (1 - homotypic))

            filterdouble <- doubletFinder_v3(
                sc_tumor,
                PCs = doubletpc,
                pN = 0.25, pK = mpk, nExp = nexp_poi,
                reuse.pANN = FALSE, sct = T
            )
            print(table(filterdouble@meta.data[, 7]))

            sc_tumor <- subset(
                sc_tumor,
                cells = rownames(filterdouble@meta.data[
                    which(filterdouble@meta.data[, 7] == "Singlet"),
                ])
            )
        }
        mat[i, 3] <- nrow(sc_tumor@meta.data)
        sc_tumor$sampleid <- sampleid
        samplelist <- append(samplelist, sc_tumor)
        print(paste("Finishing processing", sampleid, "at", Sys.time()))
    }

    combined <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), samplelist)

    saveRDS(combined, paste0(outputdir, "/read.Rds"))
    mat <- as.data.frame(mat)
    colnames(mat) <- c("sampleid", "raw", "doublet")
    write.table(
        mat,
        paste0(outputdir, "/cellnumber.txt"),
        quote = F, sep = "\t"
    )
    p1 <- ggplot(mat, aes(x = sampleid, y = raw)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        labs(y = "Cell number") +
        geom_text(
            aes(label = raw, vjust = 0.8, hjust = 0.5)
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none"
        )
    p2 <- ggplot(mat, aes(x = sampleid, y = doublet)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        labs(y = "Cell number") +
        geom_text(
            aes(label = raw, vjust = 0.8, hjust = 0.5)
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none"
        )
    final <- p1 + p2 + plot_layout(ncol = 1)
    ggsave(
        paste0(outputdir, "/raw_doublet.pdf"),
        final,
        width = wid, height = hei
    )
    saveRDS(final, paste0(outputdir, "/p1.Rds"))
    return(combined)
}
