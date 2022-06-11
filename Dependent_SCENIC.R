# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' SCENIC analysis
#'
#' This function is used to probe the transcription factor.
#' @param exp the raw counts of the incorporated object
#'
#'

RunSCENIC <- function(exp) {
    org <- "hgnc"
    dbs <- defaultDbNames[[org]]
    scenicoptions <- initializeScenic(
        org = org, dbDir = "/work/xiaxy/work/cisTarget_databases",
        dbs = dbs, nCores = 60
    )

    Genes_Keep <- geneFiltering(exp, # nolint
        scenicOptions = scenicoptions,
        minCountsPerGene = 3 * .01 * ncol(exp),
        minSamples = ncol(exp) * .01
    )
    exprmat_filtered <- exp[Genes_Keep, ]

    runCorrelation(exprmat_filtered, scenicoptions)
    exprmat_filtered <- log2(exprmat_filtered + 1)
    runGenie3(exprmat_filtered, scenicoptions)

    ### Build and score the GRN
    exprmat_log <- log2(exp + 1)
    scenicoptions <- runSCENIC_1_coexNetwork2modules(scenicoptions)
    scenicoptions <- runSCENIC_2_createRegulons(
        scenicoptions,
        coexMethod = c("top5perTarget")
    ) #** Only for toy run!!
    scenicoptions@settings$nCores <- 1
    scenicoptions <- runSCENIC_3_scoreCells(scenicoptions, exprmat_log)

    saveRDS(
        scenicoptions,
        file = "morescenicoptions.Rds"
    )

    regulons <- loadInt(scenicoptions, "regulons")
    regulons <- loadInt(scenicoptions, "aucell_regulons")
    regulontargetsinfo <- loadInt(scenicoptions, "regulonTargetsInfo")

    regulonauc <- loadInt(scenicoptions, "aucell_regulonAUC")
    regulonauc <- regulonauc[onlyNonDuplicatedExtended(rownames(regulonauc)), ]
    write.table(
        regulonauc@assays@data$AUC, "SCENIC.txt",
        quote = F, sep = "\t"
    )
}
