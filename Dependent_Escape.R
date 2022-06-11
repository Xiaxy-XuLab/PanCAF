# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Hallmarker analysis
#'
#' This function is applied to calculate the hallmarker value for each cell
#'

RunEscape <- function(object,
                      lib = "H") {
    gshallmark <- getGeneSets(library = lib)
    esseurat <- enrichIt(
        obj = object, gene.sets = gshallmark,
        groups = 1000, cores = 30
    )
    write.table(esseurat, "Hall_Marker.txt", quote = F, sep = "\t") # nolint
    object <- Seurat::AddMetaData(object, esseurat)
    return(object)
}
