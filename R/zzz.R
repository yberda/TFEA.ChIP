.onAttach <- function(libname, pkgname)
{
    packageStartupMessage( paste0(
        "Due to space limitations, the internal database in TFEA.ChIP includes only ",
        "ChIP-seq experiments from cell types classified as ENCODE tiers 1, 2, and 2.5.\n",
        "To download the full ReMap2022 database, as well as other ready-to-use databases, ",
        "visit https://github.com/yberda/ChIPDBData and https://github.com/LauraPS1/TFEA.ChIP_downloads"
    )
    )
}

utils::globalVariables(c("chip_metadata", "MetaData", "tf.name", "chip.name", "EntrezID", "Accession", "TF", "pval", "pval.adj", "OR", "pval.adj.safe", "rank_score", "log10.adj.pVal", "log2.OR", "highlight", "pointText", "Arg.ES", "ES", "symbol", "Treatment", "Cell", "x", "y", "RES", "arg.ES"))
