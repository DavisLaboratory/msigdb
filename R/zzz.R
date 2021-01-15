#' The Molecular Signatures Database (MSigDB)
#'
#' This ExperimentHub package contains gene expression signatures from the
#' molecular signatures database (MSigDB). The molecular signatures database
#' (MSigDB) is a collection of over 25000 gene expression signatures that are
#' grouped into collections and sub-collections. Metadata associated with
#' signatures is collected and stored in the data in this package.
#'
#' All data in this package are stored in a GeneSetCollection objects from the
#' GSEABase package. Each gene expression signature in the collection is stored
#' in a GeneSet object from the GSEABase package.
#'
#' The following datasets are included in this package:
#'
#' 1. msigdb.v7.1.SYM - The MSigDB v7.1 with gene expression signatures defined
#' using gene symbols.
#'
#' 2. msigdb.v7.1.EZID - The MSigDB v7.1 with gene expression signatures defined
#' using Entrez IDs.
#'
#' @format A GeneSetCollection object composed of GeneSet objects representing
#'   all non-empty gene expression signatures from the molecular signatures
#'   database (MSigDB).
#' @references Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert,
#'   B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment
#'   analysis: a knowledge-based approach for interpreting genome-wide
#'   expression profiles. Proceedings of the National Academy of Sciences,
#'   102(43), 15545-15550.
#'
#'   Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo,
#'   P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0.
#'   Bioinformatics, 27(12), 1739-1740.
#'
#'   Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P.,
#'   & Tamayo, P. (2015). The molecular signatures database hallmark gene set
#'   collection. Cell systems, 1(6), 417-425.
#'   
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' msigdb_datasets <- query(eh, "msigdb")
#'
#' @name msigdb
#' @aliases msigdb-package
#' 
NULL

.onLoad <- function(libname, pkgname) {
  fl = system.file("extdata", "metadata.csv", package = pkgname)
  titles = utils::read.csv(fl, stringsAsFactors = FALSE)$Title
  
  ExperimentHub::createHubAccessors(pkgname, 'msigdb.v7.1.SYM')
  ExperimentHub::createHubAccessors(pkgname, 'msigdb.v7.1.EZID')
}

