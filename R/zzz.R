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
#' in a GeneSet object from the GSEABase package. This data does not include
#' KEGG gene sets due to copyrights. Users can download this data using
#' functions provided in the package (see Details).
#'
#' The following datasets are included in this package:
#'
#' 1. msigdb.hs.SYM - The MSigDB v7.2 for human with gene expression signatures
#' defined using gene symbols.
#'
#' 2. msigdb.hs.EZID - The MSigDB v7.2 for human with gene expression signatures
#' defined using Entrez IDs.
#'
#' 3. msigdb.mm.SYM - The MSigDB v7.2 for mouse with gene expression signatures
#' defined using gene symbols.
#'
#' 4. msigdb.mm.EZID - The MSigDB v7.2 for mouse with gene expression signatures
#' defined using Entrez IDs.
#'
#' @format A GeneSetCollection object composed of GeneSet objects representing
#'   all non-empty gene expression signatures from the molecular signatures
#'   database (MSigDB).
#' @details Data in this package does not include gene sets from the KEGG
#'   database due to licensing limitations. Users can use the [appendKEGG()]
#'   function in this package to download KEGG gene sets directly from the
#'   MSigDB and append to existing data objects.
#'
#'   The mouse MSigDB is created by translating human genes to mouse homologs
#'   using annotations from the Mouse Genome Informatics (MGI) database for most
#'   gene sets. Gene sets in the collections c1 (positional gene sets) and c5
#'   (ontologies) are recreated as information in these gene sets is organism
#'   specific. Positional gene sets are created using Ensembl 102 annotations
#'   from biomaRt. Gene sets representing gene ontologies are derived from the
#'   mouse R/Bioconductor organism database (org.Mm.eg.db).
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
#' @section Acknowledgement: MSigDB is protected by copyright © 2004-2020 Broad
#'   Institute, Inc., Massachusetts Institute of Technology, and Regents of the
#'   University of California. Use of MSigDB is subject to the terms and
#'   conditions of the Creative Commons Attribution 4.0 International License -
#'   <https://creativecommons.org/licenses/by/4.0/>.
#'
#'   MSigDB gene sets derived from BioCarta pathways are the subject of
#'   copyright © 2000-2017 BioCarta, and are subject to Biocarta's Disclaimer of
#'   Liability and of Warranties -
#'   <https://data.broadinstitute.org/gsea-msigdb/msigdb/biocarta/biocarta_disclaimer_of_liability_and_of_warranties.txt>
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
  
  ExperimentHub::createHubAccessors(pkgname, 'msigdb.hs.SYM')
  ExperimentHub::createHubAccessors(pkgname, 'msigdb.hs.EZID')
  ExperimentHub::createHubAccessors(pkgname, 'msigdb.mm.SYM')
  ExperimentHub::createHubAccessors(pkgname, 'msigdb.mm.EZID')
}
