#' Add KEGG pathway derived gene sets to a collection
#'
#' This function adds gene sets derived from KEGG pathways to the MSigDB data
#' stored in this package. Direct user-end download from the MSigDB is required
#' to ensure KEGG licenses are adhered to.
#'
#' @param gsc a GeneSetCollection object, storing GeneSets from the MSigDB
#' @param id a character, either 'sym' or 'ezid' representing the gene
#'   identifier to download (symbol or entrez id respectively).
#' @param version a character, stating the version of msigdb to download
#'   (>=7.1).
#'
#' @return a GeneSetCollection object, storing gene sets from the MSigDB
#'   including the downloaded KEGG gene sets.
#' @export
#'
#' @examples
#' library(msigdb)
#' library(GSEABase)
#'
#' gsc <- msigdb.hs.v7.1.SYM()
#' gsc <- appendKEGG(gsc, id = 'sym', version = '7.1')
#' gsc
#' 
appendKEGG <- function(gsc, id = c('sym', 'ezid'), version = '7.1') {
  id = match.arg(id)
  stopifnot(version %in% c('7.1'))
  idtype = NULL
  
  #create file path
  if (id %in% 'sym') {
    id = 'symbols'
    idtype = GSEABase::SymbolIdentifier()
  } else{
    id = 'entrez'
    idtype = GSEABase::EntrezIdentifier()
  }
  fname = paste0('c2.cp.kegg.v', version, '.', id, '.gmt')
  link = paste('https://data.broadinstitute.org/gsea-msigdb/msigdb/release/', version, fname, sep = '/')
  
  #download and process KEGG genesets
  gsc_kegg = getGmt(link,
                    geneIdType = idtype,
                    collectionType = BroadCollection('c2', 'CP:KEGG'))
  gsc = GSEABase::GeneSetCollection(c(gsc, gsc_kegg))
  
  return(gsc)
}