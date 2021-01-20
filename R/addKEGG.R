#' Add KEGG pathway derived gene sets to a collection
#'
#' This function adds gene sets derived from KEGG pathways to the MSigDB data
#' stored in this package. Direct user-end download from the MSigDB is required
#' to ensure KEGG licenses are adhered to.
#'
#' @param gsc a GeneSetCollection object, storing GeneSets from the MSigDB
#' @param id a character, either 'sym' or 'ezid' representing the gene
#'   identifier to download (symbol or entrez id respectively).
#'
#' @return a GeneSetCollection object, storing gene sets from the MSigDB
#'   including the downloaded KEGG gene sets.
#' @export
#'
#' @examples
#' library(GSEABase)
#'
#' gs1 <- GeneSet(setName = 'gs1', geneIdType = SymbolIdentifier())
#' gsc <- GeneSetCollection(gs1)
#' gsc <- appendKEGG(gsc)
#' gsc
#' 
appendKEGG <- function(gsc) {
  version = '7.2'
  
  #get IdType
  idType = sapply(gsc, function(gs) class(GSEABase::geneIdType(gs))[1])
  idType = unique(idType)
  if (length(idType) != 1) {
    stop('Gene identifiers should be the same in a collection.')
  }
  
  #create file path
  if (idType %in% 'SymbolIdentifier') {
    id = 'symbols'
    idtype = GSEABase::SymbolIdentifier()
  } else if (idType %in% 'EntrezIdentifier'){
    id = 'entrez'
    idtype = GSEABase::EntrezIdentifier()
  } else{
    stop('Only Symbol and Entrez identifiers supported.')
  }
  
  #create URL
  fname = paste0('c2.cp.kegg.v', version, '.', id, '.gmt')
  link = paste('https://data.broadinstitute.org/gsea-msigdb/msigdb/release', version, fname, sep = '/')
  
  #download and process KEGG genesets
  gsc_kegg = GSEABase::getGmt(
    link,
    geneIdType = idtype,
    collectionType = GSEABase::BroadCollection('c2', 'CP:KEGG')
  )
  
  #remove empty gene sets
  gsc_kegg = GSEABase::GeneSetCollection(gsc_kegg[sapply(lapply(gsc_kegg, GSEABase::geneIds), length) > 0])
  
  #gene set modifications
  gsc_kegg = GSEABase::GeneSetCollection(lapply(gsc_kegg, function(gs) {
    #add URLs
    GSEABase::urls(gs) = link
    return(gs)
  }))
  
  gsc = GeneSetCollection(c(gsc, gsc_kegg))
  
  return(gsc)
}