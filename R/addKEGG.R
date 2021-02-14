#' @importFrom  org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
NULL

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
  
  #get GeneSetCollection information
  idType = getIdType(gsc)
  org = getOrganism(gsc, idType)
  id = ifelse(class(idType) %in% 'SymbolIdentifier' & org %in% 'hs', 'symbols', 'entrez')
  
  #create URL
  fname = paste0('c2.cp.kegg.v', version, '.', id, '.gmt')
  link = paste('https://data.broadinstitute.org/gsea-msigdb/msigdb/release', version, fname, sep = '/')
  
  #download and process KEGG genesets
  gsc_kegg = GSEABase::getGmt(
    link,
    geneIdType = idType,
    collectionType = GSEABase::BroadCollection('c2', 'CP:KEGG')
  )
  
  if (class(idType) %in% 'EntrezIdentifier') {
    gsc_kegg = lapply(gsc_kegg, function(gs) {
      gs@geneIdType = idType
      return(gs)
    })
    gsc_kegg = GSEABase::GeneSetCollection(gsc_kegg)
  }
  
  #conversions for mouse
  if (org %in% 'mm') {
    #convert Hs to Mm (Entrez IDs)
    allg = unique(unlist(GSEABase::geneIds(gsc_kegg)))
    gmap = convertMouseGeneList(allg)
    
    gsc_kegg = lapply(gsc_kegg, function(gs) {
      gids = na.omit(unique(gmap[GSEABase::geneIds(gs)]))
      GSEABase::geneIds(gs) = gids
      return(gs)
    })
    
    #convert to symbols if needed
    if (class(idType) %in% 'SymbolIdentifier') {
      gmap = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = gmap,
                                   column = 'SYMBOL',
                                   keytype = 'ENTREZID')
      gsc_kegg = lapply(gsc_kegg, function(gs) {
        GSEABase::geneIds(gs) = na.omit(unique(gmap[GSEABase::geneIds(gs)]))
        gs@geneIdType = GSEABase::SymbolIdentifier()
        return(gs)
      })
    }
  }
  
  #remove empty gene sets
  gsc_kegg = GSEABase::GeneSetCollection(gsc_kegg[sapply(lapply(gsc_kegg, GSEABase::geneIds), length) > 0])
  
  #annotate gene set with metadata
  gsc_kegg = GSEABase::GeneSetCollection(lapply(gsc_kegg, function(gs) {
    #add URLs
    GSEABase::urls(gs) = link
    return(gs)
  }))
  
  gsc = GSEABase::GeneSetCollection(c(gsc, gsc_kegg))
  
  return(gsc)
}

getOrganism <- function(gsc, idType) {
  #ensure ID types are the same in the collection
  keytype = ifelse(class(idType) %in% 'SymbolIdentifier', 'SYMBOL', 'ENTREZID')
  
  #check gene IDs against organism databases
  allg = unique(unlist(GSEABase::geneIds(gsc)))
  
  if (all(allg %in% AnnotationDbi::keys(org.Hs.eg.db, keytype))) {
    return(c('hs'))
  } else if (mean(allg %in% AnnotationDbi::keys(org.Hs.eg.db, keytype)) > 0.5) {
    warning('Assuming the organism to be human.')
    return(c('hs'))
  } else if (all(allg %in% AnnotationDbi::keys(org.Mm.eg.db, keytype))) {
    return(c('mm'))
  } else if (mean(allg %in% AnnotationDbi::keys(org.Mm.eg.db, keytype)) > 0.5) {
    warning('Assuming the organism to be mouse.')
    return(c('mm'))
  } else {
    stop('Cannot infer organism.')
  }
}

getIdType <- function(gsc) {
  idType = sapply(gsc, function(gs) class(GSEABase::geneIdType(gs)))
  idType = unique(idType)
  if (length(idType) != 1) {
    stop('Gene identifiers should be the same in a collection.')
  }
  
  #create GSEABase ID type
  if (idType %in% 'SymbolIdentifier') {
    idType = GSEABase::SymbolIdentifier()
  } else if (idType %in% 'EntrezIdentifier'){
    idType = GSEABase::EntrezIdentifier()
  } else{
    stop('Only Symbol and Entrez identifiers supported.')
  }
  
  return(idType)
}

convertMouseGeneList <- function(genes){
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  #convert symbols using homology
  mgi_hgnc = biomaRt::getLDS(
    attributes = c("entrezgene_id"),
    filters = "entrezgene_id",
    values = genes,
    mart = human,
    attributesL = c("entrezgene_id"),
    martL = mouse,
    uniqueRows = TRUE
  )
  
  #create a named vector to simplify mapping
  gmap = as.character(mgi_hgnc[, 2])
  names(gmap) = as.character(mgi_hgnc[, 1])
  
  return(gmap)
}
