#' Add KEGG pathway derived gene sets to a collection
#'
#' This function adds gene sets derived from KEGG pathways to the MSigDB data
#' stored in this package. Direct user-end download from the MSigDB is required
#' to ensure KEGG licenses are adhered to.
#'
#' @inheritParams subsetCollection
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
  idType = getMsigIdType(gsc)
  org = getMsigOrganism(gsc, idType)
  id = ifelse(is(idType, 'SymbolIdentifier') & org %in% 'hs', 'symbols', 'entrez')
  
  #create URL
  fname = paste0('c2.cp.kegg.v', version, '.', id, '.gmt')
  link = paste('https://data.broadinstitute.org/gsea-msigdb/msigdb/release', version, fname, sep = '/')
  
  #download and process KEGG genesets
  gsc_kegg = GSEABase::getGmt(
    link,
    geneIdType = idType,
    collectionType = GSEABase::BroadCollection('c2', 'CP:KEGG')
  )
  
  if (is(idType, 'EntrezIdentifier')) {
    gsc_kegg = lapply(gsc_kegg, function(gs) {
      gs@geneIdType = idType
      return(gs)
    })
    gsc_kegg = GSEABase::GeneSetCollection(gsc_kegg)
  }
  
  #conversions for mouse
  if (org %in% 'mm') {
    #convert Hs to Mm (Entrez IDs)
    gsc_kegg = lapply(gsc_kegg, function(gs) {
      gids = GSEABase::geneIds(gs)
      gids = hcop$mouse_entrez_gene[hcop$human_entrez_gene %in% gids]
      gids = stats::na.omit(unique(gids))
      GSEABase::geneIds(gs) = gids
      return(gs)
    })
    
    #convert to symbols if needed
    allg = unique(unlist(lapply(gsc_kegg, GSEABase::geneIds)))
    if (methods::is(idType, 'SymbolIdentifier')) {
      gmap = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                   keys = allg,
                                   column = 'SYMBOL',
                                   keytype = 'ENTREZID')
      gsc_kegg = lapply(gsc_kegg, function(gs) {
        GSEABase::geneIds(gs) = stats::na.omit(unique(gmap[GSEABase::geneIds(gs)]))
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


#' Infer organism type for the gene set collection
#'
#' Since both Human and Mouse MSigDB collections are hosted in this package,
#' this function infers the type of organism represented in a gene set
#' collection based on the gene IDs present. If not all gene IDs belong to the
#' same organism, the organism with more than 50% gene IDs present in the
#' collection is returned. In any other case, the function returns an error.
#'
#' @param idType a GSEABase::SymbolIdentifier or GSEABASE::EntrezIdentifier
#'   object, representing the ID type inferred from the [getMsigIdType()]
#'   function. Avoid providing this manually.
#' @inheritParams getMsigIdType
#'
#' @return a character, either "mm" (representing Mus musculus - mouse) or "hs"
#'   (representing Homo sapiens - human).
#' @export
#'
#' @examples
#' msigdb.v7.2.hs.SYM <- msigdb.v7.2.hs.SYM()
#' id <- getMsigIdType(msigdb.v7.2.hs.SYM)
#' getMsigOrganism(msigdb.v7.2.hs.SYM(), id)
#' 
getMsigOrganism <- function(gsc, idType) {
  #ensure ID types are the same in the collection
  keytype = ifelse(methods::is(idType, 'SymbolIdentifier'), 'SYMBOL', 'ENTREZID')
  
  #check gene IDs against organism databases
  allg = unique(unlist(GSEABase::geneIds(gsc)))
  
  if (all(allg %in% AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype))) {
    return(c('hs'))
  } else if (mean(allg %in% AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype)) > 0.5) {
    warning('Assuming the organism to be human.')
    return(c('hs'))
  } else if (all(allg %in% AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype))) {
    return(c('mm'))
  } else if (mean(allg %in% AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype)) > 0.5) {
    warning('Assuming the organism to be mouse.')
    return(c('mm'))
  } else {
    stop('Cannot infer organism.')
  }
}

#' Infer gene identifier type for the gene set collection
#'
#' The gene identifier (Symbol or Entrez ID) of a gene set collection is
#' inferred from the IDs present in the data. A collection should ideally store
#' gene sets using a single identifier type. This function returns the
#' identifier type (either SymbolIdentifier or EntrezIdentifier) of the
#' collection. It returns an error if the identifier is neither of these.
#'
#' @inheritParams subsetCollection
#'
#' @return a GSEABase::SymbolIdentifier or GSEABASE::EntrezIdentifier object,
#'   specifying the gene identifier type (gene symbols or Entrez IDs
#'   respectively).
#' @export
#'
#' @examples
#' msigdb.v7.2.hs.SYM <- msigdb.v7.2.hs.SYM()
#' id <- getMsigIdType(msigdb.v7.2.hs.SYM)
#' 
getMsigIdType <- function(gsc) {
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
