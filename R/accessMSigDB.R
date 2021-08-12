#' Retrieve MSigDB data hosted on the hub
#'
#' Download molecular signatures database (MSigDB) hosted on the ExperimentHub
#' or retrieve pre-downloaded version from cache. This package currently hosts
#' versions greater than 7.2 for human and mouse with both symbol and Entrez
#' identifiers.
#'
#' @param org a character, representing the organism whose signature database
#'   needs to be retrieved ("hs" for human and "mm" for mouse).
#' @param id a character, representing the ID type to use ("SYM" for gene
#'   symbols and "EZID" for Entrez IDs).
#' @param version a character, stating the version of MSigDB to be retrieved
#'   (should be >= 7.2). See `getMsigdbVersions()`.
#'
#' @return a GeneSetCollection, containing GeneSet objects from the specified
#'   version of the molecular signatures database (MSigDB).
#' @export
#'
#' @examples
#' gsc = getMsigdb('hs', 'SYM')
#' 
getMsigdb <- function(org = c('hs', 'mm'), id = c('SYM', 'EZID'), version = getMsigdbVersions()) {
  org = match.arg(org)
  id = match.arg(id)
  version = match.arg(version)
  checkMsigdbVersion(version)
  
  #create object name
  obj_name = paste0('msigdb.v', version, '.', org, '.', id)
  gsc = getMSigdbObject(obj_name)
  
  return(gsc)
}

getMSigdbObject <- function(obj_name) {
  #load object
  eh = ExperimentHub::ExperimentHub()
  info = AnnotationHub::mcols(AnnotationHub::query(eh, 'msigdb'))
  id = rownames(info)[info$title %in% obj_name]
  if (length(id) != 1)
    stop('Data not found')
  
  return(suppressWarnings(eh[[id]]))
}

#' Subset collections and sub-collections of MSigDB
#'
#' The molecular signatures database (MSigDB) is composed of collections and
#' sub-collection. Many analyses (e.g. gene-set enrichment using limma::fry) are
#' best carried out within specific collections rather than across the entire
#' database of signatures. This function allows subsetting of MSigDB data
#' objects within this package using collection and sub-collection types.
#'
#' @param collection a character, stating the collection(s) to be retrieved. The
#'   collection(s) must be one from the [listCollections()] function.
#' @param subcollection a character, stating the sub-collection(s) to be
#'   retrieved. The sub-collection(s) must be one from the
#'   [listSubCollections()] function. If NULL, all sub-collections are
#'   retrieved.
#' 
#' @inheritParams getMsigdb
#'
#' @return a GeneSetCollection object, containing gene sets belonging to the
#'   queries collection and/or sub-collection.
#' @export
#'
#' @examples
#' gsc = getMsigdb('hs', 'SYM')
#' subsetCollection(gsc, collection = "h")
#' 
subsetCollection <- function(gsc, collection, subcollection = NULL) {
  stopifnot(length(gsc) > 0)
  stopifnot(collection %in% listCollections(gsc))
  stopifnot(is.null(subcollection) | subcollection %in% listSubCollections(gsc))
  
  if (is.null(subcollection))
    subcollection = c(listSubCollections(gsc), NA)
  
  #filter collection & sub-collection
  ctype = lapply(gsc, GSEABase::collectionType)
  gsc = gsc[sapply(ctype, GSEABase::bcCategory) %in% collection &
              sapply(ctype, GSEABase::bcSubCategory) %in% subcollection]
  
  return(gsc)
}

#' List all collection types within a MSigDB gene set collection
#'
#' This function lists all the collection types present in a MSigDB gene set
#' collection. Descriptions of collections can be found at the MSigDB website.
#'
#' @inheritParams getMsigdb
#'
#' @return a character vector, containing character codes for all collections
#'   present in the GeneSetCollection object.
#' @export
#'
#' @examples
#' gsc = getMsigdb('hs', 'SYM')
#' listCollections(gsc)
#' 
listCollections <- function(gsc) {
  cat = unique(sapply(lapply(gsc, GSEABase::collectionType), GSEABase::bcCategory))
  cat = as.character(stats::na.omit(cat))
  return(cat)
}

#' List all sub-collection types within a MSigDB gene set collection
#'
#' This function lists all the sub-collection types present in a MSigDB gene set
#' collection. Descriptions of sub-collections can be found at the MSigDB
#' website.
#'
#' @inheritParams getMsigdb
#'
#' @return a character vector, containing character codes for all
#'   sub-collections present in the GeneSetCollection object.
#' @export
#'
#' @examples
#' gsc = getMsigdb('hs', 'SYM')
#' listSubCollections(gsc)
#' 
listSubCollections <- function(gsc) {
  subcat = unique(sapply(lapply(gsc, GSEABase::collectionType), GSEABase::bcSubCategory))
  subcat = as.character(stats::na.omit(subcat))
  return(subcat)
}

checkMsigdbVersion <- function(version) {
  stopifnot(version %in% getMsigdbVersions())
}