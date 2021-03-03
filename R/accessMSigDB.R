#' Subset collections and sub-collections of MSigDB
#'
#' The molecular signatures database (MSigDB) is composed of collections and
#' sub-collection. Many analyses (e.g. gene-set enrichment using limma::fry) are
#' best carried out within specific collections rather than across the entire
#' database of signatures. This function allows subsetting of MSigDB data
#' objects within this package using collection and sub-collection types.
#'
#' @param gsc a GeneSetCollection object, storing GeneSets from the MSigDB
#' @param collection a character, stating the collection(s) to be retrieved. The
#'   collection(s) must be one from the [listCollections()] function.
#' @param subcollection a character, stating the sub-collection(s) to be
#'   retrieved. The sub-collection(s) must be one from the
#'   [listSubCollections()] function. If NULL, all sub-collections are
#'   retrieved.
#'
#' @return a GeneSetCollection object, containing gene sets belonging to the
#'   queries collection and/or sub-collection.
#' @export
#'
#' @examples
#' msigdb.hs.SYM <- msigdb.hs.SYM()
#' subsetCollection(msigdb.hs.SYM, collection = "h")
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
#' @inheritParams subsetCollection
#'
#' @return a character vector, containing character codes for all collections
#'   present in the GeneSetCollection object.
#' @export
#'
#' @examples
#' msigdb.hs.SYM <- msigdb.hs.SYM()
#' listCollections(msigdb.hs.SYM)
#' 
listCollections <- function(gsc) {
  cat = unique(sapply(lapply(gsc, GSEABase::collectionType), GSEABase::bcCategory))
  cat = as.character(na.omit(cat))
  return(cat)
}

#' List all sub-collection types within a MSigDB gene set collection
#'
#' This function lists all the sub-collection types present in a MSigDB gene set
#' collection. Descriptions of sub-collections can be found at the MSigDB
#' website.
#'
#' @inheritParams subsetCollection
#'
#' @return a character vector, containing character codes for all
#'   sub-collections present in the GeneSetCollection object.
#' @export
#'
#' @examples
#' msigdb.hs.SYM <- msigdb.hs.SYM()
#' listSubCollections(msigdb.hs.SYM)
#' 
listSubCollections <- function(gsc) {
  subcat = unique(sapply(lapply(gsc, GSEABase::collectionType), GSEABase::bcSubCategory))
  subcat = as.character(na.omit(subcat))
  return(subcat)
}

