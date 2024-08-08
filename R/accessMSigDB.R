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
  
  #create object name
  obj_name = paste0('msigdb.v', version, '.', org, '.', id)
  gsc = getMSigdbObject(obj_name)
  
  return(gsc)
}

#' Retrieve IMEx PPI hosted on the hub
#'
#' Download International Molecular Exchange (IMEx) protein-protein interaction
#' (PPI) hosted on the ExperimentHub or retrieve pre-downloaded version from
#' cache. This package currently hosts versions for human and mouse with both
#' symbol and Entrez identifiers.
#'
#' @param org a character, representing the organism whose PPI database needs to
#'   be retrieved ("hs" for human and "mm" for mouse).
#' @param inferred a logical, indicating whether inference from other organisms
#'   should be included in the PPI.
#' @param version a character, stating the version of IMEX to be retrieved. See
#'   `getMsigdbVersions()`.
#'
#' @return a data.frame, containing the IMEx PPI.
#' @export
#'
#' @examples
#' imex = getIMEX("hs")
#'
getIMEX <- function(org = c("hs", "mm"), inferred = FALSE, version = getIMEXVersions()) {
  org = match.arg(org)
  version = match.arg(version)
  org = c("hs" = "9606", "mm" = "10090")[org]

  # create object name
  version = as.Date(version)
  obj_name = paste0("imex_hsmm_", format(version, "%m"), format(version, "%y"))
  imex = getMSigdbObject(obj_name)
  imex = imex[imex$Taxid %in% org & (imex$Inferred | inferred), ]
  imex = as.data.frame(imex)

  return(imex)
}

#' Retrieve MSigDB data hosted on the hub
#'
#' Download molecular signatures database (MSigDB) hosted on the ExperimentHub
#' or retrieve pre-downloaded version from cache. This package currently hosts
#' versions greater than 7.2 for human and mouse with both symbol and Entrez
#' identifiers.
#'
#' @param org a character, representing the organism whose signature database
#'   needs to be retrieved ("hs" for human and "mm" for mouse).
#' @param version a character, stating the version of MSigDB to be retrieved
#'   (should be >= 7.2). See `getMsigdbVersions()`.
#'
#' @return a list of named numeric vectors, containing inverse document frequency (IDF) weights. Names represent terms that the IDF is computed for. IDFs are computed using gene-set names ("Name") and short descriptions ("Short").
#' @export
#'
#' @examples
#' gsc = getMsigdbIDF("hs")
#'
getMsigdbIDF <- function(org = c('hs', 'mm'), version = getMsigdbVersions()) {
  org = match.arg(org)
  version = match.arg(version)
  
  #create object name
  obj_name = paste0('msigdb.v', version, '.', org, '.idf')
  idf = getMSigdbObject(obj_name)
  
  return(idf)
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
#' @param gsc a GeneSetCollection object, containing MSigDB genesets in the form
#'   of GeneSet objects.
#' @param collection a character, stating the collection(s) to be retrieved. The
#'   collection(s) must be one from the [listCollections()] function.
#' @param subcollection a character, stating the sub-collection(s) to be
#'   retrieved. The sub-collection(s) must be one from the
#'   [listSubCollections()] function.
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
subsetCollection <- function(gsc, collection = c(), subcollection = c()) {
  stopifnot(length(gsc) > 0)
  stopifnot(all(collection %in% listCollections(gsc)))
  stopifnot(all(subcollection %in% listSubCollections(gsc)))
  
  #filter collection & sub-collection
  ctype = lapply(gsc, GSEABase::collectionType)
  if (length(subcollection) == 0){
    gsc = gsc[sapply(ctype, GSEABase::bcCategory) %in% collection]
  } else {
    gsc = gsc[sapply(ctype, GSEABase::bcCategory) %in% collection & 
                sapply(ctype, GSEABase::bcSubCategory) %in% subcollection]
  }
  
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
