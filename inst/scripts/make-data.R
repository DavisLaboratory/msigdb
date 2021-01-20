library(GSEABase)
library(BiocFileCache)

fpath = tempfile()
bfc = BiocFileCache(fpath, ask = FALSE)

getMsigdbData <- function(version) {
  msigdb_url = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/___/msigdb_v___.xml'
  msigdb_url = gsub('___', version, msigdb_url)
  
  #download file
  msigpath = bfcrpath(bfc, msigdb_url)
  
  #----Symbols----
  #read genesets into a geneset collection
  msigdb = getBroadSets(msigpath, membersId = 'MEMBERS_SYMBOLIZED')
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% 'CP:KEGG'])
  
  #gene set modifications
  msigdb = GeneSetCollection(lapply(msigdb, function(gs) {
    #add URLs
    urls(gs) = msigdb_url
    return(gs)
  }))
  
  #remove empty gene sets
  msigdb_lengths = sapply(lapply(msigdb, geneIds), length)
  msigdb = GeneSetCollection(msigdb[msigdb_lengths > 0])
  msigdb.sym = msigdb
  
  #----Entrez IDs----
  #read genesets into a geneset collection
  msigdb = getBroadSets(msigpath, membersId = 'MEMBERS_EZID')
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% 'CP:KEGG'])
  
  #gene set modifications
  msigdb = GeneSetCollection(lapply(msigdb, function(gs) {
    #convert idType
    gs@geneIdType = EntrezIdentifier()
    #add URLs
    urls(gs) = msigdb_url
    return(gs)
  }))
  
  #remove empty gene sets
  msigdb_lengths = sapply(lapply(msigdb, geneIds), length)
  msigdb = GeneSetCollection(msigdb[msigdb_lengths > 0])
  msigdb.ezid = msigdb
  
  return(list(msigdb.sym, msigdb.ezid))
}

v7.1 = getMsigdbData('7.1')
msigdb.hs.v7.1.SYM = v7.1[[1]]
msigdb.hs.v7.1.EZID = v7.1[[2]]
save(msigdb.hs.v7.1.SYM, file = 'msigdb.hs.v7.1.SYM.rda')
save(msigdb.hs.v7.1.EZID, file = 'msigdb.hs.v7.1.EZID.rda')

v7.2 = getMsigdbData('7.2')
msigdb.hs.v7.2.SYM = v7.2[[1]]
msigdb.hs.v7.2.EZID = v7.2[[2]]
save(msigdb.hs.v7.2.SYM, file = 'msigdb.hs.v7.2.SYM.rda')
save(msigdb.hs.v7.2.EZID, file = 'msigdb.hs.v7.2.EZID.rda')
