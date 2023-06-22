library(GSEABase)
library(BiocFileCache)
library(GO.db)
library(stringr)
library(reshape2)
library(limma)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(foreach)
library(doSNOW)
library(iterators)

#file cache to download files to
bfc = BiocFileCache()

#getBroadSets that can deal with duplicate geneset names
getBroadSets2 <- function (uri, ..., membersId = c("MEMBERS_SYMBOLIZED", "MEMBERS_EZID")) 
{
  membersId <- match.arg(membersId)
  factories <- sapply(uri, GSEABase:::.BroadXMLNodeToGeneSet_factory, 
                      membersId = membersId)
  tryCatch({
    geneSets <- unlist(mapply(GSEABase:::.fromXML, uri, "//GENESET", 
                              factories, SIMPLIFY = FALSE, USE.NAMES = FALSE))
  }, error = function(err) {
    stop("'getBroadSets' failed to create gene sets:\n  ", 
         conditionMessage(err), call. = FALSE)
  })
  
  snames = sapply(geneSets, setName)
  if (any(duplicated(snames))) {
    #identify dups
    dups = snames[duplicated(snames)]
    uq = geneSets[!snames %in% dups]
    #merge dups
    dupmerged = lapply(dups, function(x) {
      x = geneSets[snames %in% x]
      geneIds(x[[1]]) = do.call(union, lapply(x, geneIds))
      return(x[[1]])
    })
    geneSets = c(uq, dupmerged)
  }
  
  GeneSetCollection(geneSets)
}

#function to download any given version of MSigDB
getMsigdbData <- function(msigdb_ver, old = FALSE, is.zip = TRUE) {
  msigdb_url = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/___.Hs/msigdb_v___.Hs.xml'
  msigdb_url = gsub('___', msigdb_ver, msigdb_url)
  msigdb_url = ifelse(is.zip, paste0(msigdb_url, '.zip'), msigdb_url)
  
  #download file
  msigpath = bfcrpath(bfc, msigdb_url)
  msigpath = ifelse(is.zip, unzip(msigpath, exdir = tempdir()), msigpath)
  
  #replace C5:GO
  c5 = createC5OrgDb(org.Hs.egGO2ALLEGS, 'Homo sapiens', old)
  c5 = ez2symMsigdb(c5, org.Hs.eg.db)
  
  #----Symbols----
  #read genesets into a geneset collection
  msigdb = getBroadSets2(msigpath, membersId = 'MEMBERS_SYMBOLIZED')
  #remove KEGG (due to licenses)
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% 'CP:KEGG'])
  #remove archived
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcCategory) %in% 'archived'])
  #remove GO
  gosubcat = c('GO:BP', 'GO:MF', 'GO:CC')
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% gosubcat])
  
  #gene set modifications
  msigdb = GeneSetCollection(lapply(msigdb, function(gs) {
    #add URLs
    urls(gs) = msigdb_url
    setVersion(gs) = new('Versions', msigdb_ver)
    return(gs)
  }))
  
  #add GO
  msigdb = c(msigdb, c5$sym)
  
  #remove empty gene sets
  msigdb_lengths = sapply(lapply(msigdb, geneIds), length)
  msigdb = GeneSetCollection(msigdb[msigdb_lengths > 0])
  msigdb.sym = msigdb
  
  #----Entrez IDs----
  #read genesets into a geneset collection
  msigdb = getBroadSets2(msigpath, membersId = 'MEMBERS_EZID')
  #remove KEGG (due to licenses)
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% 'CP:KEGG'])
  #remove archived
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcCategory) %in% 'archived'])
  #remove GO
  gosubcat = c('GO:BP', 'GO:MF', 'GO:CC')
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% gosubcat])
  
  #gene set modifications
  msigdb = GeneSetCollection(lapply(msigdb, function(gs) {
    #convert idType
    gs@geneIdType = EntrezIdentifier()
    #add URLs
    urls(gs) = msigdb_url
    setVersion(gs) = new('Versions', msigdb_ver)
    return(gs)
  }))
  
  #add GO
  msigdb = c(msigdb, c5$ezid)
  
  #remove empty gene sets
  msigdb_lengths = sapply(lapply(msigdb, geneIds), length)
  msigdb = GeneSetCollection(msigdb[msigdb_lengths > 0])
  msigdb.ezid = msigdb
  
  return(list('sym' = msigdb.sym, 'ezid' = msigdb.ezid))
}

#function to convert the human MSigDB to mouse
createMmMsigdbData <- function(hsdb, old = FALSE) {
  #remove c1 and c5 genesets (these need to be replaced completely); except for HPO in c5
  rmgs = sapply(lapply(hsdb, collectionType), bcCategory) %in% c('c1', 'c5')
  rmgs[sapply(lapply(hsdb, collectionType), bcSubCategory) %in% 'HPO'] = FALSE
  mmdb = hsdb[!rmgs]
  
  #convert IDs using Ensembl homology annotations
  mmdb = lapply(mmdb, function(gs) {
    gids = geneIds(gs)
    gids = msigdb:::hcop$mouse_entrez_gene[msigdb:::hcop$human_entrez_gene %in% gids]
    gids = na.omit(unique(gids))
    geneIds(gs) = gids
    return(gs)
  })
  
  #create c5 category
  c5 = createC5OrgDb(org.Mm.egGO2ALLEGS, 'Mus musculus', old)
  
  #create c1 category
  c1 = createC1MmNCBI()
  
  #combine all and create Mm MSigDB
  mmdb = c(mmdb, c1, c5)
  mmdb = mmdb[order(sapply(mmdb, setName))]
  mmdb.all = ez2symMsigdb(mmdb, org.Mm.eg.db)
  
  return(mmdb.all)
}

ez2symMsigdb <- function(sigdb, orgdb) {
  #only retain genes with ids in the OrgDb
  allg = unique(keys(orgdb, 'ENTREZID'))
  sigdb = lapply(sigdb, function(gs) {
    geneIds(gs) = intersect(geneIds(gs), allg)
    return(gs)
  })
  
  #remove empty genesets
  sigdb = sigdb[sapply(lapply(sigdb, geneIds), length) > 0]
  sigdb = GeneSetCollection(sigdb)
  
  #convert to symbols
  gids = unique(unlist(geneIds(sigdb)))
  gmap = mapIds(orgdb, keys = gids, column = 'SYMBOL', keytype = 'ENTREZID')
  sigdb.sym = lapply(sigdb, function(gs) {
    geneIds(gs) = na.omit(unique(gmap[geneIds(gs)]))
    gs@geneIdType = SymbolIdentifier()
    return(gs)
  })
  sigdb.sym = GeneSetCollection(sigdb.sym)
  
  return(list('sym' = sigdb.sym, 'ezid' = sigdb))
}

#code adapted from Gordon K. Smyth and Alexandra Garnham (WEHI)
createC1MmNCBI <- function() {
  ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz'
  
  #download file
  geneinfo_path = bfcrpath(bfc, ncbi_url)
  x = read.delim(geneinfo_path)
  
  # Extract cytoband string
  Loc = strsplit2(x$map_location, split = '\\|')
  i1 = grep(' [A-Z]', Loc[, 1])
  i2 = grep(' [A-Z]', Loc[, 2])
  Cytoband = rep('', nrow(Loc))
  Cytoband[i1] = Loc[i1, 1]
  Cytoband[i2] = Loc[i2, 2]
  
  # Remove empty strings
  i = Cytoband != ''
  GeneID = as.character(x$GeneID[i])
  Cytoband = Cytoband[i]
  
  # Split the cytoband string by space and '-'
  CB2 = strsplit2(Cytoband, split = '-')
  CB1 = strsplit2(CB2[, 1], split = ' ')
  Cytoband1 = paste(CB1[, 1], CB1[, 2], sep = 'q')
  Tab1 = data.frame(GeneID = GeneID, Cytoband = Cytoband1)
  
  # Genes in >1 cytoband
  i = grep('-', Cytoband)
  Cytoband2 = paste(CB1[i, 1], str_trim(CB2[i, 2]), sep = 'q')
  Tab2 = data.frame(GeneID = GeneID[i], Cytoband = Cytoband2)
  
  # Special case
  Tab3 = data.frame(GeneID = c('170942', '654820'), Cytoband = c('YqE', 'YqE'))
  
  # Make sets
  Tab = rbind(Tab1, Tab2, Tab3)
  
  #create geneset
  posgenes = split(Tab, Tab$Cytoband)
  c1 = lapply(posgenes, function (x) {
    gs = GeneSet(
      as.character(unique(x$GeneID)),
      setName = x$Cytoband[1],
      collectionType = BroadCollection(category = 'c1'),
      shortDescription = paste('Ensembl Genes in Cytogenetic Band', x$Cytoband[1]),
      organism = 'Mus musculus'
    )
    geneIdType(gs) = EntrezIdentifier()
    return(gs)
  })
  
  return(c1)
}

createC5OrgDb <- function(db, org, old) {
  gomap = as.list(db)
  gomap = lapply(gomap, as.character)
  
  #create c5 genesets using OrgDb
  c5 = mapply(
    function(genes, gsname, gstype, gsdesc) {
      go_prefix = ifelse(old, 'GO_', paste0('GO', gstype, '_'))
      gs = GeneSet(
        as.character(unique(genes)),
        setName = paste0(go_prefix, gsub(' ', '_', str_to_upper(gsname))),
        collectionType = BroadCollection(category = 'c5', subCategory = paste0('GO:', gstype)),
        shortDescription = gsdesc,
        organism = org
      )
      geneIdType(gs) = EntrezIdentifier()
      return(gs)
    },
    gomap,
    mapIds(GO.db, names(gomap), column = 'TERM', keytype = 'GOID'),
    mapIds(GO.db, names(gomap), column = 'ONTOLOGY', keytype = 'GOID'),
    mapIds(GO.db, names(gomap), column = 'DEFINITION', keytype = 'GOID'),
    SIMPLIFY = FALSE
  )
  
  return(c5)
}

createHCOPmap <- function() {
  hcop_url = 'https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz'
  hcop_path = bfcrpath(bfc, hcop_url)
  
  #read mappings
  hcop = read.delim(hcop_path)
  hcop = unique(hcop[, c('mouse_entrez_gene', 'human_entrez_gene')])
  hcop = hcop[!hcop$mouse_entrez_gene %in% '-' & !hcop$human_entrez_gene %in% '-', ]
  
  return(hcop)
}

computeIdf <- function(msigGsc) {
  rmwords = vissE:::getMsigBlacklist()
  
  signames = sapply(msigGsc, GSEABase::setName)
  sigdesc_s = sapply(msigGsc, GSEABase::description)
  docs = list('Name' = signames, 'Short' = sigdesc_s)
  docs = lapply(docs, unique)
  
  #text-mining
  docs = lapply(docs, function(d) tm::Corpus(tm::VectorSource(d)))
  toSpace <- tm::content_transformer(function (x, pattern) gsub(pattern, " ", x))
  docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "_"))
  docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "/"))
  docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "@"))
  docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "\\|"))
  docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "\\("))
  docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "\\)"))
  
  # Convert the text to lower case
  docs = lapply(docs, function(d) tm::tm_map(d, tm::content_transformer(tolower)))
  # Remove numbers
  # docs = lapply(docs, function(d) tm_map(d, removeNumbers))
  # Remove english common stopwords
  docs = lapply(docs, function(d) tm::tm_map(d, tm::removeWords, tm::stopwords('english')))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs = lapply(docs, function(d) tm::tm_map(d, tm::removeWords, rmwords))
  # Remove punctuations
  docs = lapply(docs, function(d) tm::tm_map(d, tm::removePunctuation))
  # Eliminate extra white spaces
  docs = lapply(docs, function(d) tm::tm_map(d, tm::stripWhitespace))
  # Remove full numbers
  docs = lapply(docs, function(d) tm::tm_filter(d, function(x) !grepl('\\b[0-9]+\\b', x)))
  # Text lemmatisation
  docs = lapply(docs, function(d) tm::tm_map(d, textstem::lemmatize_strings))
  
  #compute idf
  dtms = lapply(docs, tm::TermDocumentMatrix)
  dtms = lapply(dtms, as.matrix)
  
  #compute IDF
  idfs = lapply(dtms, function(x) {
    idf = log(ncol(x) / rowSums(x != 0))
    return(idf)
  })
  #sort names to quicken searches
  idfs = lapply(idfs, function(x) {
    x[sort(names(x))]
  })
  
  return(idfs)
}

processMsigdbData <- function(msigdb_ver, old = FALSE, is.zip = TRUE) {
  bname = paste0('msigdb.v', msigdb_ver)
  
  #download data
  msigdb = getMsigdbData(msigdb_ver, old, is.zip)
  
  #create human data
  msigdb.hs.SYM = msigdb[[1]]
  msigdb.hs.EZID = msigdb[[2]]
  saveRDS(msigdb.hs.SYM, file = paste0(bname, '.hs.SYM.rds'))
  saveRDS(msigdb.hs.EZID, file = paste0(bname, '.hs.EZID.rds'))
  
  #create mouse data
  msigdb.mm = createMmMsigdbData(msigdb.hs.EZID, old)
  msigdb.mm.SYM = msigdb.mm[[1]]
  msigdb.mm.EZID = msigdb.mm[[2]]
  saveRDS(msigdb.mm.SYM, file = paste0(bname, '.mm.SYM.rds'))
  saveRDS(msigdb.mm.EZID, file = paste0(bname, '.mm.EZID.rds'))
  
  #compute IDF
  idf_hs = computeIdf(msigdb.hs.EZID)
  idf_mm = computeIdf(msigdb.mm.EZID)
  saveRDS(idf_hs, file = paste0(bname, '.hs.idf.rds'))
  saveRDS(idf_mm, file = paste0(bname, '.mm.idf.rds'))
}

#----human-mouse homologs----
hcop = createHCOPmap()
usethis::use_data(hcop, internal = TRUE, overwrite = TRUE)

#----prepare MSigDB for the package----
cl = makeSOCKcluster(10, outfile = '')
registerDoSNOW(cl)

processMsigdbData('2022.1', is.zip = FALSE)
processMsigdbData('2023.1')

stopCluster(cl)

