library(GSEABase)
library(BiocFileCache)
library(biomaRt)
library(org.Mm.eg.db)
library(GO.db)
library(stringr)

#file cache to download files to
fpath = tempfile()
bfc = BiocFileCache(fpath, ask = FALSE)

#function to download any given version of MSigDB
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

#function to convert the human MSigDB to mouse
# hsdb needs to be the db formed of symbols, not EZIDs
createMmMsigdbData <- function(hsdb, isSym = TRUE) {
  mousechr = c(as.character(1:19), 'MT', 'X', 'Y')
  
  #remove c1 and c5 genesets (these need to be replaced completely)
  mmdb = hsdb[!sapply(lapply(hsdb, collectionType), bcCategory) %in% c('c1', 'c5')]
  
  #convert IDs using MGI
  allg = unique(unlist(lapply(mmdb, geneIds)))
  allg = convertMouseGeneList(allg)
  gmap = allg$MGI.symbol
  names(gmap) = allg$HGNC.symbol
  
  mmdb = lapply(mmdb, function(gs) {
    gids = na.omit(unique(gmap[geneIds(gs)]))
    geneIds(gs) = gids
    return(gs)
  })
  
  #create c5 category
  gomap = as.list(org.Mm.egGO2ALLEGS)
  gomap = lapply(gomap, as.character)
  
  if (isSym) {
    #convert to symbols
    idmap = mapIds(org.Mm.eg.db, unique(unlist(gomap)), 'SYMBOL', keytype = 'ENTREZID')
    gomap = lapply(gomap, function (x) {
      x = as.character(na.omit(idmap[x]))
      return(x)
    })
  }
  
  #create c5 genesets
  c5 = mapply(
    function(genes, gsname, gstype, gsdesc) {
      gs = GeneSet(
        unique(genes),
        setName = paste0('GO_', gsub(' ', '_', str_to_upper(gsname))),
        collectionType = BroadCollection(category = 'c5', subCategory = paste0('GO:', gstype)),
        shortDescription = gsdesc,
        organism = 'Mus musculus'
      )
      if (isSym) {
        geneIdType(gs) = SymbolIdentifier()
      } else {
        geneIdType(gs) = EntrezIdentifier()
      }
      return(gs)
    },
    gomap,
    mapIds(GO.db, names(gomap), column = 'TERM', keytype = 'GOID'),
    mapIds(GO.db, names(gomap), column = 'ONTOLOGY', keytype = 'GOID'),
    mapIds(GO.db, names(gomap), column = 'DEFINITION', keytype = 'GOID'),
    SIMPLIFY = FALSE
  )
  
  #create c1 category
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  idtype = ifelse(isSym, 'mgi_symbol', 'entrezgene_id')
  posgenes = select(
    mouse,
    keys = mousechr,
    columns = c('chromosome_name', 'band', idtype),
    keytype = 'chromosome_name'
  )
  posgenes$band = gsub('\\..*', '', posgenes$band)
  posgenes = posgenes[posgenes[, 3] != '' & !is.na(posgenes[, 3]), ]
  posgenes$band = paste(posgenes$chromosome_name, posgenes$band, sep = 'q')
  posgenes$band = paste0('chr', posgenes$band)
  posgenes = split(posgenes, posgenes$band)
  c1 = lapply(posgenes, function (x) {
    gs = GeneSet(
      unique(x[, 3]),
      setName = x$band[1],
      collectionType = BroadCollection(category = 'c1'),
      shortDescription = paste('Ensembl Genes in Cytogenetic Band', x$band[1]),
      organism = 'Mus musculus'
    )
    if (isSym) {
      geneIdType(gs) = SymbolIdentifier()
    } else {
      geneIdType(gs) = EntrezIdentifier()
    }
    return(gs)
  })
  
  #combine all and create Mm MSigDB
  mmdb = c(mmdb, c1, c5)
  mmdb = mmdb[order(sapply(mmdb, setName))]
  #only retain genes on the primary scaffold
  allg = select(
    mouse,
    keys = mousechr,
    columns = ifelse(isSym, 'mgi_symbol', 'entrezgene_id'),
    keytype = 'chromosome_name'
  )[, 1]
  mmdb = lapply(mmdb, function(gs) {
    geneIds(gs) = intersect(geneIds(gs), allg)
    return(gs)
  })
  #remove empty genesets
  mmdb = mmdb[sapply(lapply(mmdb, geneIds), length) > 0]
  
  ####
  ## What should we do about the Human Phenotype Ontology?
  ##  1. Replace with mammalian phenotype ontology
  ##  2. Translate to mouse
  ####
  
  mmdb = GeneSetCollection(mmdb)
  
  return(mmdb)
}

## https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
# Basic function to convert human to mouse gene names
convertMouseGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  return(genesV2)
}

#----create human data----
msigdb = getMsigdbData('7.2')
msigdb.hs.SYM = msigdb[[1]]
msigdb.hs.EZID = msigdb[[2]]
save(msigdb.hs.SYM, file = 'msigdb.hs.SYM.rda')
save(msigdb.hs.EZID, file = 'msigdb.hs.EZID.rda')

#----create mouse data----
msigdb.mm.SYM = createMmMsigdbData(msigdb.hs.SYM)
msigdb.mm.EZID = createMmMsigdbData(msigdb.hs.SYM, isSym = FALSE)



