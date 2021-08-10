library(GSEABase)
library(BiocFileCache)
library(biomaRt)
library(org.Mm.eg.db)
library(GO.db)
library(stringr)
library(reshape2)
library(limma)

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
  #remove KEGG (due to licenses)
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcSubCategory) %in% 'CP:KEGG'])
  #remove archived
  msigdb = GeneSetCollection(msigdb[!sapply(lapply(msigdb, collectionType), bcCategory) %in% 'archived'])
  
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
  c5 = createC5MmOrgDb(old)
  
  #create c1 category
  c1 = createC1MmNCBI()
  
  ####
  ## What should we do about the Human Phenotype Ontology?
  ##  1. Replace with mammalian phenotype ontology
  ##  2. Translate to mouse
  ####

  #combine all and create Mm MSigDB
  mmdb = c(mmdb, c1, c5)
  mmdb = mmdb[order(sapply(mmdb, setName))]
  
  #only retain genes with ids in the OrgDb
  allg = unique(keys(org.Mm.eg.db, 'ENTREZID'))
  mmdb = lapply(mmdb, function(gs) {
    geneIds(gs) = intersect(geneIds(gs), allg)
    return(gs)
  })
  
  #remove empty genesets
  mmdb = mmdb[sapply(lapply(mmdb, geneIds), length) > 0]
  mmdb = GeneSetCollection(mmdb)
  
  #convert to symbols
  gids = unique(unlist(geneIds(mmdb)))
  gmap = mapIds(org.Mm.eg.db, keys = gids, column = 'SYMBOL', keytype = 'ENTREZID')
  mmdb.sym = lapply(mmdb, function(gs) {
    geneIds(gs) = na.omit(unique(gmap[geneIds(gs)]))
    gs@geneIdType = SymbolIdentifier()
    return(gs)
  })
  mmdb.sym = GeneSetCollection(mmdb.sym)
  
  return(list('sym' = mmdb.sym, 'ezid' = mmdb))
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

createC5MmOrgDb <- function(old) {
  gomap = as.list(org.Mm.egGO2ALLEGS)
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
        organism = 'Mus musculus'
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

#----human-mouse homologs----
hcop = createHCOPmap()
usethis::use_data(hcop, internal = TRUE, overwrite = TRUE)

#----Version 7.2----
#create human data
msigdb = getMsigdbData('7.2')
msigdb.v7.2.hs.SYM = msigdb[[1]]
msigdb.v7.2.hs.EZID = msigdb[[2]]
save(msigdb.v7.2.hs.SYM, file = 'msigdb.v7.2.hs.SYM.rda')
save(msigdb.v7.2.hs.EZID, file = 'msigdb.v7.2.hs.EZID.rda')

#create mouse data
msigdb.mm = createMmMsigdbData(msigdb.v7.2.hs.EZID, old = TRUE)
msigdb.v7.2.mm.SYM = msigdb.mm[[1]]
msigdb.v7.2.mm.EZID = msigdb.mm[[2]]
save(msigdb.v7.2.mm.SYM, file = 'msigdb.v7.2.mm.SYM.rda')
save(msigdb.v7.2.mm.EZID, file = 'msigdb.v7.2.mm.EZID.rda')

#----Version 7.4----
#create human data
msigdb = getMsigdbData('7.4')
msigdb.v7.4.hs.SYM = msigdb[[1]]
msigdb.v7.4.hs.EZID = msigdb[[2]]
save(msigdb.v7.4.hs.SYM, file = 'msigdb.v7.4.hs.SYM.rda')
save(msigdb.v7.4.hs.EZID, file = 'msigdb.v7.4.hs.EZID.rda')

#create mouse data
msigdb.mm = createMmMsigdbData(msigdb.v7.4.hs.EZID)
msigdb.v7.4.mm.SYM = msigdb.mm[[1]]
msigdb.v7.4.mm.EZID = msigdb.mm[[2]]
save(msigdb.v7.4.mm.SYM, file = 'msigdb.v7.4.mm.SYM.rda')
save(msigdb.v7.4.mm.EZID, file = 'msigdb.v7.4.mm.EZID.rda')



