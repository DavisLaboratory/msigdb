library(GSEABase)
library(BiocFileCache)
library(biomaRt)
library(org.Mm.eg.db)
library(GO.db)
library(stringr)
library(reshape2)

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
  
  return(list('sym' = msigdb.sym, 'ezid' = msigdb.ezid))
}

#function to convert the human MSigDB to mouse
createMmMsigdbData <- function(hsdb) {
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
  c5 = createC5MmOrgDb()
  
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

createC1MmbiomaRt <- function() {
  mousechr = c(as.character(1:19), 'MT', 'X', 'Y')
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  #retrieve band information
  posgenes = select(
    mouse,
    keys = mousechr,
    columns = c('chromosome_name', 'band', 'entrezgene_id'),
    keytype = 'chromosome_name'
  )
  #discard sub-banding (decimals)
  posgenes$band = gsub('\\..*', '', posgenes$band)
  #discard annotations where gene id missing
  posgenes = posgenes[posgenes[, 3] != '' & !is.na(posgenes[, 3]), ]
  #combine chr and band to create locus
  posgenes$band = paste(posgenes$chromosome_name, posgenes$band, sep = 'q')
  posgenes$band = paste0('chr', posgenes$band)
  #create geneset
  posgenes = split(posgenes, posgenes$band)
  c1 = lapply(posgenes, function (x) {
    gs = GeneSet(
      as.character(unique(x[, 3])),
      setName = x$band[1],
      collectionType = BroadCollection(category = 'c1'),
      shortDescription = paste('Ensembl Genes in Cytogenetic Band', x$band[1]),
      organism = 'Mus musculus'
    )
    geneIdType(gs) = EntrezIdentifier()
    return(gs)
  })
  
  return(c1)
}

createC5MmOrgDb <- function() {
  gomap = as.list(org.Mm.egGO2ALLEGS)
  gomap = lapply(gomap, as.character)
  
  #create c5 genesets using OrgDb
  c5 = mapply(
    function(genes, gsname, gstype, gsdesc) {
      gs = GeneSet(
        as.character(unique(genes)),
        setName = paste0('GO_', gsub(' ', '_', str_to_upper(gsname))),
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

#----create human data----
msigdb = getMsigdbData('7.2')
msigdb.hs.SYM = msigdb[[1]]
msigdb.hs.EZID = msigdb[[2]]
save(msigdb.hs.SYM, file = 'msigdb.v7.2.hs.SYM.rda')
save(msigdb.hs.EZID, file = 'msigdb.v7.2.hs.EZID.rda')

#----create mouse data----
hcop = createHCOPmap()
usethis::use_data(hcop, internal = TRUE, overwrite = TRUE)

msigdb.mm = createMmMsigdbData(msigdb.hs.EZID)
msigdb.mm.SYM = msigdb.mm[[1]]
msigdb.mm.EZID = msigdb.mm[[2]]
save(msigdb.mm.SYM, file = 'msigdb.v7.2.mm.SYM.rda')
save(msigdb.mm.EZID, file = 'msigdb.v7.2.mm.EZID.rda')



