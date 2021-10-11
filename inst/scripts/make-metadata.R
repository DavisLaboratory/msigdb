#retrieve file list that need to be updated
allfiles = list.files(pattern = '.rds$', recursive = TRUE)
allfiles = setdiff(allfiles, basename(old_meta$RDataPath))
old_meta = read.csv('inst/extdata/metadata.csv')
msigdb_files = basename(grep('msigdb', allfiles, value = TRUE))
imex_files = basename(grep('imex', allfiles, value = TRUE))
meta = c()

#----msigdb metadata----
#params
org_name = c('hs' = 'Homo sapiens', 'mm' = 'Mus musculus')
org_taxid = c('hs' = 9606, 'mm' = 10090)

#build description for each file
buildDescription <- function(ver, org, type) {
  if (type %in% 'SYM') {
    sprintf('Gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets. Signatures are represented using gene symbols.', org_name[org], ver)
  } else if (type %in% 'EZID') {
    sprintf('Gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets. Signatures are represented using Entrez IDs.', org_name[org], ver)
  } else {
    sprintf('This data is for internal use within vissE. Inverse document frequencies of words from gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets.', org_name[org], ver)
  }
}

if (length(msigdb_files) > 0) {
  msigdb_meta = plyr::ldply(msigdb_files, function(x) {
    #create regex
    fregex = 'msigdb.v([0-9]\\.[0-9])\\.([a-z]+)\\.(.*).rds$'
    
    #extract file info
    ver = gsub(fregex, '\\1', x)
    org = gsub(fregex, '\\2', x)
    type = gsub(fregex, '\\3', x)
    
    #build record
    c(
      Title = gsub('.rds$', '', x),
      Description = buildDescription(ver, org, type),
      BiocVersion = as.numeric(as.character(BiocManager::version())) + 0.01,
      Genome = NA,
      SourceType = 'XML',
      SourceUrl = sprintf('https://data.broadinstitute.org/gsea-msigdb/msigdb/release/%s/msigdb_v%s.xml', ver, ver),
      SourceVersion = ver,
      Species = as.character(org_name[org]),
      TaxonomyId = as.numeric(org_taxid[org]),
      Coordinate_1_based = TRUE,
      DataProvider = 'Broad Institute',
      Maintainer = 'Dharmesh D. Bhuva <bhuva.d@wehi.edu.au>',
      RDataClass = 'GSEABase::GeneSetCollection',
      DispatchClass = 'Rds',
      RDataPath = file.path('msigdb', paste0('msigdb.v', ver), x)
    )
  })
  
  meta = rbind(meta, msigdb_meta)
}

#----IMEX metadata----
##build record
meta = rbind(meta, c(
  Title = 'imex_hsmm_0721',
  Description = 'Protein-protein interaction (PPI) network for human and mouse obtained from the international molecular exchange (IMEX) in July 2021.',
  BiocVersion = as.numeric(as.character(BiocManager::version())) + 0.01,
  Genome = NA,
  SourceType = 'TXT',
  SourceUrl = 'ftp://ftp.ebi.ac.uk/pub/databases/intact/2021-07-06/psimitab/intact-micluster.txt',
  SourceVersion = '2021-07-06',
  Species = 'Homo sapiens',
  TaxonomyId = '9606',
  Coordinate_1_based = TRUE,
  DataProvider = 'EBI',
  Maintainer = 'Dharmesh D. Bhuva <bhuva.d@wehi.edu.au>',
  RDataClass = 'data.frame',
  DispatchClass = 'Rds',
  RDataPath = file.path('msigdb', 'IMEx', 'imex_hsmm_0721.rds')
))

#write and test metadata file
meta = rbind(old_meta, meta)
write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)
ExperimentHubData::makeExperimentHubMetadata('../msigdb')
