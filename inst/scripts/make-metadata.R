#params
org_name = c('hs' = 'Homo sapiens', 'mm' = 'Mus musculus')
org_taxid = c('hs' = 9606, 'mm' = 10090)

#build description for each file
buildDescription <- function(ver, org, type) {
  if (type %in% 'SYM') {
    sprintf('Gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets. Signatures are represented using gene symbols.', org_name[org], ver)
  } else if (type %in% 'EZID') {
    sprintf('Gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets. Signatures are represented using Entrez IDs.', org_name[org], ver)
  } else if (type %in% 'adj') {
    sprintf('This data is for internal use within vissE. Adjacency matrix representation of gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets. Signatures are represented using Entrez IDs.', org_name[org], ver)
  } else {
    sprintf('This data is for internal use within vissE. Inverse document frequencies of words from gene expression signatures (%s) from the Molecular Signatures Database (v%s) excluding KEGG gene sets.', org_name[org], ver)
  }
}

#retrieve file list that need to be updated
allfiles = list.files(pattern = '.rds$')
old_meta = read.csv('inst/extdata/metadata.csv')
allfiles = setdiff(allfiles, basename(old_meta$RDataPath))

if (length(allfiles) > 0) {
  meta = plyr::ldply(allfiles, function(x) {
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
      RDataPath = file.path('msigdb', x)
    )
  })
  
  meta = rbind(meta, old_meta)
  write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)
}

#test metadata file
ExperimentHubData::makeExperimentHubMetadata('../msigdb')
