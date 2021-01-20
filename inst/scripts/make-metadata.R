meta = data.frame(
  Title = c('msigdb.hs.v7.1.SYM', 'msigdb.hs.v7.1.EZID', 'msigdb.hs.v7.2.SYM', 'msigdb.hs.v7.2.EZID'),
  Description = c(
    'Gene expression signatures from the Molecular Signatures Database (v7.1) excluding KEGG gene sets. Signatures are represented using gene symbols.',
    'Gene expression signatures from the Molecular Signatures Database (v7.1) excluding KEGG gene sets. Signatures are represented using Entrez IDs.',
    'Gene expression signatures from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using gene symbols.',
    'Gene expression signatures from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using Entrez IDs.'
  ),
  BiocVersion = c(3.13, 3.13, 3.13, 3.13),
  Genome = NA,
  SourceType = c('XML'),
  SourceUrl = c(
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/msigdb_v7.1.xml',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/msigdb_v7.1.xml',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml'
  ),
  SourceVersion = c('7.1', '7.1', '7.2', '7.2'),
  Species = c('Homo sapiens', 'Homo sapiens', 'Homo sapiens', 'Homo sapiens'),
  TaxonomyId = c(9606, 9606, 9606, 9606),
  Coordinate_1_based = TRUE,
  DataProvider = c('Broad Institute'),
  Maintainer = 'Dharmesh D. Bhuva <bhuva.d@wehi.edu.au>',
  RDataClass = 'GSEABase::GeneSetCollection',
  DispatchClass = 'Rda',
  RDataPath = c(
    'msigdb/msigdb.hs.v7.1.SYM.rda',
    'msigdb/msigdb.hs.v7.1.EZID.rda',
    'msigdb/msigdb.hs.v7.2.SYM.rda',
    'msigdb/msigdb.hs.v7.2.EZID.rda'
  )
)

write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata('../msigdb')
