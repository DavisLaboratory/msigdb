meta = data.frame(
  Title = c('msigdb.hs.SYM', 'msigdb.hs.EZID'),
  Description = c(
    'Gene expression signatures (human) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using gene symbols.',
    'Gene expression signatures (human) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using Entrez IDs.',
    'Gene expression signatures (mouse) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using gene symbols.',
    'Gene expression signatures (mouse) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using Entrez IDs.'
  ),
  BiocVersion = c(3.13, 3.13, 3.13, 3.13),
  Genome = NA,
  SourceType = c('XML'),
  SourceUrl = c(
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml'
  ),
  SourceVersion = c('7.2', '7.2', '7.2', '7.2'),
  Species = c('Homo sapiens', 'Homo sapiens', 'Mus musculus', 'Mus musculus'),
  TaxonomyId = c(9606, 9606, 10090, 10090),
  Coordinate_1_based = TRUE,
  DataProvider = c('Broad Institute'),
  Maintainer = 'Dharmesh D. Bhuva <bhuva.d@wehi.edu.au>',
  RDataClass = 'GSEABase::GeneSetCollection',
  DispatchClass = 'Rda',
  RDataPath = c(
    'msigdb/msigdb.hs.SYM.rda',
    'msigdb/msigdb.hs.EZID.rda',
    'msigdb/msigdb.mm.SYM.rda',
    'msigdb/msigdb.mm.EZID.rda'
  )
)

write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata('../msigdb')
