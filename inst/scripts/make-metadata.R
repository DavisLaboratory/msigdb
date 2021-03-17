meta = data.frame(
  Title = c(
    'msigdb.v7.2.hs.SYM',
    'msigdb.v7.2.hs.EZID',
    'msigdb.v7.2.mm.SYM',
    'msigdb.v7.2.mm.EZID'
  ),
  Description = c(
    'Gene expression signatures (human) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using gene symbols.',
    'Gene expression signatures (human) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using Entrez IDs.',
    'Gene expression signatures (mouse) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using gene symbols.',
    'Gene expression signatures (mouse) from the Molecular Signatures Database (v7.2) excluding KEGG gene sets. Signatures are represented using Entrez IDs.'
  ),
  BiocVersion = 3.13,
  Genome = NA,
  SourceType = c('XML'),
  SourceUrl = rep(
    c(
      'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml'
    ),
    each = 4
  ),
  SourceVersion = rep(c('7.2'), each = 4),
  Species = rep(c('Homo sapiens', 'Homo sapiens', 'Mus musculus', 'Mus musculus'), 1),
  TaxonomyId = rep(c(9606, 9606, 10090, 10090), 1),
  Coordinate_1_based = TRUE,
  DataProvider = 'Broad Institute',
  Maintainer = 'Dharmesh D. Bhuva <bhuva.d@wehi.edu.au>',
  RDataClass = 'GSEABase::GeneSetCollection',
  DispatchClass = 'Rda',
  RDataPath = c(
    'msigdb/msigdb.v7.2.hs.SYM.rda',
    'msigdb/msigdb.v7.2.hs.EZID.rda',
    'msigdb/msigdb.v7.2.mm.SYM.rda',
    'msigdb/msigdb.v7.2.mm.EZID.rda'
  )
)

write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata('../msigdb')
