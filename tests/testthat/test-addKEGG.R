library(GSEABase)

test_that("appendKEGG works", {
  gs1 = GeneSet(setName = 'gs1', geneIdType = SymbolIdentifier())
  gs2 = GeneSet(setName = 'gs2', geneIdType = EntrezIdentifier())
  gsc = GeneSetCollection(c(gs1, gs2))
  
  expect_error(appendKEGG(gsc), 'identifiers should be the same')
  
  gs1 = GeneSet(setName = 'gs1', geneIdType = NullIdentifier())
  gs2 = GeneSet(setName = 'gs2', geneIdType = NullIdentifier())
  gsc = GeneSetCollection(c(gs1, gs2))
  
  expect_error(appendKEGG(gsc), 'identifiers supported')
  
  gs1 = GeneSet(setName = 'gs1', geneIdType = SymbolIdentifier())
  gs2 = GeneSet(setName = 'gs2', geneIdType = SymbolIdentifier())
  gsc = GeneSetCollection(c(gs1, gs2))
  expect_length(appendKEGG(gsc), 188)
})
