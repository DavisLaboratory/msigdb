library(GSEABase)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

test_that("appendKEGG works", {
  gs1 = GeneSet(setName = 'gs1', geneIdType = SymbolIdentifier())
  gs2 = GeneSet(setName = 'gs2', geneIdType = SymbolIdentifier())
  gsc = GeneSetCollection(c(gs1, gs2))
  expect_length(appendKEGG(gsc), 188)
})

test_that("ID type inference works", {
  gs1 = GeneSet(setName = 'gs1', geneIdType = SymbolIdentifier())
  gs1 = GeneSet(setName = 'gs1', geneIdType = SymbolIdentifier())
  gs2 = GeneSet(setName = 'gs2', geneIdType = EntrezIdentifier())
  gsc = GeneSetCollection(c(gs1, gs2))
  
  expect_s4_class(getIdType(GeneSetCollection(gs1)), 'SymbolIdentifier')
  expect_s4_class(getIdType(GeneSetCollection(gs2)), 'EntrezIdentifier')
  expect_error(getIdType(gsc), 'identifiers should be the same')
  expect_error(getIdType(gs1))
  
  gs1 = GeneSet(setName = 'gs1', geneIdType = NullIdentifier())
  gs2 = GeneSet(setName = 'gs2', geneIdType = NullIdentifier())
  gsc = GeneSetCollection(c(gs1, gs2))
  
  expect_error(getIdType(gsc), 'identifiers supported')
})

test_that("Organism inference works", {
  gsmm = GeneSet(c('Esr1'), setName = 'gsmm', geneIdType = SymbolIdentifier())
  gshs = GeneSet(c('ESR1'), setName = 'gshs', geneIdType = SymbolIdentifier())
  gsmm_e = gsmm
  geneIdType(gsmm_e) = EntrezIdentifier('org.Mm.eg.db')
  setName(gsmm_e) = 'gsmm_e'
  gshs_e = gshs
  setName(gshs_e) = 'gshs_e'
  geneIdType(gshs_e) = EntrezIdentifier('org.Hs.eg.db')
  
  expect_equal(getOrganism(GeneSetCollection(gshs), getIdType(GeneSetCollection(gshs))), 'hs')
  expect_equal(getOrganism(GeneSetCollection(gshs_e), getIdType(GeneSetCollection(gshs_e))), 'hs')
  expect_equal(getOrganism(GeneSetCollection(gsmm), getIdType(GeneSetCollection(gsmm))), 'mm')
  expect_equal(getOrganism(GeneSetCollection(gsmm_e), getIdType(GeneSetCollection(gsmm_e))), 'mm')
  
  expect_error(getOrganism(GeneSetCollection(gshs, gsmm), getIdType(GeneSetCollection(gsmm))))
})
