library(GSEABase)

test_that("list collections and sub-collections works", {
  msigdb = getMsigdb()
  
  #test empty
  expect_length(listCollections(GeneSetCollection(list())), 0)
  expect_length(listSubCollections(GeneSetCollection(list())), 0)
  
  #test non-broad sets
  nullgsc = GeneSetCollection(GeneSet(setName = 'A'))
  expect_error(listCollections(nullgsc))
  expect_error(listSubCollections(nullgsc))
  
  #test all
  expect_length(listCollections(msigdb), 9)
  expect_length(listSubCollections(msigdb), 18)
  
  #test Hallmarks
  hgsc = msigdb[grepl('HALLMARK', sapply(msigdb, setName))]
  expect_length(listCollections(hgsc), 1)
  expect_length(listSubCollections(hgsc), 0)
  
  #test GO
  gogsc = msigdb[grepl('^GO.._', sapply(msigdb, setName))]
  expect_length(listCollections(gogsc), 1)
  expect_length(listSubCollections(gogsc), 3)
})

test_that("subset collections works", {
  msigdb = getMsigdb()
  
  #test empty
  expect_error(subsetCollection(GeneSetCollection(list()), 'c1'))
  
  #test non-broad sets
  nullgsc = GeneSetCollection(GeneSet(setName = 'A'))
  expect_error(subsetCollection(nullgsc, 'c1'))
  
  expect_length(subsetCollection(msigdb, 'h'), 50)
  expect_length(subsetCollection(msigdb, 'c5'), sum(sapply(lapply(msigdb, collectionType), bcCategory) %in% 'c5'))
  expect_length(subsetCollection(msigdb, 'c5', 'GO:BP'), sum(sapply(lapply(msigdb, collectionType), bcSubCategory) %in% 'GO:BP'))
})
