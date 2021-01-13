library(GSEABase)
library(BiocFileCache)

msigdb_url = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/msigdb_v7.1.xml'

#download file
fpath = tempfile()
bfc = BiocFileCache(fpath, ask = FALSE)
msigpath = bfcrpath(bfc, msigdb_url)

#----MSigDB v7.1 Symbols----
#read genesets into a geneset collection
msigdb = getBroadSets(msigpath, membersId = 'MEMBERS_SYMBOLIZED')

#remove empty gene sets
msigdb_lengths = sapply(lapply(msigdb, geneIds), length)
msigdb = GeneSetCollection(msigdb[msigdb_lengths > 0])
save(msigdb, file = 'msigdb.v7.1.SYM.rda')

#----MSigDB v7.1 Entrez IDs----
#read genesets into a geneset collection
msigdb = getBroadSets(msigpath, membersId = 'MEMBERS_EZID')

#remove empty gene sets
msigdb_lengths = sapply(lapply(msigdb, geneIds), length)
msigdb = GeneSetCollection(msigdb[msigdb_lengths > 0])
save(msigdb, file = 'msigdb.v7.1.EZID.rda')
