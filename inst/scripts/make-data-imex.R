library(plyr)
library(stringr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tidyverse)
library(BiocFileCache)

imex_date = '2022-07-11'
bfc = BiocFileCache()
imex_path = bfcrpath(bfc, sprintf('ftp://ftp.ebi.ac.uk/pub/databases/intact/%s/psimitab/intact-micluster.txt', imex_date))

#----download PPI----
#process human and mouse PPI
ppi = read.csv(imex_path, sep = '\t', header = TRUE)
ppi = ppi |> 
  #find interactions determined within the same species
  filter(Taxid.interactor.A == Taxid.interactor.B) |> 
  #extract values
  mutate(
    Taxid = str_extract(Taxid.interactor.A, '[0-9]+')
  ) |> 
  #choose human and mouse only
  filter(Taxid %in% c('9606', '10090')) |> 
  mutate(
    InteractionType = gsub('psi-mi:MI:[0-9]+\\(|\\)', '', Interaction.type.s.),
    DetectionMethod = gsub('psi-mi:MI:[0-9]+\\(|\\)', '', Interaction.detection.method.s.)
  ) |>
  #extract values
  mutate(
    Confidence = as.numeric(str_remove(Confidence.value.s., '^.*:')),
    InteractorA = str_remove(Alt..ID.s..interactor.A, '^.*:'),
    InteractorB = str_remove(Alt..ID.s..interactor.B, '^.*:')
  ) |> 
  #map uniprot to entrez
  mutate(
    EntrezA = case_when(
      Taxid == '9606' ~ mapIds(org.Hs.eg.db, InteractorA, 'ENTREZID', 'UNIPROT'),
      Taxid == '10090' ~ mapIds(org.Mm.eg.db, InteractorA, 'ENTREZID', 'UNIPROT'),
    ),
    EntrezB = case_when(
      Taxid == '9606' ~ mapIds(org.Hs.eg.db, InteractorB, 'ENTREZID', 'UNIPROT'),
      Taxid == '10090' ~ mapIds(org.Mm.eg.db, InteractorB, 'ENTREZID', 'UNIPROT'),
    )
  ) |> 
  filter(!is.na(EntrezA) & !is.na(EntrezB)) |> 
  #map uniprot to symbol
  mutate(
    SymbolA = case_when(
      Taxid == '9606' ~ mapIds(org.Hs.eg.db, EntrezA, 'SYMBOL', 'ENTREZID'),
      Taxid == '10090' ~ mapIds(org.Mm.eg.db, EntrezA, 'SYMBOL', 'ENTREZID'),
    ),
    SymbolB = case_when(
      Taxid == '9606' ~ mapIds(org.Hs.eg.db, EntrezB, 'SYMBOL', 'ENTREZID'),
      Taxid == '10090' ~ mapIds(org.Mm.eg.db, EntrezB, 'SYMBOL', 'ENTREZID'),
    )
  ) |> 
  distinct() |> 
  #summarise information for each pair
  select(
    InteractorA,
    InteractorB,
    EntrezA,
    EntrezB,
    SymbolA,
    SymbolB,
    Taxid,
    InteractionType,
    DetectionMethod,
    Confidence
  ) |> 
  mutate(Inferred = FALSE)

#split human and mouse
ppi_hs = ppi[ppi$Taxid %in% '9606', ]
ppi_mm = ppi[ppi$Taxid %in% '10090', ]

#----translate using homologs----
ppi_hs_hcop = ppi_mm |>
  left_join(msigdb:::hcop, by = c('EntrezA' = 'mouse_entrez_gene')) |> 
  left_join(msigdb:::hcop, by = c('EntrezB' = 'mouse_entrez_gene')) |> 
  ungroup() |>
  select(!(EntrezA:EntrezB)) |>
  rename(EntrezA = human_entrez_gene.x, EntrezB = human_entrez_gene.y) |>
  mutate(
    SymbolA = mapIds(org.Hs.eg.db, EntrezA, 'SYMBOL', 'ENTREZID'),
    SymbolB = mapIds(org.Hs.eg.db, EntrezB, 'SYMBOL', 'ENTREZID'),
    Taxid = '9606',
    Inferred = TRUE
  ) |>
  distinct() |> 
  as.data.frame()
ppi_hs = rbind(ppi_hs, ppi_hs_hcop)

ppi_mm_hcop = ppi_hs |>
  left_join(msigdb:::hcop, by = c('EntrezA' = 'human_entrez_gene')) |> 
  left_join(msigdb:::hcop, by = c('EntrezB' = 'human_entrez_gene')) |> 
  ungroup() |>
  select(!(EntrezA:EntrezB)) |>
  rename(EntrezA = mouse_entrez_gene.x, EntrezB = mouse_entrez_gene.y) |>
  mutate(
    SymbolA = mapIds(org.Mm.eg.db, EntrezA, 'SYMBOL', 'ENTREZID'),
    SymbolB = mapIds(org.Mm.eg.db, EntrezB, 'SYMBOL', 'ENTREZID'),
    Taxid = '10090',
    Inferred = TRUE
  ) |>
  distinct() |> 
  as.data.frame()
ppi_mm = rbind(ppi_mm, ppi_mm_hcop)
rm(ppi_mm_hcop, ppi_hs_hcop)

#----remove redundancry----
ppi = rbind(ppi_hs, ppi_mm) |> distinct()

ppi = ppi |> 
  mutate(switch = as.numeric(EntrezA) > as.numeric(EntrezB)) |> 
  mutate(
    tmpB = EntrezB,
    EntrezB = if_else(switch, EntrezA, EntrezB),
    EntrezA = if_else(switch, tmpB, EntrezA),
    tmpB = SymbolB,
    SymbolB = if_else(switch, SymbolA, SymbolB),
    SymbolA = if_else(switch, tmpB, SymbolA),
    tmpB = InteractorB,
    InteractorB = if_else(switch, InteractorA, InteractorB),
    InteractorA = if_else(switch, tmpB, InteractorA)
  ) |> 
  distinct() |> 
  group_by(EntrezA, EntrezB, Taxid) |> 
  summarise(
    InteractorA = paste(unique(InteractorA), collapse = '|'),
    InteractorB = paste(unique(InteractorB), collapse = '|'),
    SymbolA = paste(unique(SymbolA), collapse = '|'),
    SymbolB = paste(unique(SymbolB), collapse = '|'),
    InteractionType = paste(InteractionType, collapse = '|'),
    DetectionMethod = paste(DetectionMethod, collapse = '|'),
    Confidence = max(Confidence),
    Inferred = all(Inferred)
  ) |> 
  ungroup()

#remove duplicate evidence
ppi$InteractionType = sapply(str_split(ppi$InteractionType, '\\|'), function(x) {
  paste(unique(x), collapse = '|')
})
ppi$DetectionMethod = sapply(str_split(ppi$DetectionMethod, '\\|'), function(x) {
  paste(unique(x), collapse = '|')
})

saveRDS(ppi, file = sprintf('imex_hsmm_%s.rds', format(as.Date(imex_date), '%m%y')))
