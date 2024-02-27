# Pathway Analysis STK - Su Data

library(tidyverse)
library(enrichR)
library(writexl)

# Load the reference data
stk_id <- readRDS("reference_data/stk_id_map.Rds")
stk_hgnc <- readRDS("reference_data/stk_hgnc_map.Rds")

stk_map <- stk_id |>
  inner_join(stk_hgnc)

# Copy this section for each comparison
#################

# Change the name of the file you are loading
sn_comparison <- read_csv(
  file.path("results", "krsa-dpp_SN_CN-STK.csv")
) |>
  group_by(Peptide) |>
  filter(abs(LFC) == max(abs(LFC))) |>
  ungroup() |>
  select(Peptide, LFC) |>
  inner_join(stk_map) |>
  # Change the name of the file you are writing.
  write_csv(
    file.path("results", "annotated_dpp_SN_CN-STK.csv")
  )

sn_genes <- sn_comparison |>
  select(Gene, LFC) |>
  filter(
    LFC >= quantile(sn_comparison[["LFC"]], 0.90) |
      LFC <= quantile(sn_comparison[["LFC"]], 0.10)
  ) |>
  pull(Gene)

dbs <- c(
  "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022",
  "WikiPathways_2023_Human"
)

sn_enriched <- enrichr(sn_genes, dbs) |>
  imap(~ write_csv(.x, str_glue(
    file.path("results", "CSM-SN_CN-STK-{.y}-Pathways.csv")
  )))

# Change the name of the file you are writing
write_xlsx(
  sn_enriched,
  file.path("results", "SN-STK-Pathways.xlsx")
)
#################


# Change the name of the file you are loading
cbt_comparison <- read_csv(
  file.path("results", "krsa-dpp_CBT_CN-STK.csv")
) |>
  group_by(Peptide) |>
  filter(abs(LFC) == max(abs(LFC))) |>
  ungroup() |>
  select(Peptide, LFC) |>
  inner_join(stk_map) |>
  # Change the name of the file you are writing.
  write_csv(
    file.path("results", "annotated_dpp_CBT_CN-STK.csv")
  )

cbt_genes <- cbt_comparison |>
  select(Gene, LFC) |>
  filter(
    LFC >= quantile(cbt_comparison[["LFC"]], 0.90) |
      LFC <= quantile(cbt_comparison[["LFC"]], 0.10)
  ) |>
  pull(Gene)

dbs <- c(
  "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022",
  "WikiPathways_2023_Human"
)

cbt_enriched <- enrichr(cbt_genes, dbs) |>
  imap(~ write_csv(.x, str_glue(
    file.path("results", "CSM-CBT_CN-STK-{.y}-Pathways.csv")
  )))

# Change the name of the file you are writing
write_xlsx(
  cbt_enriched,
  file.path("results", "CBT-STK-Pathways.xlsx")
)
#################

# Change the name of the file you are loading
sbt_comparison <- read_csv(
  file.path("results", "krsa-dpp_SBT_CN-STK.csv")
) |>
  group_by(Peptide) |>
  filter(abs(LFC) == max(abs(LFC))) |>
  ungroup() |>
  select(Peptide, LFC) |>
  inner_join(stk_map) |>
  # Change the name of the file you are writing.
  write_csv(
    file.path("results", "annotated_dpp_SBT_CN-STK.csv")
  )

sbt_genes <- sbt_comparison |>
  select(Gene, LFC) |>
  filter(
    LFC >= quantile(sbt_comparison[["LFC"]], 0.90) |
      LFC <= quantile(sbt_comparison[["LFC"]], 0.10)
  ) |>
  pull(Gene)

dbs <- c(
  "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022",
  "WikiPathways_2023_Human"
)

sbt_enriched <- enrichr(sbt_genes, dbs) |>
  imap(~ write_csv(.x, str_glue(
    file.path("results", "CSM-SBT_CN-STK-{.y}-Pathways.csv")
  )))

# Change the name of the file you are writing
write_xlsx(
  sbt_enriched,
  file.path("results", "SBT-STK-Pathways.xlsx")
)
#################
