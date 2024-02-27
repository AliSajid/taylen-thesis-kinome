# Creedenzymatic Analysis

library(tidyverse)
library(creedenzymatic)

process_creedenzymatic <- function(krsa_path, uka_path, peptide_path) {
  krsa_data <- read_csv(krsa_path, show_col_types = FALSE) |>
    select(Kinase, Score = AvgZ) |>
    read_krsa(trns = "abs", sort = "desc")

  uka_data <- read_tsv(uka_path, show_col_types = FALSE) |>
    select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
    read_uka(trns = "abs", sort = "desc")

  peptide_data <- read_csv(peptide_path, show_col_types = FALSE) |>
    select(Peptide, Score = totalMeanLFC)

  kea3_data <- read_kea(
    peptide_data,
    sort = "asc",
    trns = "abs",
    method = "MeanRank",
    lib = "kinase-substrate"
  )

  ptmsea_data <- read_ptmsea(peptide_data)

  combined <- combine_tools(
    KRSA_df = krsa_data,
    UKA_df = uka_data,
    KEA3_df = kea3_data,
    PTM_SEA_df = ptmsea_data
  )

  combined
}

krsa_files <- c(
  "results//krsa-krsa_table_full_SN_CN_STK.csv",
  "results//krsa-krsa_table_full_CBT_CN_STK.csv",
  "results//krsa-krsa_table_full_SBT_CN_STK.csv"
)

uka_files <- c(
  "kinome_data/UKA/SN-CN/Summaryresults 20231024-1011.txt",
  "kinome_data/UKA/CBT-CN/Summaryresults 20231024-1014.txt",
  "kinome_data/UKA/CN-SBT/Summaryresults 20231024-1016.txt"
)

peptide_files <- c(
  "results/krsa-dpp_SN_CN-STK.csv",
  "results/krsa-dpp_CBT_CN-STK.csv",
  "results/krsa-dpp_SBT_CN-STK.csv"
)

result <-
  list(
    krsa_path = krsa_files,
    uka_path = uka_files,
    peptide_path = peptide_files
  ) |>
  pmap(process_creedenzymatic) |>
  set_names(c("SN-CN", "CBT-CN", "SBT-CN")) |>
  imap_dfr(~ write_csv(.x, str_glue("results/{.y}_creedenzymatic.csv")), .id = "Comparison")
