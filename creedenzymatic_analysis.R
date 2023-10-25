# Creedenzymatic Analysis

library(tidyverse)
library(creedenzymatic)

process_creedenzymatic <-
  function(krsa_path, uka_path, peptide_path) {
    krsa_data <- read_csv(krsa_path, show_col_types = FALSE) |>
      select(Kinase, Score = AvgZ) |>
      read_krsa(trns = "abs", sort = "desc")

    uka_data <- read_tsv(uka_path, show_col_types = FALSE) |>
      select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
      read_uka(trns = "abs", sort = "desc")

    peptide_data <-
      read_csv(peptide_path, show_col_types = FALSE) |>
      select(Peptide, Score = totalMeanLFC)

    kea3_data <-
      read_kea(
        peptide_data,
        sort = "asc",
        trns = "abs",
        method = "MeanRank",
        lib = "kinase-substrate"
      )

    ptmsea_data <-
      read_ptmsea(peptide_data)

    combined <- combine_tools(
      KRSA_df = krsa_data,
      UKA_df = uka_data,
      KEA3_df = kea3_data,
      PTM_SEA_df = ptmsea_data
    )

    combined
  }

krsa_files <- c(
  "results//krsa-krsa_table_full_CMS_CTL_STK.csv",
  "results//krsa-krsa_table_full_CPP_CTL_STK.csv",
  "results//krsa-krsa_table_full_SPP_CTL_STK.csv"
)

uka_files <- c(
  "kinome_data/UKA/CMS-CTL/Summaryresults 20231024-1011.txt",
  "kinome_data/UKA/CPP-CTL/Summaryresults 20231024-1014.txt",
  "kinome_data/UKA/CTL-SPP/Summaryresults 20231024-1016.txt"
)

peptide_files <- c(
  "results/krsa-dpp_CMS_CTL-STK.csv",
  "results/krsa-dpp_CPP_CTL-STK.csv",
  "results/krsa-dpp_SPP_CTL-STK.csv"
)

result <-
  list(
    krsa_path = krsa_files,
    uka_path = uka_files,
    peptide_path = peptide_files
  ) |>
  pmap(process_creedenzymatic) |>
  set_names(c("CMS-CTL", "CPP-CTL", "SPP-CTL")) |>
  imap_dfr(~ write_csv(.x, str_glue("results/{.y}_creedenzymatic.csv")), .id = "Comparison")
