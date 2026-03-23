library(dplyr)
library(tidyr)
library(WorldFlora)
library(data.table) # fread 用に追加

# =========================================================
# 初期設定
# =========================================================
##### 20250812
dir.create("Data2", showWarnings = FALSE) # Re-analysis for revise
dir.create("Data3", showWarnings = FALSE) # 最終出力用
file.copy("Data/GPT4o_all_fin.csv", "Data2/GPT4o_all_fin.csv", overwrite = TRUE) 
file.copy("Data/TRY_data_categ_fin.csv", "Data2/TRY_data_categ_fin.csv", overwrite = TRUE) 
Out_dir <- "Output2/"

# WFO (World Flora Online) データの設定
Taxa_file <- "Data2/classification.tsv" # version : v.2025.06

# =========================================================
# ヘルパー関数 (first_non_na_col, aggregate_by_key_chr2)
# =========================================================

first_non_na_col <- function(df, cols, fallback = NULL) {
  avail <- cols[cols %in% names(df)]
  out <- rep(NA_character_, nrow(df))
  if (length(avail) > 0) {
    for (nm in avail) {
      pick <- is.na(out) | out == ""
      v <- as.character(df[[nm]])
      out[pick] <- v[pick]
    }
  }
  if (!is.null(fallback) && fallback %in% names(df)) {
    pick <- is.na(out) | out == ""
    out[pick] <- as.character(df[[fallback]])[pick]
  }
  out
}

aggregate_by_key_chr2 <- function(df, key = "AccSpeciesName", str_cols = character(), bin_cols = character(), sep = ":") {
  missing <- setdiff(c(key, str_cols, bin_cols), names(df))
  if (length(missing)) stop("見つからない列があります: ", paste(missing, collapse = ", "))
  dup <- intersect(str_cols, bin_cols)
  if (length(dup)) stop("str_cols と bin_cols の両方に含まれる列: ", paste(dup, collapse = ", "))

  rest_cols <- setdiff(names(df), c(key, str_cols, bin_cols))
  string_agg <- function(x, sep = ":") {
    vals <- unique(stats::na.omit(as.character(x)))
    if (length(vals) == 0) NA_character_ else paste(vals, collapse = sep)
  }
  bin_agg <- function(x) as.integer(any(x %in% c(1L, "1", TRUE), na.rm = TRUE))

  df_first <- df %>% dplyr::group_by(.data[[key]]) %>% dplyr::slice(1L) %>% dplyr::ungroup()
  df_first <- df_first[, c(key, rest_cols), drop = FALSE]
  out <- df_first

  if (length(str_cols) > 0) {
    str_summ <- df %>% dplyr::group_by(.data[[key]]) %>%
      dplyr::summarise(dplyr::across(tidyselect::all_of(str_cols), ~ string_agg(.x, sep = sep)), .groups = "drop")
    out <- dplyr::left_join(out, str_summ, by = key)
  }
  if (length(bin_cols) > 0) {
    bin_summ <- df %>% dplyr::group_by(.data[[key]]) %>%
      dplyr::summarise(dplyr::across(tidyselect::all_of(bin_cols), bin_agg), .groups = "drop")
    out <- dplyr::left_join(out, bin_summ, by = key)
  }
  out
}

# =========================================================
# 前処理・マッチング関数
# =========================================================

prepare_gpt_data_exact <- function(GPT_file) {
  GPT_data <- read.csv(GPT_file, header = TRUE, row.names = 1)
  GPT_data[grepl("Psilotrichopsis curtisii", GPT_data$Species), ]$Species <- "Psilotrichopsis curtisii var. hainanensis"
  GPT_data[grepl("Arabidopsis halleri", GPT_data$Species), ]$Species <- "Arabidopsis halleri subsp. gemmifera"
  GPT_data[grepl("Arabidopsis lyrata", GPT_data$Species), ]$Species <- "Arabidopsis lyrata subsp. kamchatica"
  GPT_data$Species <- gsub("[[:space:]][^[:space:]][[:space:]]", " x ", GPT_data$Species)
  GPT_data
}

process_wfo_matching_exact <- function(input_df, spec_col, output_path, taxa_file) {
  clean_df <- WFO.match.fuzzyjoin(spec.data = input_df[[spec_col]], WFO.file = taxa_file, stringdist.method = "lv", fuzzydist.max = 4, spec.name.nobrackets = TRUE, spec.name.sub = TRUE)
  clean_df[clean_df$Hybrid != "×", ]$Hybrid <- "NotHybrid"
  clean_df[clean_df$Hybrid == "×", ]$Hybrid <- "Hybrid"
  clean_df$spec.name <- gsub("×", " x ", clean_df$spec.name)
  clean_df$scientificName <- gsub("×", "x", clean_df$scientificName)

  match_all <- clean_df %>% mutate(spec.name = as.character(spec.name)) %>%
    left_join(input_df %>% mutate(!!spec_col := as.character(.data[[spec_col]])), by = setNames(spec_col, "spec.name.ORIG"))

  df_matched <- match_all %>% filter((Matched == TRUE) & (taxonRank == "species"))
  dfm <- df_matched
  dfm$wfo_key <- first_non_na_col(dfm, c("acceptedNameUsageID", "scientificNameID", "taxonID"), fallback = "scientificName")
  dfm$wfo_name <- if ("scientificName" %in% names(dfm)) dfm$scientificName else dfm$wfo_key

  by_spec <- dfm %>% group_by(spec.name) %>% summarise(n_keys = n_distinct(wfo_key, na.rm = TRUE), .groups = "drop")
  by_key  <- dfm %>% group_by(wfo_key)   %>% summarise(n_specs = n_distinct(spec.name), .groups = "drop")

  df_one_to_one <- dfm %>% inner_join(filter(by_spec, n_keys == 1), by = "spec.name") %>% inner_join(filter(by_key, n_specs == 1), by = "wfo_key") %>% distinct(spec.name, wfo_key, .keep_all = TRUE)
  df_one_to_many <- dfm %>% semi_join(filter(by_spec, n_keys > 1), by = "spec.name")
  df_many_to_one <- dfm %>% semi_join(filter(by_key, n_specs > 1), by = "wfo_key")

  phen_cols <- intersect(names(df_many_to_one), c("White", "Yellow", "Red", "Blue", "Purple", "Green"))
  df_many_to_one_collapsed <- if (length(phen_cols) > 0) {
    df_many_to_one %>% group_by(wfo_key, wfo_name) %>%
      summarise(across(all_of(phen_cols), ~ any(. %in% c(1, TRUE, "1", "TRUE"), na.rm = TRUE)), species_merged = paste(sort(unique(spec.name)), collapse = "; "), n_collapsed = n_distinct(spec.name), .groups = "drop")
  } else { tibble() }

  one_to_one_full <- dfm %>% semi_join(df_one_to_one, by = c("spec.name", "wfo_key", "wfo_name")) %>% mutate(Combine_status = "one_to_one", Original_species = spec.name)
  one_to_many_full <- dfm %>% semi_join(df_one_to_many, by = c("spec.name", "wfo_key", "wfo_name")) %>% mutate(Combine_status = "one_to_many", Original_species = spec.name)

  if (nrow(df_many_to_one_collapsed) > 0) {
    collapsed_join <- df_many_to_one_collapsed[, unique(c("wfo_name", setdiff(names(df_many_to_one_collapsed), "wfo_id"))), drop = FALSE]
    many_to_one_full <- dfm %>% select(-any_of(setdiff(names(collapsed_join), "wfo_name"))) %>%
      inner_join(collapsed_join, by = "wfo_name") %>% mutate(Combine_status = "many_to_one", Original_species = species_merged)
  } else {
    many_to_one_full <- mutate(dfm[0,], Combine_status=character(), Original_species=character())
  }

  pick_many <- many_to_one_full %>% distinct(wfo_name, .keep_all = TRUE)
  pick_1_1  <- one_to_one_full %>% anti_join(pick_many, by = "wfo_name") %>% distinct(wfo_name, .keep_all = TRUE)
  pick_1_N  <- one_to_many_full %>% anti_join(pick_many, by = "wfo_name") %>% anti_join(pick_1_1, by = "wfo_name") %>% distinct(wfo_name, .keep_all = TRUE)
  
  combined_df <- bind_rows(pick_many, pick_1_1, pick_1_N) %>% relocate(Combine_status, Original_species, .after = wfo_name)
  write.csv(combined_df, output_path)
  
  return(list(summary_tbl = tibble::tibble(total = n_distinct(dfm$spec.name)), combined_df = combined_df))
}

process_try_dataset_exact <- function(TRY_data, dataset_name, output_path, taxa_file, str_cols, bin_cols) {
  sub_df <- TRY_data[TRY_data$Dataset == dataset_name, ]
  sub_df <- aggregate_by_key_chr2(sub_df, key = "AccSpeciesName", str_cols = str_cols, bin_cols = bin_cols)
  process_wfo_matching_exact(sub_df, "AccSpeciesName", output_path, taxa_file)
}

# =========================================================
# 実行セクション (GPT & TRY)
# =========================================================

GPT_data <- prepare_gpt_data_exact("Data2/GPT4o_all_fin.csv")
res_gpt <- process_wfo_matching_exact(GPT_data, "Species", "Data2/GPT4o_wfochecked.csv", Taxa_file)

TRY_data <- read.csv("Data2/TRY_data_categ_fin.csv", header = TRUE, row.names = 1)
str_cols <- c("OrigValueStr", "GPT_return2"); bin_cols <- c("White", "Yellow", "Red", "Blue", "Purple", "Green", "NoDescription")
try_targets <- c("PLANTSdata USDA" = "Data2/TRY_USDA_wfochecked.csv", "BiolFlor Database" = "Data2/TRY_BiolFlor_wfochecked.csv", "Orchid Trait Dataset" = "Data2/TRY_Orchid_wfochecked.csv", "Trait Data for African Plants - a Photo Guide" = "Data2/TRY_African_wfochecked.csv", "Functional Flowering Plant Traits" = "Data2/TRY_Functional_wfochecked.csv", "Traits of urban species from Ibagu\x81EColombia" = "Data2/TRY_Urban_wfochecked.csv")

res_try <- lapply(names(try_targets), function(ds) process_try_dataset_exact(TRY_data, ds, try_targets[[ds]], Taxa_file, str_cols, bin_cols))
names(res_try) <- names(try_targets)

# =========================================================
# Data Check (被子植物判定と最終保存)
# =========================================================

###
# Data Check
###

## Data was checked by wfo
# TRY data
Data_Urban <- fread("Data2/TRY_Urban_wfochecked.csv", data.table = F)
Data_Funct <- fread("Data2/TRY_Functional_wfochecked.csv", data.table = F)
Data_Africa <- fread("Data2/TRY_African_wfochecked.csv", data.table = F)
Data_Orchid <- fread("Data2/TRY_Orchid_wfochecked.csv", data.table = F)
Data_Biol <- fread("Data2/TRY_BiolFlor_wfochecked.csv", data.table = F)
Data_USDA <- fread("Data2/TRY_USDA_wfochecked.csv", data.table = F)

# GPT data
Data_GPT <- fread("Data2/GPT4o_wfochecked.csv", data.table = F)

# Combine
SpeciesKey_list <-unique(c(Data_Urban$wfo_key, Data_Funct$wfo_key, Data_Africa$wfo_key, Data_Orchid$wfo_key, Data_Biol$wfo_key, Data_USDA$wfo_key, Data_GPT$wfo_key))
Species_list <-unique(c(Data_Urban$wfo_name, Data_Funct$wfo_name, Data_Africa$wfo_name, Data_Orchid$wfo_name, Data_Biol$wfo_name, Data_USDA$wfo_name, Data_GPT$wfo_name))

# re-check taxa
Taxa_data <- fread("Data2/classification.tsv", data.table = F)
Taxa_data_filt_key <- Taxa_data[Taxa_data$scientificNameID %in% SpeciesKey_list, ]
Taxa_data_filt_name <- Taxa_data[Taxa_data$scientificName %in% Species_list, ]

Taxa_data_filt_key_A <- Taxa_data_filt_key[Taxa_data_filt_key$majorGroup == "A", ]
Taxa_data_filt_name_A <- Taxa_data_filt_name[Taxa_data_filt_name$majorGroup == "A", ]

SpeciesKey_list_A <- SpeciesKey_list[SpeciesKey_list %in% Taxa_data_filt_key_A$scientificNameID]
Species_list_A <- Species_list[Species_list %in% Taxa_data_filt_name_A$scientificName]

# Group detection
datasets <- list(Data_GPT, Data_Urban, Data_Funct, Data_Africa, Data_Orchid, Data_Biol, Data_USDA)
names(datasets) <- c("GPT", "Urban", "Funct", "Africa", "Orchid", "Biol", "USDA")

# 各データセットに Group 列を追加
processed_datasets <- lapply(datasets, function(df) {
  df$Group <- "OtherOrNotInWfo"
  df[df$wfo_name %in% Species_list_A | df$wfo_key %in% SpeciesKey_list_A, ]$Group <- "Angiosperms"
  return(df)
})

# 個別の変数に書き戻し（write.csv 用）
Data_GPT <- processed_datasets$GPT
Data_Urban <- processed_datasets$Urban
Data_Funct <- processed_datasets$Funct
Data_Africa <- processed_datasets$Africa
Data_Orchid <- processed_datasets$Orchid
Data_Biol <- processed_datasets$Biol
Data_USDA <- processed_datasets$USDA

# 最終保存
write.csv(Data_GPT, "Data3/GPT4o_fin.csv")
write.csv(Data_Urban, "Data3/TRY_Urban_fin.csv")
write.csv(Data_Funct, "Data3/TRY_Functional_fin.csv")
write.csv(Data_Africa, "Data3/TRY_Africa_fin.csv")
write.csv(Data_Orchid, "Data3/TRY_Orchid_fin.csv")
write.csv(Data_Biol, "Data3/TRY_Biol_fin.csv")
write.csv(Data_USDA, "Data3/TRY_USDA_fin.csv")
