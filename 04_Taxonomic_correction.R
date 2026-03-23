library(dplyr)
library(tidyr)
library(WorldFlora)

# Comments have been added by ChatGPT 5.4 Thinking
# =========================================================
# 初期設定
# =========================================================
# - Data2 ディレクトリを作成し、元データを再解析用にコピーする
# - 解析結果は test1 と同じ挙動になるように Data2 / Output2 を前提とする

##### 20250812
dir.create("Data2", showWarnings = FALSE) # Re-analysis for revise
file.copy("Data/GPT4o_all_fin.csv", "Data2/GPT4o_all_fin.csv", overwrite = TRUE) # From Flora of China data
file.copy("Data/TRY_data_categ_fin.csv", "Data2/TRY_data_categ_fin.csv", overwrite = TRUE) # From TRY data
Out_dir <- "Output2/"

# =========================================================
# WFO (World Flora Online) データの設定
# =========================================================
# taxa information was arranged with World flora online data.
# classification.tsv は World Flora Online から取得した分類データを想定
Taxa_file <- "Data2/classification.tsv" # Downloaded from https://www.worldfloraonline.org/downloadData (version : v.2025.06)

# =========================================================
# ヘルパー関数
# =========================================================

# ---------------------------------------------------------
# first_non_na_col()
# ---------------------------------------------------------
# 目的:
#   指定した複数列の中から、各行について最初に見つかった非 NA / 非空文字の値を返す。
# 用途:
#   WFO のキー列候補（acceptedNameUsageID, scientificNameID, taxonID など）から
#   優先順位付きで代表値を1つ選ぶために使う。
#
# 引数:
#   df       : 対象データフレーム
#   cols     : 優先順に並べた候補列名ベクトル
#   fallback : cols で埋まらなかったときに使う予備列名
#
# 戻り値:
#   各行に対して 1 つの文字列値を持つ character ベクトル
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

# ---------------------------------------------------------
# aggregate_by_key_chr2()
# ---------------------------------------------------------
# 目的:
#   同じ種名キーを持つ複数行を 1 行に集約する。
#
# 集約ルール:
#   1) key ごとに「その他の列」は先頭行の値を採用
#   2) str_cols は NA を除き、重複を除いて sep で結合
#   3) bin_cols は any() により 0/1 に集約
#
# 用途:
#   TRY データでは同一種に複数レコードがあるため、WFO マッチング前に
#   test1 と同様の形で種単位へまとめる。
aggregate_by_key_chr2 <- function(df,
                                  key = "AccSpeciesName",
                                  str_cols = character(),
                                  bin_cols = character(),
                                  sep = ":") {

  # 入力チェック
  missing <- setdiff(c(key, str_cols, bin_cols), names(df))
  if (length(missing)) stop("見つからない列があります: ", paste(missing, collapse = ", "))

  dup <- intersect(str_cols, bin_cols)
  if (length(dup)) stop("str_cols と bin_cols の両方に含まれる列: ", paste(dup, collapse = ", "))

  # “その他の列” = 集約対象ではなく、先頭行の値を保持する列
  rest_cols <- setdiff(names(df), c(key, str_cols, bin_cols))

  # 文字列列の連結（重複・NA 除去、初出順を維持）
  string_agg <- function(x, sep = ":") {
    vals <- unique(stats::na.omit(as.character(x)))
    if (length(vals) == 0) NA_character_ else paste(vals, collapse = sep)
  }

  # 0/1 列の集約（1 つでも TRUE/1 があれば 1）
  bin_agg <- function(x) as.integer(any(x %in% c(1L, "1", TRUE), na.rm = TRUE))

  # 1) “その他の列”は key ごとに先頭行を採用
  df_first <- df %>%
    dplyr::group_by(.data[[key]]) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()

  # select() を使わず base subset を使用（名前マスク回避）
  df_first <- df_first[, c(key, rest_cols), drop = FALSE]

  out <- df_first

  # 2) 文字列列を結合
  if (length(str_cols) > 0) {
    str_summ <- df %>%
      dplyr::group_by(.data[[key]]) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::all_of(str_cols), ~ string_agg(.x, sep = sep)),
        .groups = "drop"
      )
    out <- dplyr::left_join(out, str_summ, by = key)
  }

  # 3) 0/1 列を集約
  if (length(bin_cols) > 0) {
    bin_summ <- df %>%
      dplyr::group_by(.data[[key]]) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::all_of(bin_cols), bin_agg),
        .groups = "drop"
      )
    out <- dplyr::left_join(out, bin_summ, by = key)
  }

  out
}

# ---------------------------------------------------------
# prepare_gpt_data_exact()
# ---------------------------------------------------------
# 目的:
#   GPT 側データに対して、test1 と同じ手順で種名表記を前処理する。
#
# 処理内容:
#   - 手動修正が必要な Species 表記を置換
#   - hybrid 表記ゆれを " x " に統一
#   - 途中確認のため print() を残す（test1 互換）
prepare_gpt_data_exact <- function(GPT_file) {
  GPT_data <- read.csv(GPT_file, header = TRUE, row.names = 1)

  ## Steps to manually address variations in Species name notation in the GPT4o file
  print(GPT_data[grepl("\\(", GPT_data$Species), ]$Species)
  ### "Psilotrichopsis curtisii (Oliver) C. C. Townsend var. hainanensis"
  ### "Arabidopsis halleri (Linnaeus) O’Kane & Al-Shehbaz subsp. gemmifera"
  ### "Arabidopsis lyrata (Linnaeus) O’Kane & Al-Shehbaz subsp. kamchatica"
  GPT_data[grepl("Psilotrichopsis curtisii", GPT_data$Species), ]$Species <- "Psilotrichopsis curtisii var. hainanensis"
  GPT_data[grepl("Arabidopsis halleri", GPT_data$Species), ]$Species <- "Arabidopsis halleri subsp. gemmifera"
  GPT_data[grepl("Arabidopsis lyrata", GPT_data$Species), ]$Species <- "Arabidopsis lyrata subsp. kamchatica"

  ### Hybrid species name variant
  GPT_data$Species <- gsub("[[:space:]][^[:space:]][[:space:]]", " x ", GPT_data$Species)
  print(GPT_data[grepl("[[:space:]][^[:space:]][[:space:]]", GPT_data$Species), ]) # Check

  GPT_data
}

# ---------------------------------------------------------
# process_wfo_matching_exact()
# ---------------------------------------------------------
# 目的:
#   1つの入力データフレームに対して WFO マッチングを行い、
#   test1 と同じロジックで matched / unmatched / combined_df を作成する。
#
# 主な処理:
#   1) WFO.match.fuzzyjoin() により種名マッチング
#   2) Hybrid 表記や scientificName の表記ゆれを補正
#   3) 元データを join して match_all を作成
#   4) species / non-species / unmatched に分割
#   5) 1対1, 1対多, 多対1 の関係を分類
#   6) test1 と同じ優先順位で combined_df を構築
#   7) CSV を test1 と同じ write.csv() 仕様で出力
process_wfo_matching_exact <- function(input_df, spec_col, output_path, taxa_file) {
  clean_df <- WFO.match.fuzzyjoin(
    spec.data = input_df[[spec_col]],
    WFO.file  = taxa_file,
    stringdist.method = "lv",
    fuzzydist.max = 4,
    spec.name.nobrackets = TRUE,
    spec.name.sub = TRUE
  )

  # Hybrid 列の表記を test1 と同じラベルへ揃える
  clean_df[clean_df$Hybrid != "×", ]$Hybrid <- "NotHybrid"
  clean_df[clean_df$Hybrid == "×", ]$Hybrid <- "Hybrid"

  # 表記ゆれの統一
  clean_df$spec.name <- gsub("×", " x ", clean_df$spec.name)
  clean_df$scientificName <- gsub("×", "x", clean_df$scientificName)

  # 元入力データを spec.name.ORIG ベースで再結合
  match_all <- clean_df %>%
    mutate(spec.name = as.character(spec.name)) %>%
    left_join(
      input_df %>% mutate(!!spec_col := as.character(.data[[spec_col]])),
      by = setNames(spec_col, "spec.name.ORIG")
    )

  # マッチ結果の区分け
  df_matched   <- match_all %>% filter((Matched == TRUE) & (taxonRank == "species"))
  df_unmatched <- match_all %>% filter(is.na(Matched) | Matched == FALSE) %>%
    distinct(spec.name, .keep_all = TRUE)
  df_nospecies <- match_all %>% filter(taxonRank != "species") %>%
    distinct(spec.name, .keep_all = TRUE)

  # matched のみを対象に、後続の整理用列を追加
  dfm <- df_matched
  dfm$wfo_key     <- first_non_na_col(dfm,
                                      c("acceptedNameUsageID", "scientificNameID", "taxonID"),
                                      fallback = "scientificName")
  dfm$wfo_name    <- if ("scientificName" %in% names(dfm)) dfm$scientificName else dfm$wfo_key
  dfm$genus_input <- sub(" .*$", "", as.character(dfm$spec.name))
  dfm$genus_wfo   <- if ("scientificName" %in% names(dfm)) sub(" .*$", "", dfm$scientificName) else NA_character_
  dfm$fuzzy_dist  <- if ("Fuzzy.dist" %in% names(dfm)) suppressWarnings(as.numeric(dfm$Fuzzy.dist)) else NA_real_

  # species -> key, key -> species の対応数を集計
  by_spec <- dfm %>% group_by(spec.name) %>% summarise(n_keys = n_distinct(wfo_key, na.rm = TRUE), .groups = "drop")
  by_key  <- dfm %>% group_by(wfo_key)   %>% summarise(n_specs = n_distinct(spec.name), .groups = "drop")

  # 1対1: species も key も一意
  df_one_to_one <- dfm %>%
    inner_join(filter(by_spec, n_keys == 1), by = "spec.name") %>%
    inner_join(filter(by_key,  n_specs == 1), by = "wfo_key") %>%
    distinct(spec.name, wfo_key, .keep_all = TRUE)

  # 1対多: 1 species に複数 key がぶら下がる
  df_one_to_many <- dfm %>%
    semi_join(filter(by_spec, n_keys > 1), by = "spec.name") %>%
    arrange(spec.name, wfo_key)

  # 多対1: 複数 species が同じ key に集約される
  df_many_to_one <- dfm %>%
    semi_join(filter(by_key, n_specs > 1), by = "wfo_key") %>%
    arrange(wfo_key, spec.name)

  # 色カテゴリなどの phenotype 列を many_to_one 用に統合
  phen_cols <- intersect(names(df_many_to_one), c("White", "Yellow", "Red", "Blue", "Purple", "Green"))
  df_many_to_one_collapsed <- if (length(phen_cols) > 0) {
    df_many_to_one %>%
      group_by(wfo_key, wfo_name) %>%
      summarise(
        across(all_of(phen_cols), ~ any(. %in% c(1, TRUE, "1", "TRUE"), na.rm = TRUE)),
        species_merged = paste(sort(unique(spec.name)), collapse = "; "),
        n_collapsed    = n_distinct(spec.name),
        .groups = "drop"
      )
  } else {
    tibble()
  }

  # 主要件数のサマリー
  summary_tbl <- tibble::tibble(
    total_input          = dplyr::n_distinct(dfm$spec.name),
    one_to_one_pairs     = nrow(df_one_to_one),
    one_to_many_species  = dplyr::n_distinct(df_one_to_many$spec.name),
    many_to_one_targets  = dplyr::n_distinct(df_many_to_one$wfo_key)
  )

  # 以降の結合で必要な列が存在することを確認
  stopifnot(all(c("spec.name", "wfo_name", "wfo_key") %in% names(dfm)))

  # 1対1データにステータス列を付与
  one_to_one_full <- dfm %>%
    semi_join(df_one_to_one %>% distinct(spec.name, wfo_key, wfo_name),
              by = c("spec.name", "wfo_key", "wfo_name")) %>%
    mutate(
      Combine_status   = "one_to_one",
      Original_species = spec.name
    )

  # 1対多データにステータス列を付与
  one_to_many_full <- dfm %>%
    semi_join(df_one_to_many %>% distinct(spec.name, wfo_key, wfo_name),
              by = c("spec.name", "wfo_key", "wfo_name")) %>%
    mutate(
      Combine_status   = "one_to_many",
      Original_species = spec.name
    )

  # 多対1データは species_merged を持つ collapsed テーブルと結合し、
  # 代表行を 1 つにまとめるための準備を行う
  if (nrow(df_many_to_one_collapsed) > 0) {
    extra_cols <- setdiff(names(df_many_to_one_collapsed), "wfo_id")
    keep_cols  <- unique(c("wfo_name", extra_cols))
    collapsed_join <- df_many_to_one_collapsed[, keep_cols, drop = FALSE]

    dup_cols <- setdiff(intersect(names(dfm), names(collapsed_join)), "wfo_name")
    dfm_trim <- dfm[, setdiff(names(dfm), dup_cols), drop = FALSE]

    many_to_one_full <- dfm_trim %>%
      dplyr::inner_join(collapsed_join, by = "wfo_name") %>%
      dplyr::mutate(
        Combine_status   = "many_to_one",
        Original_species = if ("species_merged" %in% names(.)) species_merged else NA_character_
      )
  } else {
    many_to_one_full <- dplyr::slice_head(dfm, n = 0) %>%
      dplyr::mutate(Combine_status = character(), Original_species = character())
  }

  # test1 と同じ優先順位で combined_df を作る
  # 優先順位: many_to_one -> one_to_one -> one_to_many
  # さらに wfo_name ごとに代表行を 1 行にする
  pick_many <- many_to_one_full %>%
    dplyr::distinct(wfo_name, .keep_all = TRUE)

  pick_one_to_one <- one_to_one_full %>%
    dplyr::anti_join(pick_many %>% dplyr::select(wfo_name), by = "wfo_name") %>%
    dplyr::distinct(wfo_name, .keep_all = TRUE)

  pick_one_to_many <- one_to_many_full %>%
    dplyr::anti_join(pick_many %>% dplyr::select(wfo_name), by = "wfo_name") %>%
    dplyr::anti_join(pick_one_to_one %>% dplyr::select(wfo_name), by = "wfo_name") %>%
    dplyr::distinct(wfo_name, .keep_all = TRUE)

  combined_df <- dplyr::bind_rows(pick_many, pick_one_to_one, pick_one_to_many) %>%
    dplyr::relocate(Combine_status, Original_species, .after = wfo_name)

  # 件数確認（test1 互換で値をコンソールに出す）
  unique(combined_df$Combine_status)
  sum(combined_df$Combine_status == "one_to_one")
  sum(combined_df$Combine_status == "one_to_many")
  sum(combined_df$Combine_status == "many_to_one")

  # 出力は test1 と同じ write.csv() 仕様
  write.csv(combined_df, output_path)

  list(
    clean_df = clean_df,
    match_all = match_all,
    df_matched = df_matched,
    df_unmatched = df_unmatched,
    df_nospecies = df_nospecies,
    summary_tbl = summary_tbl,
    combined_df = combined_df
  )
}

# ---------------------------------------------------------
# process_try_dataset_exact()
# ---------------------------------------------------------
# 目的:
#   TRY データのうち特定 Dataset 名に対応する部分集合を取り出し、
#   種単位に集約してから WFO マッチングへ流す。
#
# 処理内容:
#   1) Dataset 名でフィルタ
#   2) AccSpeciesName 単位に集約
#   3) process_wfo_matching_exact() を呼び出して出力
process_try_dataset_exact <- function(TRY_data, dataset_name, output_path, taxa_file, str_cols, bin_cols) {
  sub_df <- TRY_data[TRY_data$Dataset == dataset_name, ]
  sub_df <- aggregate_by_key_chr2(
    df       = sub_df,
    key      = "AccSpeciesName",
    str_cols = str_cols,
    bin_cols = bin_cols
  )

  process_wfo_matching_exact(
    input_df    = sub_df,
    spec_col    = "AccSpeciesName",
    output_path = output_path,
    taxa_file   = taxa_file
  )
}

# =========================================================
# GPT result (from Flora of China)
# =========================================================
# 1) GPT 出力 CSV を読み込む
# 2) Species 表記を test1 と同じ前処理で補正
# 3) WFO マッチング結果を CSV と list で出力

GPT_file <- "Data2/GPT4o_all_fin.csv"
GPT_data <- prepare_gpt_data_exact(GPT_file)
res_gpt <- process_wfo_matching_exact(
  input_df    = GPT_data,
  spec_col    = "Species",
  output_path = "Data2/GPT4o_wfochecked.csv",
  taxa_file   = Taxa_file
)

# =========================================================
# For TRY data
# =========================================================
# 1) TRY 元データを読み込む
# 2) 集約対象の文字列列 / 0-1 列を定義する
# 3) 対象 Dataset ごとに WFO マッチングを実行する

TRY_file <- "Data2/TRY_data_categ_fin.csv"
TRY_data <- read.csv(TRY_file, header = TRUE, row.names = 1)

# 同一 species に複数値がある場合は文字列連結する列
str_cols <- c("OrigValueStr", "GPT_return2")

# 同一 species に複数値がある場合は any() でまとめる列
bin_cols <- c("White", "Yellow", "Red", "Blue", "Purple", "Green", "NoDescription")

# 読み込まれた Dataset 名の確認
print(unique(TRY_data$Dataset))

# Dataset 名と出力先ファイルの対応表
try_targets <- c(
  "PLANTSdata USDA" = "Data2/TRY_USDA_wfochecked.csv",
  "BiolFlor Database" = "Data2/TRY_BiolFlor_wfochecked.csv",
  "Orchid Trait Dataset" = "Data2/TRY_Orchid_wfochecked.csv",
  "Trait Data for African Plants - a Photo Guide" = "Data2/TRY_African_wfochecked.csv",
  "Functional Flowering Plant Traits" = "Data2/TRY_Functional_wfochecked.csv",
  "Traits of urban species from Ibagu\x81EColombia" = "Data2/TRY_Urban_wfochecked.csv"
)

# 各 Dataset について順番に同じ処理を行う
res_try <- lapply(names(try_targets), function(dataset_name) {
  process_try_dataset_exact(
    TRY_data     = TRY_data,
    dataset_name = dataset_name,
    output_path  = unname(try_targets[[dataset_name]]),
    taxa_file    = Taxa_file,
    str_cols     = str_cols,
    bin_cols     = bin_cols
  )
})

# list の名前を Dataset 名にして参照しやすくする
names(res_try) <- names(try_targets)
