# ---- Packages ---------------------------------------------------
library(data.table)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggraph)
library(igraph)
library(gridExtra)

# ---- Paths and shared settings ---------------------------------
output_dir <- "Output_20260120"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

color_cols <- c("White", "Yellow", "Red", "Blue", "Purple", "Green")
node_order <- c("Yellow", "White", "Red", "Purple", "Blue", "Green")

# ---- Input data -------------------------------------------------
Data_GPT        <- read.csv("Data3/GPT4o_fin.csv", header = TRUE, row.names = 1)
Data_Biol       <- read.csv("Data3/TRY_Biol_fin.csv", header = TRUE, row.names = 1)
Data_Orchid     <- read.csv("Data3/TRY_Orchid_fin.csv", header = TRUE, row.names = 1)
Data_Urban      <- read.csv("Data3/TRY_Urban_fin.csv", header = TRUE, row.names = 1)
Data_USDA       <- read.csv("Data3/TRY_USDA_fin.csv", header = TRUE, row.names = 1)
Data_African    <- read.csv("Data3/TRY_Africa_fin.csv", header = TRUE, row.names = 1)
Data_Functional <- read.csv("Data3/TRY_Functional_fin.csv", header = TRUE, row.names = 1)
Data_TRYALL     <- read.csv("Data3/TRY_combined.csv", header = TRUE, row.names = 1)

Data_all <- fread("Data3/All_data.csv", data.table = FALSE)
Species_angiosperms <- Data_all[Data_all$Group == "Angiosperms", ]$wfo_name

# ---- Shared preprocessing --------------------------------------
filter_to_angiosperms <- function(df, species_names) {
  df[df$wfo_name %in% species_names, ]
}

filter_nonzero_colors <- function(df, cols) {
  df[rowSums(df[cols]) != 0, ]
}

Data_GPT        <- filter_nonzero_colors(filter_to_angiosperms(Data_GPT,        Species_angiosperms), color_cols)
Data_Biol       <- filter_nonzero_colors(filter_to_angiosperms(Data_Biol,       Species_angiosperms), color_cols)
Data_Orchid     <- filter_nonzero_colors(filter_to_angiosperms(Data_Orchid,     Species_angiosperms), color_cols)
Data_Urban      <- filter_nonzero_colors(filter_to_angiosperms(Data_Urban,      Species_angiosperms), color_cols)
Data_USDA       <- filter_nonzero_colors(filter_to_angiosperms(Data_USDA,       Species_angiosperms), color_cols)
Data_African    <- filter_nonzero_colors(filter_to_angiosperms(Data_African,    Species_angiosperms), color_cols)
Data_Functional <- filter_nonzero_colors(filter_to_angiosperms(Data_Functional, Species_angiosperms), color_cols)
Data_TRYALL     <- filter_nonzero_colors(filter_to_angiosperms(Data_TRYALL,     Species_angiosperms), color_cols)

# ---- Output 1: Consistency_results.csv -------------------------
consistency_df_list <- list(
  GPT4o       = Data_GPT,
  TRY_all     = Data_TRYALL,
  TRY_Biol    = Data_Biol,
  TRY_USDA    = Data_USDA,
  TRY_Orchid  = Data_Orchid,
  TRY_Urban   = Data_Urban,
  TRY_African = Data_African,
  TRY_Funct   = Data_Functional
)

consistency_combinations <- combn(names(consistency_df_list), 2, simplify = FALSE)

consistency_results_detailed <- map_dfr(consistency_combinations, ~ {
  df1_name <- .x[1]
  df2_name <- .x[2]
  df1 <- consistency_df_list[[df1_name]]
  df2 <- consistency_df_list[[df2_name]]

  comparison_df <- inner_join(df1, df2, by = "wfo_name", suffix = c(".df1", ".df2"))

  if (nrow(comparison_df) == 0) {
    return(tibble(
      comparison = paste(df1_name, "vs", df2_name),
      common_wfo_names = 0
    ))
  }

  analysis_df <- comparison_df %>%
    mutate(
      matched_color_count = rowSums(
        dplyr::select(., all_of(paste0(color_cols, ".df1"))) ==
          dplyr::select(., all_of(paste0(color_cols, ".df2"))),
        na.rm = TRUE
      ),
      similarity_per_wfo = (matched_color_count / length(color_cols)) * 100,
      is_full_match = (matched_color_count == length(color_cols))
    )

  summary_stats <- analysis_df %>%
    summarise(
      common_wfo_names = n(),
      full_match_rate = mean(is_full_match, na.rm = TRUE) * 100,
      avg_similarity = mean(similarity_per_wfo, na.rm = TRUE),
      sd_similarity = sd(similarity_per_wfo, na.rm = TRUE),
      min_similarity = min(similarity_per_wfo, na.rm = TRUE),
      max_similarity = max(similarity_per_wfo, na.rm = TRUE)
    )

  tibble(
    comparison = paste(df1_name, "vs", df2_name),
    !!!summary_stats
  )
})

write.csv(consistency_results_detailed, file.path(output_dir, "Consistency_results.csv"))

# ---- Output 2: Mismatch_report.csv -----------------------------
prepare_mismatch_df <- function(df, df_name, cols) {
  if (df_name == "GPT4o") {
    df %>%
      transmute(
        wfo_name,
        across(all_of(cols)),
        Proc_Description = Description,
        Proc_GPT_return1 = GPT_return1,
        Proc_GPT_return2 = GPT_return2
      )
  } else {
    df %>%
      transmute(
        wfo_name,
        across(all_of(cols)),
        Proc_Description = OrigValueStr,
        Proc_GPT_return1 = NA_character_,
        Proc_GPT_return2 = GPT_return2
      )
  }
}

mismatch_df_list <- list(
  GPT4o       = Data_GPT,
  TRY_Biol    = Data_Biol,
  TRY_USDA    = Data_USDA,
  TRY_Orchid  = Data_Orchid,
  TRY_Urban   = Data_Urban,
  TRY_African = Data_African,
  TRY_Funct   = Data_Functional
)

mismatch_combinations <- combn(names(mismatch_df_list), 2, simplify = FALSE)

mismatch_report <- map_dfr(mismatch_combinations, ~ {
  df1_name <- .x[1]
  df2_name <- .x[2]
  df1 <- mismatch_df_list[[df1_name]]
  df2 <- mismatch_df_list[[df2_name]]

  df1_proc <- prepare_mismatch_df(df1, df1_name, color_cols)
  df2_proc <- prepare_mismatch_df(df2, df2_name, color_cols)

  comparison_df <- inner_join(df1_proc, df2_proc, by = "wfo_name", suffix = c(".df1", ".df2"))
  if (nrow(comparison_df) == 0) return(NULL)

  mismatched_df <- comparison_df %>%
    mutate(
      matched_color_count = rowSums(
        dplyr::select(., all_of(paste0(color_cols, ".df1"))) ==
          dplyr::select(., all_of(paste0(color_cols, ".df2"))),
        na.rm = TRUE
      )
    ) %>%
    filter(matched_color_count != length(color_cols))

  if (nrow(mismatched_df) == 0) return(NULL)

  mismatched_df %>%
    rowwise() %>%
    mutate(
      df1_colors = paste0(ifelse(c_across(all_of(paste0(color_cols, ".df1"))), "T", "F"), collapse = ""),
      df2_colors = paste0(ifelse(c_across(all_of(paste0(color_cols, ".df2"))), "T", "F"), collapse = "")
    ) %>%
    ungroup() %>%
    transmute(
      df1_name = df1_name,
      df2_name = df2_name,
      wfo_name,
      df1_colors,
      df2_colors,
      df1_Description = Proc_Description.df1,
      df1_GPT_return1 = Proc_GPT_return1.df1,
      df1_GPT_return2 = Proc_GPT_return2.df1,
      df2_Description = Proc_Description.df2,
      df2_GPT_return1 = Proc_GPT_return1.df2,
      df2_GPT_return2 = Proc_GPT_return2.df2
    )
})

temp <- c()
for (i in mismatch_report$wfo_name) {
  temp <- c(temp, Data_all[Data_all$wfo_name == i, ]$Combine_status)
}

mismatch_report$Conbine_status <- temp
write.csv(mismatch_report, file.path(output_dir, "Mismatch_report.csv"))

# ---- Outputs 3-6: RR network figure and edge tables -----------
gpt_vs_try <- mismatch_report %>%
  filter(df1_name == "GPT4o" | df2_name == "GPT4o")

try_vs_try <- mismatch_report %>%
  filter(df1_name != "GPT4o" & df2_name != "GPT4o")

analyze_pairs_binary <- function(df, group_name,
                                 colors = c("White", "Yellow", "Red", "Blue", "Purple", "Green"),
                                 src_prefix = "df1_", tgt_prefix = "df2_",
                                 src_bits = "df1_colors", tgt_bits = "df2_colors") {
  have_cols <- all(paste0(src_prefix, colors) %in% names(df)) &&
    all(paste0(tgt_prefix, colors) %in% names(df))
  have_bits <- all(c(src_bits, tgt_bits) %in% names(df))

  if (!have_cols && !have_bits) {
    stop(
      "想定列が見つかりません。ワイド列: ",
      paste0(src_prefix, colors, collapse = ", "), " / ",
      paste0(tgt_prefix, colors, collapse = ", "),
      " もしくはビット列: ", src_bits, ", ", tgt_bits, " を用意してください。"
    )
  }

  bit_at <- function(x, k) {
    ch <- substring(as.character(x), k, k)
    as.integer(ch %in% c("1", "T", "t", "Y", "y"))
  }

  get_triplet <- function(i, j, direction) {
    if (direction == "df1_vs_df2") {
      src_pref <- src_prefix
      tgt_pref <- tgt_prefix
      source_name <- "df1"
      src_b <- src_bits
      tgt_b <- tgt_bits
    } else {
      src_pref <- tgt_prefix
      tgt_pref <- src_prefix
      source_name <- "df2"
      src_b <- tgt_bits
      tgt_b <- src_bits
    }

    if (have_cols) {
      si <- as.integer(df[[paste0(src_pref, colors[i])]])
      ti <- as.integer(df[[paste0(tgt_pref, colors[i])]])
      tj <- as.integer(df[[paste0(tgt_pref, colors[j])]])
    } else {
      si <- bit_at(df[[src_b]], i)
      ti <- bit_at(df[[tgt_b]], i)
      tj <- bit_at(df[[tgt_b]], j)
    }

    list(si = si, ti = ti, tj = tj, source_name = source_name)
  }

  res <- list()

  for (direction in c("df1_vs_df2", "df2_vs_df1")) {
    for (i in seq_along(colors)) {
      for (j in seq_along(colors)) {
        if (i == j) next

        trip <- get_triplet(i, j, direction)
        si <- trip$si
        ti <- trip$ti
        tj <- trip$tj

        A <- (ti == 0)
        B <- (si == 1)
        J <- (tj == 1)

        valid <- !is.na(A) & !is.na(B) & !is.na(J)
        if (!any(valid)) next

        a <- sum(valid &  B &  J)
        b <- sum(valid &  B & !J)
        c <- sum(valid & !B &  J)
        d <- sum(valid & !B & !J)

        n_treated <- a + b
        n_control <- c + d
        n_A <- n_treated + n_control
        if (n_treated == 0 || n_control == 0) next

        p <- a / n_treated
        p0_A <- (a + c) / n_A
        p0_c <- c / n_control

        lift_A <- if (p0_A == 0) NA_real_ else p / p0_A
        rr <- if (p0_c == 0) NA_real_ else p / p0_c

        m <- matrix(c(a, b, c, d), nrow = 2)
        ft <- suppressWarnings(fisher.test(m))
        or <- unname(ft$estimate)
        pval <- ft$p.value
        ci <- ft$conf.int

        res[[length(res) + 1]] <- tibble(
          comparison_group = group_name,
          direction = direction,
          source_name = trip$source_name,
          i = colors[i],
          j = colors[j],
          a = a,
          b = b,
          c = c,
          d = d,
          n_treated = n_treated,
          n_control = n_control,
          n_A = n_A,
          p = p,
          baseline_A = p0_A,
          baseline_control = p0_c,
          lift_A = lift_A,
          rr = rr,
          or = or,
          or_low = ci[1],
          or_high = ci[2],
          pval = pval
        )
      }
    }
  }

  out <- dplyr::bind_rows(res)
  if (nrow(out) == 0) return(out)
  out$qval <- p.adjust(out$pval, method = "BH")
  out
}

summarize_edges <- function(res, weight = c("rr", "or")) {
  weight <- match.arg(weight)

  res %>%
    filter(is.finite(.data[[weight]])) %>%
    group_by(i, j) %>%
    summarise(
      w = weighted.mean(.data[[weight]], n_treated, na.rm = TRUE),
      n_treated = sum(n_treated),
      n_control = sum(n_control),
      qmin = min(qval),
      .groups = "drop"
    ) %>%
    rename(from = i, to = j)
}

create_network_plot_rr <- function(edges, title, node_order) {
  edges_filtered <- edges %>%
    filter(n_treated >= 30, n_control >= 30, w >= 1.2, qmin < 0.1)

  all_nodes <- tibble(name = factor(node_order, levels = node_order))
  g <- graph_from_data_frame(edges_filtered, vertices = all_nodes, directed = TRUE)

  ggraph(g, layout = "circle") +
    geom_edge_fan(
      aes(width = w),
      arrow = grid::arrow(length = grid::unit(2, "mm")),
      end_cap = circle(3, "mm")
    ) +
    scale_edge_width(
      range = c(1, 4),
      limits = c(1, 4),
      name = "RR (mean)"
    ) +
    geom_node_point(size = 10, color = "skyblue") +
    geom_node_text(aes(label = name), fontface = "bold") +
    theme_graph(base_family = "sans") +
    labs(title = title) +
    coord_fixed()
}

res_gpt <- analyze_pairs_binary(
  gpt_vs_try,
  "GPT_vs_TRY",
  src_prefix = "df1_",
  tgt_prefix = "df2_"
)

res_try <- analyze_pairs_binary(
  try_vs_try,
  "TRY_vs_TRY",
  src_prefix = "df1_",
  tgt_prefix = "df2_"
)

edges_gpt_df1 <- summarize_edges(res_gpt %>% filter(source_name == "df1"), weight = "rr")
edges_gpt_df2 <- summarize_edges(res_gpt %>% filter(source_name == "df2"), weight = "rr")
edges_try     <- summarize_edges(res_try, weight = "rr")

plot_gpt1 <- create_network_plot_rr(edges_gpt_df1, "GPT adds TRY (RR)", node_order)
plot_gpt2 <- create_network_plot_rr(edges_gpt_df2, "TRY adds GPT (RR)", node_order)
plot_try  <- create_network_plot_rr(edges_try,     "TRY adds TRY (RR)", node_order)

plot_network <- gridExtra::grid.arrange(plot_gpt1, plot_gpt2, plot_try, ncol = 3)

ggsave(
  file = file.path(output_dir, "Data_network_RR.png"),
  plot = plot_network,
  dpi = 300,
  width = 15,
  height = 5
)

write.csv(edges_gpt_df1, file.path(output_dir, "Data_network_gptAddtrye_RR.csv"))
write.csv(edges_gpt_df2, file.path(output_dir, "Data_network_tryAddgpt_RR.csv"))
write.csv(edges_try,     file.path(output_dir, "Data_network_trytry_RR.csv"))
