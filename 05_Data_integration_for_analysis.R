# ============================================================
# test4_readable_same_output.R
# ------------------------------------------------------------
# Purpose:
#   Reorganized, easier-to-read version of test4.
#
# Important:
#   This script is intentionally kept output-compatible with test3.
#   To preserve strict output equivalence, potentially awkward behaviors
#   from test3/test4 are also retained as-is.
# ============================================================


# ============================================================
# 1) Load input datasets
# ============================================================
Data_GPT        <- read.csv("Data3/GPT4o_fin.csv",          header = T, row.names = 1)
Data_Biol       <- read.csv("Data3/TRY_Biol_fin.csv",       header = T, row.names = 1)
Data_Orchid     <- read.csv("Data3/TRY_Orchid_fin.csv",     header = T, row.names = 1)
Data_Urban      <- read.csv("Data3/TRY_Urban_fin.csv",      header = T, row.names = 1)
Data_USDA       <- read.csv("Data3/TRY_USDA_fin.csv",       header = T, row.names = 1)
Data_African    <- read.csv("Data3/TRY_Africa_fin.csv",     header = T, row.names = 1)
Data_Functional <- read.csv("Data3/TRY_Functional_fin.csv", header = T, row.names = 1)


# ============================================================
# 2) Shared column settings
# ============================================================
extract_columns <- c(
  "wfo_name", "taxonID", "family", "genus", "Group", "Hybrid",
  "Combine_status", "White", "Yellow", "Red", "Blue", "Purple",
  "Green", "NoDescription"
)

color_cols <- c("White", "Yellow", "Red", "Blue", "Purple", "Green")
meta_cols  <- setdiff(extract_columns, c("wfo_name", color_cols, "Hybrid", "Combine_status"))


# ============================================================
# 3) Helper: collapse Combine_status within a species
# ------------------------------------------------------------
# Rule retained exactly from test3/test4:
#   - If all values are one_to_one -> one_to_one
#   - Otherwise prefer many_to_one, then one_to_many
# ============================================================
pick_combine <- function(x) {
  x <- na.omit(as.character(x))
  if (length(x) == 0) return(NA_character_)
  if (any(x != "one_to_one")) {
    if ("many_to_one" %in% x) return("many_to_one")
    if ("one_to_many" %in% x) return("one_to_many")
  }
  "one_to_one"
}


# ============================================================
# 4) Integrate GPT + TRY datasets
# ============================================================
df_list <- list(
  GPT     = Data_GPT,
  Biol    = Data_Biol,
  Orchid  = Data_Orchid,
  Urban   = Data_Urban,
  USDA    = Data_USDA,
  African = Data_African,
  Funct   = Data_Functional
)

# Keep only required columns from each dataset.
df_list_selected <- purrr::map(df_list, ~ {
  dplyr::select(.x, any_of(extract_columns))
})

# Stack all datasets into a long table.
all_data_long <- bind_rows(df_list_selected, .id = "source")

# Summarize by wfo_name.
integrated_data <- all_data_long %>%
  filter(!is.na(wfo_name) & wfo_name != "") %>%
  group_by(wfo_name) %>%
  summarise(
    # Any dataset marking a color as present makes the integrated value 1.
    across(all_of(color_cols), ~ as.integer(any(.x == 1, na.rm = TRUE))),

    # If any source marks the species as Hybrid, keep Hybrid.
    Hybrid = if (any(Hybrid == "Hybrid", na.rm = TRUE)) "Hybrid" else "NotHybrid",

    # Resolve Combine_status with the helper above.
    Combine_status = pick_combine(Combine_status),

    # For metadata columns, keep the first non-NA value.
    across(all_of(meta_cols), ~ first(na.omit(.))),

    .groups = "drop"
  )

# Track whether each integrated species appears in GPT / TRY sources.
GPT_species <- Data_GPT$wfo_name
TRY_species <- unique(c(
  Data_Biol$wfo_name,
  Data_Orchid$wfo_name,
  Data_Urban$wfo_name,
  Data_USDA$wfo_name,
  Data_African$wfo_name,
  Data_Functional$wfo_name
))

integrated_data$GPT_data <- 0
integrated_data[integrated_data$wfo_name %in% GPT_species, ]$GPT_data <- 1
integrated_data$TRY_data <- 0
integrated_data[integrated_data$wfo_name %in% TRY_species, ]$TRY_data <- 1

# Detect hybrid names from spec.name when available.
hybrid_wfo_names_list <- purrr::map(df_list, ~ {
  if ("spec.name" %in% colnames(.x)) {
    .x %>%
      filter(grepl("_x|×", spec.name)) %>%
      pull(wfo_name)
  } else {
    character(0)
  }
})

hybrid_wfo_vector <- hybrid_wfo_names_list %>%
  unlist() %>%
  unique()

integrated_data[integrated_data$wfo_name %in% hybrid_wfo_vector, ]$Hybrid <- "Hybrid"
# integrated_data[!integrated_data$wfo_name %in% hybrid_wfo_vector, ]$Hybrid <- "NotHybrid"

# Important: row.names is intentionally not specified,
# because test3 writes with the default behavior.
write.csv(integrated_data, "Data3/Integrated_data.csv")


# ============================================================
# 5) Integrate TRY-only datasets
# ============================================================
df_list <- list(
  Biol    = Data_Biol,
  Orchid  = Data_Orchid,
  Urban   = Data_Urban,
  USDA    = Data_USDA,
  African = Data_African,
  Funct   = Data_Functional
)

# Keep only required columns from each TRY dataset.
df_list_selected <- purrr::map(df_list, ~ {
  dplyr::select(.x, any_of(extract_columns))
})

# Stack all TRY datasets into a long table.
all_data_long <- bind_rows(df_list_selected, .id = "source")

# Summarize by wfo_name.
integrated_data <- all_data_long %>%
  filter(!is.na(wfo_name) & wfo_name != "") %>%
  group_by(wfo_name) %>%
  summarise(
    # Any dataset marking a color as present makes the integrated value 1.
    across(all_of(color_cols), ~ as.integer(any(.x == 1, na.rm = TRUE))),

    # For metadata columns, keep the first non-NA value.
    across(all_of(meta_cols), ~ first(na.omit(.x))),

    .groups = 'drop'
  )

# Important: row.names is intentionally not specified,
# because test3 writes with the default behavior.
write.csv(integrated_data, "Data3/TRY_combined.csv")


# ============================================================
# 6) Merge integrated Group labels into previous master table
# ------------------------------------------------------------
# The explicit for-loop is retained from test3/test4 to preserve behavior.
# ============================================================
Prev_data  <- fread("Data2/All_data.csv",        data.table = F)
Integ_data <- fread("Data3/Integrated_data.csv", data.table = F)

group_list <- c()
for (i in Prev_data$wfo_name) {
  group_temp <- Integ_data[Integ_data$wfo_name == i, ]$Group
  group_list <- c(group_list, group_temp)
}

Prev_data$Group <- group_list
write.csv(Prev_data, "Data3/All_data.csv")  # Table S2 all data


# ============================================================
# 7) Basic summaries
# ============================================================
Data <- fread("Data3/All_data.csv", data.table = F)

nrow(Data)                              # All entry 26905
table(Data$Group)                       # Angiosperms 26766, OtherOrNotInWfo 139
table(Data[Data$GPT_data, ]$Group)      # Angiosperms 17098, Other 4
table(Data[Data$TRY_data, ]$Group)      # Angiosperms 10756, Other 136

Data_angiosperms <- Data[Data$Group == "Angiosperms", ]
nrow(Data_angiosperms)

table(Data_angiosperms[Data_angiosperms$GPT_data, ]$Hybrid)
table(Data_angiosperms[Data_angiosperms$TRY_data, ]$Hybrid)
sum(Data_angiosperms$Hybrid == "NotHybrid")


# ============================================================
# 8) Filter analysis dataset and save outputs
# ============================================================
Data <- fread("Data3/All_data.csv", data.table = F)
color_cols <- c("White", "Yellow", "Red", "Purple", "Blue", "Green")

Data <- Data[rowSums(Data[color_cols]) > 0, ]

Data[is.na(Data$bio1), ]$GBIF_counts_Env <- 0
Data[is.na(Data$bio1), ]$GBIF_counts     <- 0

Data <- Data[(Data$Group == "Angiosperms") &
               (Data$Hybrid == "NotHybrid") &
               (Data$Combine_status == "one_to_one") &
               (Data$GBIF_counts_Env > 5) &
               (rowSums(Data[color_cols]) > 0), ]

nrow(Data)
sum(Data$GPT_data)
sum(Data$TRY_data)

# Important: row.names is intentionally not specified,
# because test3 writes with the default behavior.
write.csv(Data, "Data3/Data_for_withoutPhylogeny.csv")

Data_withPhylo <- Data[Data$Phylogeny == T, ]

nrow(Data_withPhylo)
sum(Data_withPhylo$GPT_data)
sum(Data_withPhylo$TRY_data)

# Important: row.names is intentionally not specified,
# because test3 writes with the default behavior.
write.csv(Data_withPhylo, "Data3/Data_for_Analysis.csv")


# ============================================================
# 9) Build phylogeny input and save tree object
# ------------------------------------------------------------
# NOTE:
#   family = Data$family is retained exactly from test3/test4 to preserve
#   identical behavior/output, even though Data_withPhylo$family might look
#   more natural at first glance.
# ============================================================
Data_withPhylo <- fread("Data3/Data_for_Analysis.csv", data.table = F)

sp <- data.frame(
  species = Data_withPhylo$wfo_name,
  genus   = Data_withPhylo$genus,
  family  = Data$family,
  stringsAsFactors = FALSE
)

res_lcvp <- phylo.maker(
  sp.list   = sp,
  tree      = GBOTB.extended.LCVP,
  nodes     = nodes.info.1.LCVP,
  scenarios = c("S1")
)

save(res_lcvp, file = "Data3/Phylogeny_all_S1.rda")


# ============================================================
# 10) Plot basic color-overlap information
# ============================================================
Data <- fread("Data3/All_data.csv", data.table = F)
color_cols <- c("White", "Yellow", "Red", "Purple", "Blue", "Green")

upset_plot <- upset(
  Data,
  sets            = color_cols,
  mainbar.y.label = "Number of Species (Intersection Size > 6)",
  sets.x.label    = "Total Species per Color",
  order.by        = "freq",  # Order the intersections by frequency
  nsets           = 6,        # Number of sets to include
  nintersects     = 64,
  number.angles   = 30,
  point.size      = 3.5,
  line.size       = 2,
  text.scale      = c(1.3, 1.3, 1, 1, 2, 1.2)
)

png(file = "Output_20260120/Data_Color_Intersections.png", width = 12, height = 7, units = 'in', res = 300)
upset_plot
dev.off()

pdf(file = "Output_20260120/Data_Color_Intersections.pdf", width = 12, height = 7)
upset_plot
dev.off()


# ============================================================
# 11) Calculate and plot lift values
# ============================================================
library(pheatmap)

color_matrix         <- as.matrix(Data[color_cols])
total_species        <- nrow(color_matrix)
color_totals         <- colSums(color_matrix)
co_occurrence_matrix <- t(color_matrix) %*% color_matrix

lift_matrix <- matrix(
  0,
  nrow = length(color_cols),
  ncol = length(color_cols),
  dimnames = list(color_cols, color_cols)
)

for (i in 1:length(color_cols)) {
  for (j in 1:length(color_cols)) {
    if (i == j) {
      lift_matrix[i, j] <- NA  # Do not calculate the diagonal.
    } else {
      # Lift(A, B) = P(A and B) / (P(A) * P(B))
      #            = Count(A and B) * N / (Count(A) * Count(B))
      prob_A_and_B <- co_occurrence_matrix[i, j] / total_species
      prob_A       <- color_totals[i] / total_species
      prob_B       <- color_totals[j] / total_species

      # Avoid division by zero.
      if (prob_A == 0 | prob_B == 0) {
        lift_matrix[i, j] <- 0
      } else {
        lift_matrix[i, j] <- prob_A_and_B / (prob_A * prob_B)
      }
    }
  }
}

print(round(lift_matrix, 2))

pdf_filename <- "Output_20260120/lift_heatmap.pdf"
pdf(pdf_filename, width = 8, height = 7)
pheatmap(
  lift_matrix,
  display_numbers = TRUE,
  number_format   = "%.2f",
  main            = "Co-occurrence Lift between Colors",
  fontsize_number = 12,
  cluster_rows    = F,
  cluster_cols    = F,
  color           = colorRampPalette(c("lightblue", "ivory", "coral"))(100),
  na_col          = "grey70"
)
dev.off()

average_lift <- mean(lift_matrix, na.rm = TRUE)
