#############################
# Unified Processing Script #
#############################

### 事前準備：必要なライブラリ・関数の読み込み
source("R_function.R")        # Lc(), Search_GBIF() などの関数を含む
source("GBIF_function.R")       # GBIF用の補助関数
library(terra)                # rast()等の関数のため
library(raster)               # WorldClimデータ取得用
library(taxize)               # 分類学情報取得用
library(FactoMineR)           # PCA解析用
library(stringr)              # 文字列処理用

#############################
### Step 1. GBIF検索      ###
#############################
# ※初期のデータはGPT4o_all_fin.csvとTRY_data_categ_fin.csvから読み込む

Data_dir <- "InputData/"
# 2つのファイルから読み込み
file1 <- Lc(Data_dir, "GPT4o_all_fin.csv")
file2 <- Lc(Data_dir, "TRY_data_categ_fin.csv")
data1 <- read.csv(file1, header = TRUE)
data2 <- read.csv(file2, header = TRUE)

# 両ファイルの"Species"列を統合し重複を除く
species_vec <- unique(c(data1$Species, data2$Species))
Species_df <- data.frame(Species = species_vec, stringsAsFactors = FALSE)
Species_df$GBIF_search <- "NotYet"

# GBIF結果出力ディレクトリ
GBIF_Out_dir <- "GBIF_Output/"

# 各種についてGBIFから出現データを取得
for(i in 1:nrow(Species_df)) {
  sp_name <- Species_df[i, "Species"]
  cat(i, " : ", sp_name, "\n")
  
  # すでにデータが存在していればスキップ
  if (file.exists(Lc(GBIF_Out_dir, sp_name, ".csv"))) {
    Species_df[i, "GBIF_search"] <- "Success"
    cat("  skip (既存ファイル)\n")
    next
  }
  
  # GBIFからデータ取得（Search_GBIF関数を利用）
  Sp_loc <- Search_GBIF(sp_name)
  if (mode(Sp_loc) == "list") {
    Species_df[i, "GBIF_search"] <- "Success"
    write.csv(Sp_loc, Lc(GBIF_Out_dir, sp_name, ".csv"))
  } else {
    Species_df[i, "GBIF_search"] <- Sp_loc  # 失敗理由等を記録
    next
  }
}
# GBIF検索結果一覧を保存
write.csv(Species_df, Lc(GBIF_Out_dir, "Processed_SpeciesData.csv"))

#############################
### Step 2. 環境データ付加  ###
#############################
# ※ここでは土壌パラメータなどのGeoTiff(.vrt)ファイルを読み込み，
#    各種GBIF出現データに環境情報を付加します。

# --- 土壌パラメータ（GeoTiff）の読み込み ---
GeoTiff_dir1 <- "Environmental_Data/Soil/"
GeoTiff1_list <- list.files(GeoTiff_dir1, pattern = ".vrt")
GeoTiff1 <- list()
for(i in 1:length(GeoTiff1_list)) {
  GeoTiff1[[i]] <- rast(Lc(GeoTiff_dir1, GeoTiff1_list[i]))
}

# --- WRB土壌分類データの読み込み ---
GeoTiff_dir2 <- "Environmental_Data/Soil/WRB/"
GeoTiff2_list <- list.files(GeoTiff_dir2, pattern = ".vrt")
GeoTiff2 <- list()
for(i in 1:length(GeoTiff2_list)) {
  GeoTiff2[[i]] <- rast(Lc(GeoTiff_dir2, GeoTiff2_list[i]))
}
cat("TIFFファイルをロードしました。\n")

# --- 各種GBIF出現データに土壌環境データを追加 ---
Loc_dir <- "GBIF_Processed_Data/"  # GBIF出現データ格納ディレクトリ
Loc_list <- list.files(Loc_dir, pattern = ".csv")
endnum <- length(Loc_list)
pb <- txtProgressBar(min = 0, max = endnum, style = 3)
startnum <- 1

for (filename in Loc_list) {
  gbif_data <- read.csv(Lc(Loc_dir, filename), header = TRUE, row.names = 1)
  
  # 各GeoTiffファイルから土壌環境データを抽出
  for (i in 1:length(GeoTiff1)) {
    gbif_data <- Get_soilenv(gbif_data, GeoTiff1[[i]])
  }
  for (i in 1:length(GeoTiff2)) {
    gbif_data <- Get_soilenv(gbif_data, GeoTiff2[[i]])
  }
  
  write.csv(gbif_data, Lc(Loc_dir, filename))
  setTxtProgressBar(pb, startnum)
  startnum <- startnum + 1
}

# --- 気候データ（WorldClim）の抽出 ---
Out_dir_WorldClim <- "GBIF_Processed_Data_WorldClim/"
Clim_alt <- raster::getData('worldclim', var = 'alt', res = 5)  # 標高
Clim_bio <- raster::getData('worldclim', var = 'bio', res = 5)  # 生物気候変数

Loc_list <- list.files(Loc_dir, pattern = ".csv")
endnum <- length(Loc_list)
pb <- txtProgressBar(min = 0, max = endnum, style = 3)
startnum <- 1

for (filename in Loc_list) {
  gbif_data <- read.csv(Lc(Loc_dir, filename), header = TRUE, row.names = 1)
  Sp_loc <- data.frame(lon = gbif_data$decimalLongitude,
                       lat = gbif_data$decimalLatitude)
  Clim_df <- Extract_clim(Sp_loc, Clim_bio)
  Alt_df <- Extract_clim(Sp_loc, Clim_alt)
  gbif_data <- cbind(gbif_data, Clim_df, Alt_df)
  write.csv(gbif_data, Lc(Out_dir_WorldClim, filename))
  setTxtProgressBar(pb, startnum)
  startnum <- startnum + 1
}

#############################
### Step 3. 分類学チェック ###
#############################
# ※ここではFoC（Flora of China）用データとTRY用データそれぞれに対して
#    taxizeパッケージを用いた分類学情報の取得と必要に応じた手動修正を行います。

#### ① FoCデータの分類学情報の追加 ####
Data_file <- "Data/FoC_Data.csv"
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Genus_list <- unique(Data$Genus)

# NCBIからtaxonomic IDを取得（eudicots）
id_list <- get_uid(Genus_list, ask = FALSE, division_filter = "eudicots")
id_list_result <- attributes(id_list)
Genus_found <- Genus_list[!is.na(id_list_result$uri)]
Genus_notfound <- Genus_list[is.na(id_list_result$uri)]

# 試行としてmonocotsもチェック
id_list2 <- get_uid(Genus_notfound, ask = FALSE, division_filter = "monocots")
id_list2_result <- attributes(id_list2)
Genus_monocots <- Genus_notfound[!is.na(id_list2_result$uri)]
Genus_notfound2 <- Genus_notfound[is.na(id_list2_result$uri)]

# 該当するFamilyリストを取得
Family_list <- unique(Data[Data$Genus %in% Genus_notfound2, "Family"])
id_list3 <- get_uid(Family_list, ask = FALSE)

# 各分類学階層を取得
class_1 <- classification(id_list, db = "ncbi")
class_2 <- classification(id_list2, db = "ncbi")
class_3 <- classification(id_list3, db = "ncbi")
rank_cols <- c("genus", "family", "order", "class", "phylum", "kingdom")

# 分類学情報抽出用関数
extract_taxonomy <- function(class_list, genus_list) {
  result_df <- data.frame(matrix(nrow = length(class_list), ncol = length(rank_cols)))
  rownames(result_df) <- genus_list
  for (i in seq_along(class_list)) {
    if (is.null(class_list[[i]])) next
    temp_df <- class_list[[i]]$taxa[class_list[[i]]$taxa$rank %in% rank_cols, ]
    temp_df$rank <- factor(temp_df$rank, levels = rank_cols)
    temp_df <- temp_df[order(temp_df$rank, decreasing = TRUE), ]
    if (nrow(temp_df) == length(rank_cols)) {
      result_df[i, ] <- temp_df$name
    }
  }
  colnames(result_df) <- rank_cols
  result_df <- na.omit(result_df)
  return(result_df)
}
Taxa_genus <- rbind(extract_taxonomy(class_1, Genus_found),
                    extract_taxonomy(class_2, Genus_monocots))
Taxa_family <- extract_taxonomy(class_3, Family_list)

# 必要に応じて手動修正（例：Acrocephalus）
Taxa_genus["Acrocephalus", ] <- c("Viridiplantae", "Streptophyta", "Magnoliopsida",
                                   "Lamiales", "Lamiaceae", "Acrocephalus")

# 分類学データの保存
write.csv(Taxa_genus, "Data/Taxa_Genus.csv")
write.csv(Taxa_family, "Data/Taxa_Family.csv")

# FoCデータへ分類学情報を統合
Taxa_genus$Genus <- Taxa_genus$genus
Taxa_family$Family <- Taxa_family$family
Data2 <- merge(Data, Taxa_genus, by = "Genus", all.x = TRUE)
Data2 <- merge(Data2, Taxa_family, by = "Family", all = TRUE)

# 複数の列に分かれる情報を統合する関数
consolidate_ranks <- function(df, rank_cols) {
  out_list <- vector("list", length(rank_cols))
  names(out_list) <- rank_cols
  for (rank in rank_cols) {
    out_list[[rank]] <- apply(df[, grep(rank, colnames(df)), drop = FALSE], 1, function(x) {
      x <- na.omit(x)
      if (length(x) == 0) return("NoData")
      return(x[1])
    })
  }
  return(data.frame(out_list))
}
rank_cols2 <- c("family", "order", "class", "phylum", "kingdom")
out_df_tax <- consolidate_ranks(Data2, rank_cols2)

# 必要な出力列を選択し最終データを作成
out_colnames <- c("Species", "White", "Yellow", "Red", "Blue", "Purple",
                  "Green", "NoDescription", "GBIF_data", "GBIF_datacount")
Final_Taxonomy_Data <- cbind(Data2[out_colnames], out_df_tax)
Final_Taxonomy_Data <- Final_Taxonomy_Data[!duplicated(Final_Taxonomy_Data), ]
write.csv(Final_Taxonomy_Data, "Data/Final_Taxonomy_Data.csv")

#### ② TRYデータの分類学情報の追加 ####
Data_file <- "Data/TRY_data.csv"
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Genus_list <- sapply(strsplit(Data$AccSpeciesName, " "), `[`, 1)
Genus_list_unique <- unique(Genus_list)

# NCBIからID取得
id_list <- get_uid(Genus_list_unique, ask = FALSE)
id_list_result <- attributes(id_list)
Genus_found <- Genus_list_unique[!is.na(id_list_result$uri)]
Genus_notfound <- Genus_list_unique[is.na(id_list_result$uri)]
class_1 <- classification(id_list, db = "ncbi")

# 同じ関数を利用して分類学情報抽出
rank_cols <- c("genus", "family", "order", "class", "phylum", "kingdom")
out_df1 <- extract_taxonomy(class_1, Genus_found)

# 手動で補完すべき分類群をデータフレームで作成（例）
manual_df <- data.frame(
  Atadinus = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Rosales", "Rhamnaceae", "Atadinus"),
  Schlagintweitia = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Schlagintweitia"),
  Blepharoneuron = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Poales", "Poaceae", "Blepharoneuron"),
  Endotropis = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Gentianales", "Apocynaceae", "Endotropis"),
  Adenorandia = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Gentianales", "Rubiaceae", "Adenorandia"),
  Ageratinastrum = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Ageratinastrum"),
  Afroaster = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Afroaster"),
  Byrsanthus = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Malpighiales", "Salicaceae", "Byrsanthus"),
  Congolanthus = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Gentianales", "Gentianaceae", "Congolanthus"),
  Dewildemania = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Dewildemania"),
  Dobera = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Brassicales", "Salvadoraceae", "Dobera"),
  Gerardiina = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Lamiales", "Orobanchaceae", "Gerardiina"),
  Hiernia = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Lamiales", "Orobanchaceae", "Hiernia"),
  Lepidostephium = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Lepidostephium"),
  Litogyne = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Litogyne"),
  Mariscus = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Poales", "Cyperaceae", "Mariscus"),
  Physotrichia = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Apiales", "Apiaceae", "Physotrichia"),
  Sphaerocodon = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Gentianales", "Apocynaceae", "Sphaerocodon"),
  Stenandriopsis = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Lamiales", "Acanthaceae", "Stenandriopsis"),
  Triceratorhynchus = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asparagales", "Orchidaceae", "Triceratorhynchus"),
  Uragoga = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Gentianales", "Rubiaceae", "Uragoga"),
  Pseudopegolettia = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Pseudopegolettia"),
  Oocephala = c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae", "Oocephala")
)
manual_df <- data.frame(t(manual_df))
Taxa_genus_TRY <- rbind(manual_df, out_df1)
# ※必要に応じて列名や順序を調整してください
write.csv(Taxa_genus_TRY, "Data/Taxa_TRY_Genus.csv")

Taxa_genus_TRY$Genus <- Taxa_genus_TRY$genus
Data$Genus <- Genus_list
Data2_TRY <- merge(Data, Taxa_genus_TRY, by = "Genus", all.x = TRUE)
write.csv(Data2_TRY, "Data/TRY_Final_Processed.csv", row.names = FALSE)

#############################
### Step 4. TRYデータとFoCデータの統合 ###
#############################
# ※ここでは両データセットを読み込み、種名の正規化や各種カラムの整備後に統合します

foc_file <- "Data/FoC_fin_20240723.csv"
try_file <- "Data/TRY_fin_20240729.csv"
merged_file <- "Data/Data_merged_FoC_TRY_20240730.csv"

# データ読み込み
Data1 <- read.csv(foc_file, header = TRUE, row.names = 1)
Data2 <- read.csv(try_file, header = TRUE, row.names = 1)

# Data1：種名からgenus抽出，種名の大文字小文字統一
Data1$genus <- sapply(strsplit(Data1$Species, " "), `[`, 1)
capitalize_species <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}
Data1$Species <- sapply(Data1$Species, capitalize_species)
Data2$Species <- sapply(Data2$Species, capitalize_species)

# データセット識別子の付加
Data1$Datasets <- "FoC"
Data2$Datasets <- "TRY"
write.csv(Data1, "Data/FoC_fin_20240730.csv", row.names = FALSE)
write.csv(Data2, "Data/TRY_fin_20240730.csv", row.names = FALSE)

# 両データセットの結合
Data_merge <- rbind(Data1, Data2)
color_cols <- c("White", "Yellow", "Red", "Blue", "Purple", "Green", "NoDescription")
rank_cols <- c("genus", "family", "order", "class", "phylum", "kingdom",
               "GBIF_data", "GBIF_datacount")

# 重複種リストの抽出（複数データセットに現れる種）
sp_count <- table(Data_merge$Species)
multisp_list <- names(sp_count[sp_count > 1])

# 重複種について色データを集約（存在するかを二値化）
out_df_color <- data.frame(matrix(nrow = 0, ncol = length(color_cols)))
for (species in multisp_list) {
  temp_df <- Data_merge[Data_merge$Species == species, color_cols]
  temp_df <- apply(temp_df, 2, function(x) as.numeric(sum(x) > 0))
  out_df_color <- rbind(out_df_color, temp_df)
}
colnames(out_df_color) <- color_cols
out_df_color$Species <- multisp_list
out_df_color <- out_df_color[order(out_df_color$Species), ]

# 重複種のうち1件のみを代表として処理
Data_merge_multi <- Data_merge[Data_merge$Species %in% multisp_list, ]
Data_merge_multi <- Data_merge_multi[!duplicated(Data_merge_multi$Species), ]
Data_merge_multi <- Data_merge_multi[order(Data_merge_multi$Species), ]
Data_merge_multi$FoC <- 1
Data_merge_multi$TRY <- 1
Data_merge_multi <- cbind(out_df_color, Data_merge_multi[rank_cols])

# 単独種の処理
Data_merge_single <- Data_merge[!Data_merge$Species %in% multisp_list, ]
Data_merge_single$FoC <- 0
Data_merge_single$TRY <- 0
Data_merge_single$FoC[Data_merge_single$Species %in% Data1$Species] <- 1
Data_merge_single$TRY[Data_merge_single$Species %in% Data2$Species] <- 1
Data_merge_single <- Data_merge_single[, colnames(Data_merge_single) != "Datasets"]

# 単独種と複数種のデータを統合
Data_merged <- rbind(Data_merge_single, Data_merge_multi)
Data_merged <- Data_merged[order(Data_merged$Species), ]
write.csv(Data_merged, merged_file, row.names = FALSE)

# --- 欠損している分類学情報の補完 ---
Data_merged <- read.csv(merged_file, header = TRUE, row.names = 1)
Notaxa_family <- unique(Data_merged[Data_merged$kingdom == "NoData", "family"])
for (family in Notaxa_family) {
  temp_df <- Data_merged[Data_merged$family == family, ]
  if (nrow(temp_df) == 1) {
    print(paste("Manual review:", family))
    next
  }
  temp_df_valid <- temp_df[temp_df$kingdom != "NoData", rank_cols]
  if (nrow(temp_df_valid) == 0) {
    print(paste("No valid taxonomy:", family))
    next
  }
  temp_df_valid <- temp_df_valid[1, rank_cols]
  for (rank in rank_cols) {
    Data_merged[Data_merged$family == family, rank] <- temp_df_valid[rank]
  }
}
write.csv(Data_merged, merged_file, row.names = FALSE)
Data_merged <- read.csv(merged_file, header = TRUE, row.names = 1)

#############################
### Step 5. 環境データ前処理  ###
#############################
# ※GBIF出現データに基づき，各種環境・土壌パラメータの集約値を計算します

# 入力種リスト（種ごとに1行）
Data <- read.csv("InputData/MergedSpeciesData.csv", header = TRUE, row.names = 1)
Loc_dir <- "GBIF_Processed_Data/"
Loc_list <- list.files(Loc_dir, pattern = ".csv")

# 使用する環境変数名（WorldClim等）
env_cols <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
              "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15",
              "bio16", "bio17", "bio18", "bio19", "Alt_df")
soil_cols <- c("bdod_0.5cm_mean", "cec_0.5cm_mean", "cfvo_0.5cm_mean",
               "clay_0.5cm_mean", "nitrogen_0.5cm_mean", "ocd_0.5cm_mean",
               "ocs_0.30cm_mean", "phh2o_0.5cm_mean", "sand_0.5cm_mean",
               "silt_0.5cm_mean", "soc_0.5cm_mean")

# 各種について，該当するGBIF環境データファイルが存在するかチェック
file_exists <- sapply(Data$Species, function(sp) {
  file.exists(Lc(Loc_dir, paste0(sp, ".csv")))
})
GBIF_sp <- Data$Species[file_exists]

# 各種ごとに環境情報の要約（20%トリム平均など）を計算
out_list <- list()
for (sp in GBIF_sp) {
  sp_file <- Lc(Loc_dir, paste0(sp, ".csv"))
  env_df <- read.csv(sp_file, row.names = 1)
  
  # 対象列のみ抽出，NAと重複行を除去
  env_df <- na.omit(env_df[c(env_cols, soil_cols)])
  env_df <- env_df[!duplicated(env_df), ]
  
  # 出現点が6点未満ならスキップ
  if (nrow(env_df) < 6) next
  
  env_list <- apply(env_df, 2, function(x) mean(x, trim = 0.2))
  # 土壌の環境変動（ユークリッド距離の平均）を追加
  env_list["Soil_dist"] <- mean(dist(env_df[soil_cols], method = "euclidean"))
  env_list["GBIF_accept"] <- nrow(env_df)
  env_list["Species"] <- sp
  
  out_list[[sp]] <- env_list
}
out_df_env <- data.frame(t(data.frame(out_list)))
# 元データと結合して保存
out_df_env <- merge(Data, out_df_env, by = "Species", all = TRUE)
write.csv(out_df_env, "OutputData/MergedSpeciesData_withEnvironment.csv")

#############################
### Step 6. PCAによる環境軸作成 ###
#############################
# ※必要に応じた外部関数を読み込み，環境変数（気候+土壌）を標準化後にPCAを実施します

source("Function_source.R")  # 必要な関数群の読み込み

data_file <- "Data/Data_merged_Angiosperms_withEnvironment.csv"
output_csv <- "Data/Data_for_Analysis.csv"
output_pc_df <- "Output/PC_df"
output_pca_rda <- "Output/Variable_PCA.rda"
output_pca_coord <- "Output/PCA_coordinates.csv"

Data <- read.csv(data_file, header = TRUE, row.names = 1)
# GBIF出現点数が存在する種のみを対象
Data <- Data[!is.na(Data$GBIF_accept), ]

# 環境変数のリスト
env_cols <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
              "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15",
              "bio16", "bio17", "bio18", "bio19", "Alt_df")
soil_cols <- c("bdod_0.5cm_mean", "cec_0.5cm_mean", "cfvo_0.5cm_mean",
               "clay_0.5cm_mean", "nitrogen_0.5cm_mean", "ocd_0.5cm_mean",
               "ocs_0.30cm_mean", "phh2o_0.5cm_mean", "sand_0.5cm_mean",
               "silt_0.5cm_mean", "soc_0.5cm_mean", "Soil_dist")
variable_list <- c(env_cols, soil_cols)

# 標準化してPCA実施（10主成分まで保持）
Data[variable_list] <- scale(Data[variable_list])
Variable_PCA <- PCA(Data[, variable_list], scale.unit = TRUE, ncp = 10, graph = FALSE)
print(Variable_PCA$eig)      # 各PCの固有値（分散説明率）
print(Variable_PCA$var$coord)  # 変数のPC空間上での座標

# 各サンプルのPC得点をデータフレームにまとめる
PC_names <- paste0("PC", 1:10)
PC_df <- data.frame(Variable_PCA$ind$coord[, 1:10])
colnames(PC_df) <- PC_names
PC_df$Species <- Data$Species

# 元データとPC得点を結合して保存
Data <- merge(Data, PC_df, by = "Species", all = TRUE)
write.csv(Data, output_csv, row.names = FALSE)
save(PC_df, file = output_pc_df)
save(Variable_PCA, file = output_pca_rda)
write.csv(Variable_PCA$var$coord, output_pca_coord, row.names = TRUE)

#############################
# 処理完了
#############################
cat("全ての処理が完了しました。\n")
