#############################################
# Unified Processing Script
# Comments have been added by ChatGPTo3-mini-high
#############################################

# Load required helper functions and libraries
source("R_function.R")        # Contains helper functions (e.g., Lc(), etc.)
source("GBIF_function.R")       # Contains functions for GBIF data retrieval
library(terra)                # For handling GeoTiff raster files
library(raster)               # For downloading WorldClim data
library(taxize)               # For taxonomic information retrieval
library(FactoMineR)           # For PCA analysis
library(stringr)              # For string manipulation

#############################################
# Step 0: Prepare Unique Species List
#############################################
# Read the two initial CSV files containing species data
Data_dir <- "InputData/"
file1 <- Lc(Data_dir, "GPT4o_all_fin.csv")
file2 <- Lc(Data_dir, "TRY_data_categ_fin.csv")
data1 <- read.csv(file1, header = TRUE)
data2 <- read.csv(file2, header = TRUE)

# Extract the "Species" column from both files and combine them
species_list <- unique(c(data1$Species, data2$Species))
Species_df <- data.frame(Species = species_list, stringsAsFactors = FALSE)
Species_df$GBIF_search <- "NotYet"  # Initialize GBIF search status

#############################################
# Step 1: Search GBIF Occurrence Information
#############################################
# Define output directory for GBIF results
GBIF_Out_dir <- "GBIF_Output/"

# Loop over each unique species and search for occurrence data via GBIF
for (i in 1:nrow(Species_df)) {
  sp_name <- Species_df[i, "Species"]
  cat(i, " : ", sp_name, "\n")
  
  # Check if a GBIF result file already exists; if so, skip the search
  if (file.exists(Lc(GBIF_Out_dir, sp_name, ".csv"))) {
    Species_df[i, "GBIF_search"] <- "Success"
    cat("  Skip: File exists\n")
    next
  }
  
  # Retrieve occurrence data from GBIF using the helper function
  Sp_loc <- Search_GBIF(sp_name)
  
  # If data retrieval is successful, save the results; otherwise, record the error
  if (mode(Sp_loc) == "list") {
    Species_df[i, "GBIF_search"] <- "Success"
    write.csv(Sp_loc, Lc(GBIF_Out_dir, sp_name, ".csv"))
  } else {
    Species_df[i, "GBIF_search"] <- Sp_loc  # Record error message
    next
  }
}
# Save the updated species list with GBIF search status
write.csv(Species_df, Lc(GBIF_Out_dir, "Processed_SpeciesData.csv"))

#############################################
# Step 2: Adding Environmental Information to Occurrence Data
#############################################
# --- Load Environmental Raster Data ---

# Load soil parameter GeoTiff files (.vrt)
GeoTiff_dir1 <- "Environmental_Data/Soil/"
GeoTiff1_list <- list.files(GeoTiff_dir1, pattern = ".vrt")
GeoTiff1 <- list()
for (i in 1:length(GeoTiff1_list)) {
  GeoTiff1[[i]] <- rast(Lc(GeoTiff_dir1, GeoTiff1_list[i]))
}

# Load WRB soil classification GeoTiff files (.vrt)
GeoTiff_dir2 <- "Environmental_Data/Soil/WRB/"
GeoTiff2_list <- list.files(GeoTiff_dir2, pattern = ".vrt")
GeoTiff2 <- list()
for (i in 1:length(GeoTiff2_list)) {
  GeoTiff2[[i]] <- rast(Lc(GeoTiff_dir2, GeoTiff2_list[i]))
}
cat("TIFF files have been loaded.\n")

# --- Process Each GBIF Occurrence File ---

# Directory containing species-specific GBIF occurrence data
Loc_dir <- "GBIF_Processed_Data/"
Loc_list <- list.files(Loc_dir, pattern = ".csv")
endnum <- length(Loc_list)
pb <- txtProgressBar(min = 0, max = endnum, style = 3)
startnum <- 1

# Loop over each occurrence file to add soil environmental data
for (filename in Loc_list) {
  gbif_data <- read.csv(Lc(Loc_dir, filename), header = TRUE, row.names = 1)
  
  # For each soil parameter raster, extract the corresponding environmental data
  for (i in 1:length(GeoTiff1)) {
    gbif_data <- Get_soilenv(gbif_data, GeoTiff1[[i]])
  }
  # Similarly, extract WRB soil classification data
  for (i in 1:length(GeoTiff2)) {
    gbif_data <- Get_soilenv(gbif_data, GeoTiff2[[i]])
  }
  
  # Save the updated occurrence data
  write.csv(gbif_data, Lc(Loc_dir, filename))
  setTxtProgressBar(pb, startnum)
  startnum <- startnum + 1
}

# --- Climate Data Extraction ---
# Define output directory for climate-enhanced data
Out_dir_WorldClim <- "GBIF_Processed_Data_WorldClim/"

# Download climate data from WorldClim (using a 5 arc-minute resolution here)
Clim_alt <- raster::getData('worldclim', var = 'alt', res = 5)  # Elevation data
Clim_bio <- raster::getData('worldclim', var = 'bio', res = 5)  # Bioclimatic variables

# Process each occurrence file for climate data extraction
Loc_list <- list.files(Loc_dir, pattern = ".csv")
endnum <- length(Loc_list)
pb <- txtProgressBar(min = 0, max = endnum, style = 3)
startnum <- 1

for (filename in Loc_list) {
  gbif_data <- read.csv(Lc(Loc_dir, filename), header = TRUE, row.names = 1)
  
  # Create a data frame of coordinates from occurrence data
  Sp_loc <- data.frame(lon = gbif_data$decimalLongitude,
                       lat = gbif_data$decimalLatitude)
  
  # Extract climate and elevation data based on coordinates
  Clim_df <- Extract_clim(Sp_loc, Clim_bio)
  Alt_df <- Extract_clim(Sp_loc, Clim_alt)
  
  # Combine the new climate data with the occurrence data
  gbif_data <- cbind(gbif_data, Clim_df, Alt_df)
  write.csv(gbif_data, Lc(Out_dir_WorldClim, filename))
  setTxtProgressBar(pb, startnum)
  startnum <- startnum + 1
}

#############################################
# Step 3: Taxonomy Check with taxize and Manual Editing
#############################################
# --- Part A: Flora of China (FoC) Data ---
Data_file <- "Data/FoC_Data.csv"
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Genus_list <- unique(Data$Genus)

# Fetch taxonomic IDs for eudicots from NCBI
id_list <- get_uid(Genus_list, ask = FALSE, division_filter = "eudicots")
id_list_result <- attributes(id_list)
Genus_found <- Genus_list[!is.na(id_list_result$uri)]
Genus_notfound <- Genus_list[is.na(id_list_result$uri)]

# Attempt to fetch IDs for monocots from the not found list
id_list2 <- get_uid(Genus_notfound, ask = FALSE, division_filter = "monocots")
id_list2_result <- attributes(id_list2)
Genus_monocots <- Genus_notfound[!is.na(id_list2_result$uri)]
Genus_notfound2 <- Genus_notfound[is.na(id_list2_result$uri)]

# Fetch taxonomic IDs for families for remaining genera
Family_list <- unique(Data[Data$Genus %in% Genus_notfound2, "Family"])
id_list3 <- get_uid(Family_list, ask = FALSE)

# Retrieve classification hierarchy from NCBI
class_1 <- classification(id_list, db = "ncbi")
class_2 <- classification(id_list2, db = "ncbi")
class_3 <- classification(id_list3, db = "ncbi")
rank_cols <- c("genus", "family", "order", "class", "phylum", "kingdom")

# Function to extract taxonomy information from classification data
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

# Manually correct taxonomy if needed (example provided)
Taxa_genus["Acrocephalus", ] <- c("Viridiplantae", "Streptophyta", "Magnoliopsida",
                                   "Lamiales", "Lamiaceae", "Acrocephalus")

# Save the taxonomic data for FoC
write.csv(Taxa_genus, "Data/Taxa_Genus.csv")
write.csv(Taxa_family, "Data/Taxa_Family.csv")

# Merge the taxonomic information with the FoC species data
Taxa_genus$Genus <- Taxa_genus$genus
Taxa_family$Family <- Taxa_family$family
Data2 <- merge(Data, Taxa_genus, by = "Genus", all.x = TRUE)
Data2 <- merge(Data2, Taxa_family, by = "Family", all = TRUE)

# Function to consolidate taxonomic ranks across multiple columns
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

# Define output columns and create the final FoC dataset
out_colnames <- c("Species", "White", "Yellow", "Red", "Blue", "Purple", "Green",
                  "NoDescription", "GBIF_data", "GBIF_datacount")
Final_Taxonomy_Data <- cbind(Data2[out_colnames], out_df_tax)
Final_Taxonomy_Data <- Final_Taxonomy_Data[!duplicated(Final_Taxonomy_Data), ]
write.csv(Final_Taxonomy_Data, "Data/Final_Taxonomy_Data.csv")

# --- Part B: TRY Data ---
Data_file <- "Data/TRY_data.csv"
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
# Extract genus names from the species names in TRY data
Genus_list <- sapply(strsplit(Data$AccSpeciesName, " "), `[`, 1)
Genus_list_unique <- unique(Genus_list)

# Retrieve taxonomic IDs from NCBI for TRY data
id_list <- get_uid(Genus_list_unique, ask = FALSE)
id_list_result <- attributes(id_list)
Genus_found <- Genus_list_unique[!is.na(id_list_result$uri)]
Genus_notfound <- Genus_list_unique[is.na(id_list_result$uri)]
class_1 <- classification(id_list, db = "ncbi")
rank_cols <- c("genus", "family", "order", "class", "phylum", "kingdom")
out_df1 <- extract_taxonomy(class_1, Genus_found)

# Provide manual corrections for missing TRY taxonomy
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
write.csv(Taxa_genus_TRY, "Data/Taxa_TRY_Genus.csv")

# Merge taxonomic information with TRY species data
Taxa_genus_TRY$Genus <- Taxa_genus_TRY$genus
Data$Genus <- Genus_list
Data2_TRY <- merge(Data, Taxa_genus_TRY, by = "Genus", all.x = TRUE)
write.csv(Data2_TRY, "Data/TRY_Final_Processed.csv", row.names = FALSE)

#############################################
# Step 4: Merge TRY Data and Flora of China Data
#############################################
library(stringr)

# Define file paths for the preprocessed FoC and TRY datasets
foc_file <- "Data/FoC_fin_20240723.csv"
try_file <- "Data/TRY_fin_20240729.csv"
merged_file <- "Data/Data_merged_FoC_TRY_20240730.csv"

# Read the datasets
Data1 <- read.csv(foc_file, header = TRUE, row.names = 1)
Data2 <- read.csv(try_file, header = TRUE, row.names = 1)

# Extract genus from species names in the FoC dataset and standardize species names
Data1$genus <- sapply(strsplit(Data1$Species, " "), `[`, 1)
capitalize_species <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}
Data1$Species <- sapply(Data1$Species, capitalize_species)
Data2$Species <- sapply(Data2$Species, capitalize_species)

# Add dataset identifiers
Data1$Datasets <- "FoC"
Data2$Datasets <- "TRY"

# Save preprocessed datasets (if needed)
write.csv(Data1, "Data/FoC_fin_20240730.csv", row.names = FALSE)
write.csv(Data2, "Data/TRY_fin_20240730.csv", row.names = FALSE)

# Combine the datasets
Data_merge <- rbind(Data1, Data2)
color_cols <- c("White", "Yellow", "Red", "Blue", "Purple", "Green", "NoDescription")
rank_cols <- c("genus", "family", "order", "class", "phylum", "kingdom", "GBIF_data", "GBIF_datacount")

# Identify species present in both datasets and aggregate their color data (binary)
sp_count <- table(Data_merge$Species)
multisp_list <- names(sp_count[sp_count > 1])
out_df_color <- data.frame(matrix(nrow = 0, ncol = length(color_cols)))
for (species in multisp_list) {
  temp_df <- Data_merge[Data_merge$Species == species, color_cols]
  temp_df <- apply(temp_df, 2, function(x) as.numeric(sum(x) > 0))
  out_df_color <- rbind(out_df_color, temp_df)
}
colnames(out_df_color) <- color_cols
out_df_color$Species <- multisp_list
out_df_color <- out_df_color[order(out_df_color$Species), ]

# Process species appearing in both datasets
Data_merge_multi <- Data_merge[Data_merge$Species %in% multisp_list, ]
Data_merge_multi <- Data_merge_multi[!duplicated(Data_merge_multi$Species), ]
Data_merge_multi <- Data_merge_multi[order(Data_merge_multi$Species), ]
Data_merge_multi$FoC <- 1
Data_merge_multi$TRY <- 1
Data_merge_multi <- cbind(out_df_color, Data_merge_multi[rank_cols])

# Process species appearing in only one dataset
Data_merge_single <- Data_merge[!Data_merge$Species %in% multisp_list, ]
Data_merge_single$FoC <- 0
Data_merge_single$TRY <- 0
Data_merge_single$FoC[Data_merge_single$Species %in% Data1$Species] <- 1
Data_merge_single$TRY[Data_merge_single$Species %in% Data2$Species] <- 1
Data_merge_single <- Data_merge_single[, colnames(Data_merge_single) != "Datasets"]

# Combine and save the merged dataset
Data_merged <- rbind(Data_merge_single, Data_merge_multi)
Data_merged <- Data_merged[order(Data_merged$Species), ]
write.csv(Data_merged, merged_file, row.names = FALSE)

# Resolve missing taxonomy data by borrowing known taxonomic ranks
Data_merged <- read.csv(merged_file, header = TRUE, row.names = 1)
Notaxa_family <- unique(Data_merged[Data_merged$kingdom == "NoData", "family"])
for (family in Notaxa_family) {
  temp_df <- Data_merged[Data_merged$family == family, ]
  if (nrow(temp_df) == 1) {
    print(family)
    next
  }
  temp_df_valid <- temp_df[temp_df$kingdom != "NoData", rank_cols]
  if (nrow(temp_df_valid) == 0) {
    print(family)
    next
  }
  temp_df_valid <- temp_df_valid[1, rank_cols]
  for (rank in rank_cols) {
    Data_merged[Data_merged$family == family, rank] <- temp_df_valid[rank]
  }
}
write.csv(Data_merged, merged_file, row.names = FALSE)
Data_merged <- read.csv(merged_file, header = TRUE, row.names = 1)

#############################################
# Step 5: Environmental Data Preprocessing
#############################################
# Load the merged species occurrence data
Data <- read.csv("InputData/MergedSpeciesData.csv", header = TRUE, row.names = 1)
Loc_dir <- "GBIF_Processed_Data/"
Loc_list <- list.files(Loc_dir, pattern = ".csv")

# Check for the existence of environmental data files for each species
file_exists <- c()
for (sp in Data$Species) {
  sp_file <- Lc(sp, ".csv")
  file_exists <- c(file_exists, file.exists(Lc(Loc_dir, sp_file)))
}

# Define the environmental (climate) and soil parameter columns to extract
env_cols <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10",
              "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "Alt_df")
soil_cols <- c("bdod_0.5cm_mean", "cec_0.5cm_mean", "cfvo_0.5cm_mean", "clay_0.5cm_mean", "nitrogen_0.5cm_mean",
               "ocd_0.5cm_mean", "ocs_0.30cm_mean", "phh2o_0.5cm_mean", "sand_0.5cm_mean", "silt_0.5cm_mean",
               "soc_0.5cm_mean")

# Filter species for which GBIF environmental data files exist
GBIF_sp <- Data$Species[file_exists]

# Process each species to compute summary environmental statistics
out_list <- list()
for (sp in GBIF_sp) {
  sp_file <- Lc(sp, ".csv")
  env_df <- read.csv(Lc(Loc_dir, sp_file), row.names = 1)
  
  # Remove rows with missing values and duplicates for the selected variables
  env_df <- na.omit(env_df[c(env_cols, soil_cols)])
  env_df <- env_df[!duplicated(env_df), ]
  
  # Skip species with fewer than 6 occurrence records
  if (nrow(env_df) < 6) next
  
  # Compute the trimmed mean (20% trim) for each variable
  env_list <- apply(env_df, 2, function(x) mean(x, trim = 0.2))
  
  # Compute soil heterogeneity using the average Euclidean distance among soil variables
  env_list["Soil_dist"] <- mean(dist(env_df[soil_cols], method = "euclidean"))
  
  # Record the number of accepted GBIF occurrence records and the species name
  env_list["GBIF_accept"] <- nrow(env_df)
  env_list["Species"] <- sp
  
  out_list[[sp]] <- env_list
}
out_df <- data.frame(t(data.frame(out_list)))

# Merge the computed environmental summaries with the original species data
out_df <- merge(Data, out_df, by = "Species", all = TRUE)
write.csv(out_df, "OutputData/MergedSpeciesData_withEnvironment.csv")

#############################################
# Step 6: Create Environmental Axes using PCA
#############################################
source("Function_source.R")  # Load any additional required functions

# Define input and output file paths for PCA results
data_file <- "Data/Data_merged_Angiosperms_withEnvironment.csv"
output_csv <- "Data/Data_for_Analysis.csv"
output_pc_df <- "Output/PC_df"
output_pca_rda <- "Output/Variable_PCA.rda"
output_pca_coord <- "Output/PCA_coordinates.csv"

# Read the dataset for analysis and remove rows missing GBIF occurrence information
Data <- read.csv(data_file, header = TRUE, row.names = 1)
Data <- Data[!is.na(Data$GBIF_accept), ]

# Define the list of environmental and soil variables to include in the PCA
env_cols <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", 
              "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "Alt_df")
soil_cols <- c("bdod_0.5cm_mean", "cec_0.5cm_mean", "cfvo_0.5cm_mean", "clay_0.5cm_mean", 
               "nitrogen_0.5cm_mean", "ocd_0.5cm_mean", "ocs_0.30cm_mean", "phh2o_0.5cm_mean", 
               "sand_0.5cm_mean", "silt_0.5cm_mean", "soc_0.5cm_mean", "Soil_dist")
variable_list <- c(env_cols, soil_cols)

# Standardize the variables before performing PCA
Data[variable_list] <- apply(Data[variable_list], 2, scale)

# Perform PCA and retain 10 principal components
Variable_PCA <- PCA(Data[, variable_list], scale.unit = TRUE, ncp = 10, graph = FALSE)
print(Variable_PCA$eig)      # Print eigenvalues (variance explained)
print(Variable_PCA$var$coord)  # Print variable coordinates in PCA space

# Create a data frame of PCA scores for each species
PC_names <- paste0("PC", 1:10)
PC_df <- data.frame(Variable_PCA$ind$coord[, 1:10])
colnames(PC_df) <- PC_names
PC_df$Species <- Data$Species

# Merge the PCA scores with the original dataset and save the results
Data <- merge(Data, PC_df, by = "Species", all = TRUE)
write.csv(Data, output_csv, row.names = FALSE)
save(PC_df, file = output_pc_df)
save(Variable_PCA, file = output_pca_rda)
write.csv(Variable_PCA$var$coord, output_pca_coord, row.names = TRUE)

cat("Processing complete.\n")
