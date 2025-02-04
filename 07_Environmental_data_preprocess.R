# Load the merged species occurrence data
Data <- read.csv("InputData/MergedSpeciesData.csv", header = TRUE, row.names = 1)

# Define the directory containing species-specific GBIF environmental data
Loc_dir <- "GBIF_Processed_Data/"
Loc_list <- list.files(Loc_dir, pattern=".csv")

# Check if environmental data files exist for each species in the dataset
file_exists <- c()
for (i in Data$Species){
  sp_file <- Lc(i, ".csv")
  file_exists <- c(file_exists, file.exists(Lc(Loc_dir, sp_file)))
}

# Define environmental and soil parameter columns to extract
env_cols <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10",
              "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "Alt_df")

soil_cols <- c("bdod_0.5cm_mean", "cec_0.5cm_mean", "cfvo_0.5cm_mean", "clay_0.5cm_mean", "nitrogen_0.5cm_mean",
               "ocd_0.5cm_mean", "ocs_0.30cm_mean", "phh2o_0.5cm_mean", "sand_0.5cm_mean", "silt_0.5cm_mean",
               "soc_0.5cm_mean")

# Filter species for which GBIF environmental data files exist
GBIF_sp <- Data$Species[file_exists]

# Initialize an empty list to store processed environmental summaries
out_list <- list()
temp_num <- 1

# Process each species with available environmental data
for (i in GBIF_sp){
  sp_file <- Lc(i, ".csv")  # Construct file path
  env_df <- read.csv(Lc(Loc_dir, sp_file), row.names = 1)  # Read species-specific environmental data
  
  # Remove rows with missing values in environmental and soil data
  env_df <- na.omit(env_df[c(env_cols, soil_cols)])
  
  # Remove duplicate environmental records
  env_df <- env_df[!duplicated(env_df), ]
  
  # Skip species with fewer than 6 occurrence records
  if (nrow(env_df) < 6){
    next
  }
  
  # Compute trimmed mean (excluding extreme 20%) for each environmental variable
  env_list <- apply(env_df, 2, function(x){mean(x, trim = 0.2)})
  
  # Compute soil environmental heterogeneity using Euclidean distance
  env_list["Soil_dist"] <- mean(dist(env_df[soil_cols], method = "euclidean"))
  
  # Store the number of GBIF occurrence records used
  env_list["GBIF_accept"] <- nrow(env_df)
  
  # Store species name
  env_list["Species"] <- i
  
  # Append to the output list
  out_list[[i]] <- env_list
}

# Convert the list into a data frame
out_df <- data.frame(t(data.frame(out_list)))

# Merge the new environmental data with the original species dataset
out_df <- merge(Data, out_df, by = "Species", all = TRUE)

# Save the final dataset with environmental attributes
write.csv(out_df, "OutputData/MergedSpeciesData_withEnvironment.csv")
