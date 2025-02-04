#####
source("Functions/CustomFunctions.R")  # Load external functions
#####

# Load the merged species data with environmental variables
Data_file <- "ProcessedData/MergedSpeciesData_withEnvironment.csv"
Data <- read.csv(Data_file, header = TRUE, row.names = 1)

# Remove species without valid GBIF occurrence data
Data <- Data[!is.na(Data$GBIF_accept), ]

# Define environmental and soil variable names
env_cols <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10",
              "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "Alt_df")

soil_cols <- c("bdod_0.5cm_mean", "cec_0.5cm_mean", "cfvo_0.5cm_mean", "clay_0.5cm_mean", "nitrogen_0.5cm_mean",
               "ocd_0.5cm_mean", "ocs_0.30cm_mean", "phh2o_0.5cm_mean", "sand_0.5cm_mean", "silt_0.5cm_mean",
               "soc_0.5cm_mean", "Soil_dist")

# Combine environmental and soil variables into a single list
Variable_list <- c(env_cols, soil_cols)

# Standardize environmental and soil variables (Z-score normalization)
Data[Variable_list] <- apply(Data[Variable_list], 2, scale)

# Perform Principal Component Analysis (PCA)
Variable_PCA <- PCA(Data[, Variable_list], scale.unit = TRUE, ncp = 10, graph = FALSE)

# Print eigenvalues and PCA variable coordinates
Variable_PCA$eig
Variable_PCA$var$coord

# PCA components 1 to 10 explain 95.776% of variance

# Create a dataframe containing PCA scores for species
PC_names <- interaction(rep("PC", 10), 1:10)  # Generate PC column names
PC_df <- data.frame(Variable_PCA$ind$coord[, 1:10])  # Extract first 10 principal components
colnames(PC_df) <- PC_names  # Assign column names

# Add species names to PCA dataframe
PC_df$Species <- Data$Species

# Merge PCA results with original dataset
Data <- merge(Data, PC_df, by = "Species", all = TRUE)

# Save processed dataset with PCA results
write.csv(Data, "ProcessedData/Data_for_Analysis.csv")

# Save PCA results for further analysis
save(PC_df, file = "Output/PCA_SpeciesScores.rda")
save(Variable_PCA, file = "Output/Environmental_PCA_Model.rda")

# Save PCA variable coordinates
write.csv(Variable_PCA$var$coord, "Output/PCA_VariableCoordinates.csv")
