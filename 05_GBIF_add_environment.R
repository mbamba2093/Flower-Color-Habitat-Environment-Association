source("GBIF_function.R")

# Instructions to download ISRIC soil data (for reference)
# 1. Get environmental data from ISRIC
# wget -r -np -nH https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_0-5cm_mean/
# wget -r -np -nH https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_0-5cm_mean.vrt
# 2. Move them to a suitable directory (the .vrt file and the main directory should be in the same location).

# Define directory containing GBIF occurrence data
Loc_dir = "GBIF_Processed_Data/"

##### Load environmental raster data

### Load soil parameters (GeoTiff files)
GeoTiff_dir1 = "Environmental_Data/Soil/"
GeoTiff1_list <- list.files(GeoTiff_dir1, pattern=".vrt")

GeoTiff1 <- list()
for (i2 in 1:length(GeoTiff1_list)){
  GeoTiff1[[i2]] <- rast(Lc(GeoTiff_dir1, GeoTiff1_list[i2]))
}

### Load WRB soil classification raster data
GeoTiff_dir2 = "Environmental_Data/Soil/WRB/"
GeoTiff2_list <- list.files(GeoTiff_dir2, pattern=".vrt")

GeoTiff2 <- list()
for (i3 in 1:length(GeoTiff2_list)){
  GeoTiff2[[i3]] <- rast(Lc(GeoTiff_dir2, GeoTiff2_list[i3]))
}

print("TIFF files have been loaded.")

# List all species occurrence data files
Loc_list <- list.files(Loc_dir, pattern=".csv")
endnum = length(Loc_list)
pb <- txtProgressBar(min = 0, max = endnum, style = 3)
startnum = 1

# Process each species occurrence dataset
for (i1 in Loc_list){
  gbif_data <- read.csv(Lc(Loc_dir, i1), header = TRUE, row.names = 1)
  
  # Extract soil environmental data for each occurrence point
  for (i2 in 1:length(GeoTiff1)){
    gbif_data <- Get_soilenv(gbif_data, GeoTiff1[[i2]])
  }
  
  # Extract WRB classification data
  for (i3 in 1:length(GeoTiff2)){
    gbif_data <- Get_soilenv(gbif_data, GeoTiff2[[i3]])
  }
  
  # Save updated dataset with soil parameters
  write.csv(gbif_data, Lc(Loc_dir, i1))
  
  # Update progress bar
  setTxtProgressBar(pb, startnum)
  startnum <- startnum + 1
}

##### Climate Data Extraction

# Define input and output directories
Loc_dir = "GBIF_Processed_Data/"
Out_dir = "GBIF_Processed_Data_WorldClim/"

# Download climate data from WorldClim (10 km resolution)
Clim_alt <- raster::getData('worldclim', var='alt', res=5)  # Elevation data
Clim_bio <- raster::getData('worldclim', var='bio', res=5)  # Bioclimatic variables

# List species occurrence data files
Loc_list <- list.files(Loc_dir, pattern=".csv")
endnum = length(Loc_list)
pb <- txtProgressBar(min = 0, max = endnum, style = 3)
startnum = 1

# Process each species occurrence dataset
for (i1 in Loc_list){
  gbif_data <- read.csv(Lc(Loc_dir, i1), header = TRUE, row.names = 1)
  
  # Extract coordinates for climate data extraction
  Sp_loc <- data.frame(lon = gbif_data$decimalLongitude,
                       lat = gbif_data$decimalLatitude)
  
  # Extract climate and elevation data for occurrence points
  Clim_df <- Extract_clim(Sp_loc, Clim_bio)
  Alt_df <- Extract_clim(Sp_loc, Clim_alt)
  
  # Merge climate and elevation data with species occurrence dataset
  gbif_data <- cbind(gbif_data, Clim_df, Alt_df)
  
  # Save updated dataset with climate and elevation data
  write.csv(gbif_data, Lc(Out_dir, i1))
  
  # Update progress bar
  setTxtProgressBar(pb, startnum)
  startnum <- startnum + 1
}


