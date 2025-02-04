source("R_function.R")

# Define directories and filenames
Data_dir = "InputData/"  # Directory containing the input data file
Data_file = Lc(Data_dir, "SpeciesData.csv")  # Construct the full file path
Data = read.csv(Data_file, header = TRUE)  # Read the CSV file into a dataframe
Out_dir = "GBIF_Output/"  # Directory where GBIF results will be stored

# Create a copy of the data and initialize a new column for GBIF search results
Data2 <- Data
Data2$GBIF_search <- "NotYet"

# Loop through each row of the dataset to process species data
for (i1 in 1:nrow(Data)){
  Sp_name = Data[i1, "Species"]  # Extract species name from the dataset
  print(Lc(i1, " : ", Sp_name))  # Print the progress
  
  # Check if GBIF data for this species already exists
  if (file.exists(Lc(Out_dir, Sp_name, ".csv"))) {
    Data2[i1, "GBIF_search"] <- "Success"  # Mark as successfully retrieved
    print("skip")
    next  # Skip further processing for this species
  }
  
  # Retrieve occurrence data from GBIF
  Sp_loc = Search_GBIF(Sp_name)
  
  # If GBIF data is retrieved successfully
  if (mode(Sp_loc) == "list"){
    Occurance_status <- "Success"
    Data2[i1, "GBIF_search"] <- Occurance_status
    write.csv(Sp_loc, Lc(Out_dir, Sp_name, ".csv"))  # Save the data
  } else {
    Occurance_status <- Sp_loc
    Data2[i1, "GBIF_search"] <- Occurance_status  # Mark the failure reason
    next
  }
}

# Save the updated dataset with GBIF search results
write.csv(Data2, Lc(Out_dir, "Processed_SpeciesData.csv"))
