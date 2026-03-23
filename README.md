# Flower-Color-Habitat-Environment-Association
This repository contains the Python and R scripts used for the analysis in the tentative project titled "Association between flower color and habitat environment" by Bamba et al.

## Repository Contents
- `01_Scraping_FloraOfChina.py`: This script was used to extract descriptive information from the Flora of China. Web scraping of the Flora of China is not explicitly prohibited, and the scraping has been designed to minimize server load.
- `02_Extract_flower_color_description.py`: This script was used to extract flower color description from whole descriptive information using OpenAI api. 
- `03_Categorize_Flower_color.py`: This script was used to categorize flower color description to categorical data; white, yellow, red, blue, green and no_desription.
- `04_Taxonomic_correction.R`: This script was used for taxonomic correction and standardization of species names by matching the original species list to the World Flora Online database and preparing corrected datasets for subsequent analyses.

- `05_Data_integration_for_analysis.R`: This script was used for integrating corrected species data, trait data, and phylogenetic information to generate datasets for statistical analyses with and without phylogenetic information.

- `06_Check_consistency_and_build_networks.R`: This script was used to evaluate consistency among datasets, identify mismatched species records, and generate network data and figures summarizing relationships among the compared datasets.

- `07_Run_MCMCglmm_by_color_and_replicate.R`: This script was used to perform MCMCglmm analyses for each target color and replicate setting, generating model outputs for manually specified combinations of color category and replicate value.
