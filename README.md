# Flower-Color-Habitat-Environment-Association
This repository contains the Python and R scripts used for the analysis in the tentative project titled "Association between flower color and habitat environment" by Bamba et al.

## Repository Contents
- `01_Scraping_FloraOfChina.py`: This script was used to extract descriptive information from the Flora of China. Web scraping of the Flora of China is not explicitly prohibited, and the scraping has been designed to minimize server load.
- `02_Extract_flower_color_description.py`: This script was used to extract flower color description from whole descriptive information using OpenAI api. 
- `03_Categorize_Flower_color.py`: This script was used to categorize flower color description to categorical data; white, yellow, red, blue, green and no_desription.
- `04_Data_PreProcessing.R`: This script was used to data pre-processing for statistical analysis; GBIF search, Add environmental data, and so on.
- `05_MCMCglmm.R`: This script was used to demonstrate MCMCglmm and other statistical anslyses; MCMCglmm with full model, each variable model, each family model, and so on.

