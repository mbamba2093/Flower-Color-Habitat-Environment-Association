# Load required functions from the external source file
source("Function_source.R")

#############################################
# MCMCglmm Main Results Analysis
# Comments have been added by ChatGPT o3-mini-high
#############################################

# Set the target color for the analysis (choose one: "White", "Yellow", "RedType")
Target_color <- "RedType"   # Here, we choose RedType
Prefix <- "ALL"
Out_dir <- "Output/mcmcGLMM/"
Data_file <- "Data/Data_for_Analysis.csv" 
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Tree_file <- "Data/TimeTree_all_filt.nwk"
Tree <- read.tree(Tree_file)

# Replace spaces with underscores in species names to match tree tip labels
Data$Species <- str_replace(Data$Species, " ", "_")

# Remove rows that lack environmental information (i.e., missing GBIF_accept)
Data <- Data[!is.na(Data$GBIF_accept), ]

# Create an aggregated red type variable (if species is red, blue, or purple, then 1)
Data$RedType <- as.numeric(Data$Red == 1 | Data$Blue == 1 | Data$Purple == 1)

# Construct the analysis data frame including species names, target variable, and PCA variables (assumed to be named "PC...")
analysis_df <- data.frame(
  Species = Data$Species,
  Values = Data[, Target_color],
  Data[colnames(Data)[grep("PC", colnames(Data))]]
)

# Match the species between the tree and the phenotype data:
matched_species <- intersect(Tree$tip.label, analysis_df$Species)
# Remove tips from the tree that are not present in the phenotype data
Tree <- drop.tip(Tree, setdiff(Tree$tip.label, matched_species))
analysis_df_filt <- analysis_df[analysis_df$Species %in% matched_species, ]

# Create a phylogenetic covariance matrix (Ainv) from the tree.
epsilon <- min(Tree$edge.length[Tree$edge.length != 0])
Tree$edge.length[Tree$edge.length == 0] <- epsilon
# Convert tree to an ultrametric tree using Grafen’s method
ultrametric_tree <- compute.brlen(Tree, method = "Grafen")
Ainv <- inverseA(ultrametric_tree)$Ainv

# Set up the MCMC parameters
prior <- list(
  G = list(G1 = list(V = 1, nu = 0.002)), 
  R = list(V = 1, fix = 1)
)
nitt <- 1050000   # Total number of iterations
burnin <- 50000   # Burn-in period
thin <- 100       # Thinning interval

# Convert the response variable to character (for threshold model)
analysis_df_filt$Values <- as.character(analysis_df_filt$Values)
# Define the base formula using the first 10 principal components
base_formula <- as.formula("Values ~ PC.1 + PC.2 + PC.3 + PC.4 + PC.5 + PC.6 + PC.7 + PC.8 + PC.9 + PC.10")

# Fit the phylogenetic generalized linear mixed model (PGLMM) with threshold family
model_pglmm <- MCMCglmm(
  base_formula, 
  random = ~ Species, 
  family = "threshold", 
  ginverse = list(Species = Ainv), 
  data = analysis_df_filt, 
  prior = prior, 
  nitt = nitt, 
  burnin = burnin, 
  thin = thin, 
  verbose = TRUE
)

# Save the main results and convergence diagnostics
summary_pglm <- summary(model_pglmm)
write.csv(summary_pglm$solutions, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, ".csv"))
write.csv(data.frame(model_pglmm$Sol), Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_mcmc.csv"))
effsize <- data.frame(effectiveSize(model_pglmm$Sol))
colnames(effsize) <- Target_color
write.csv(effsize, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_effsize.csv"))
save(model_pglmm, file = Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, ".rda"))







#############################################
# MCMCglmm Analysis for Each Environmental Variable
#############################################

# Set analysis parameters
Target_color <- "RedType"  # Choose from: White, Yellow, RedType, Red, Blue, Purple
Prefix <- "ALL"            # Could also be a family name (e.g., Asteraceae, Fabaceae, etc.)
Out_dir <- "Output/mcmcGLMM_eachVariable/"
Data_file <- "Data/Data_for_Analysis.csv"  # PCA variables are already included
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Tree_file <- "Data/TimeTree_all_filt.nwk"
Tree <- read.tree(Tree_file)

# Replace spaces with underscores to match tree tip labels
Data$Species <- str_replace(Data$Species, " ", "_")
# Remove rows missing environmental information
Data <- Data[!is.na(Data$GBIF_accept), ]
# Uncomment and adjust the following line if you want to analyze a single family
# Data <- Data[Data$family == "Orchidaceae", ]

# Create the aggregated red type variable
Data$RedType <- as.numeric(Data$Red == 1 | Data$Blue == 1 | Data$Purple == 1)

# Build the analysis data frame containing species names, the target variable, and each environmental variable (columns 20 to 51)
analysis_df <- data.frame(
  Species = Data$Species,
  Values = Data[, Target_color],
  Data[colnames(Data)[20:51]]  # Adjust column indices if necessary
)

# Match species between the tree and analysis data
matched_species <- intersect(Tree$tip.label, analysis_df$Species)
Tree <- drop.tip(Tree, setdiff(Tree$tip.label, matched_species))
analysis_df_filt <- analysis_df[analysis_df$Species %in% matched_species, ]

# Create a phylogenetic covariance matrix from the tree
epsilon <- min(Tree$edge.length[Tree$edge.length != 0])
Tree$edge.length[Tree$edge.length == 0] <- epsilon
ultrametric_tree <- compute.brlen(Tree, method = "Grafen")
Ainv <- inverseA(ultrametric_tree)$Ainv

# Set up MCMC parameters
prior <- list(
  G = list(G1 = list(V = 1, nu = 0.002)), 
  R = list(V = 1, fix = 1)
)
nitt <- 1050000   # Total iterations
burnin <- 50000   # Burn-in period
thin <- 100       # Thinning interval

# Convert response variable to character for threshold modeling
analysis_df_filt$Values <- as.character(analysis_df_filt$Values)

# Get the names of the environmental variables (columns 3 onward)
variable_names <- colnames(analysis_df_filt)[3:length(colnames(analysis_df_filt))]

# Loop over each environmental variable
for (i in variable_names) {
  # Standardize the predictor variable
  analysis_df_filt[i] <- scale(analysis_df_filt[i])
  
  # Create a model formula for the current variable
  base_formula <- as.formula(Lc("Values ~ ", i))
  
  # Fit the MCMCglmm model for the current variable
  model_pglmm <- MCMCglmm(
    base_formula, 
    random = ~ Species, 
    family = "threshold", 
    ginverse = list(Species = Ainv), 
    data = analysis_df_filt, 
    prior = prior, 
    nitt = nitt, 
    burnin = burnin, 
    thin = thin, 
    verbose = TRUE
  )
  
  # Save the output for the current variable
  summary_pglm <- summary(model_pglmm)
  write.csv(summary_pglm$solutions, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_", i, ".csv"))
  write.csv(data.frame(model_pglmm$Sol), Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_", i, "_mcmc.csv"))
  effsize <- data.frame(effectiveSize(model_pglmm$Sol))
  colnames(effsize) <- Target_color
  write.csv(effsize, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_", i, "_effsize.csv"))
  save(model_pglmm, file = Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_", i, ".rda"))
}







#############################################
# MCMCglmm Analysis for Each Family with Sufficient Sample Size
#############################################

# Define output directory and data files
Out_dir <- "Output/mcmcGLMM_Family8/"
Data_file <- "Data/Data_for_Analysis.csv"  # PCA variables are included
Tree_file <- "Data/TimeTree_all_filt.nwk"
# Set the list of target colors for analysis (e.g., "White", "Yellow", "RedType")
Target_color_list <- c("White", "Yellow", "RedType")
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
# Filter families with more than 50 samples
family_list <- names(table(Data$family)[table(Data$family) > 50])

# Loop over each family
for (i1 in family_list) {
  Prefix <- i1
  Data <- read.csv(Data_file, header = TRUE, row.names = 1)
  # Create aggregated red type variable
  Data$RedType <- as.numeric(Data$Red == 1 | Data$Blue == 1 | Data$Purple == 1)
  Tree <- read.tree(Tree_file)
  Data$Species <- str_replace(Data$Species, " ", "_")
  Data <- Data[!is.na(Data$GBIF_accept), ]
  # Select only the data for the current family
  Data <- Data[Data$family == i1, ]
  
  # Build analysis data frame with species names and PCA variables
  analysis_df <- data.frame(
    Species = Data$Species,
    Data[colnames(Data)[grep("PC", colnames(Data))]]
  )
  
  # Match species with the tree and filter
  matched_species <- intersect(Tree$tip.label, analysis_df$Species)
  Tree <- drop.tip(Tree, setdiff(Tree$tip.label, matched_species))
  analysis_df_filt <- analysis_df[analysis_df$Species %in% matched_species, ]
  
  # Skip the family if sample size is too small
  if (nrow(analysis_df_filt) < 50) {
    next
  }
  
  # Loop over each target color
  for (i2 in Target_color_list) {
    Target_color <- i2
    # Add the phenotype values for the target color into the analysis data
    analysis_df_filt$Values <- Data[Data$Species %in% matched_species, Target_color]
    
    # Skip if there is no variation in the phenotype data
    if (length(unique(analysis_df_filt$Values)) == 1) {
      next
    }
    
    # Create phylogenetic covariance matrix (Ainv)
    epsilon <- min(Tree$edge.length[Tree$edge.length != 0])
    Tree$edge.length[Tree$edge.length == 0] <- epsilon
    ultrametric_tree <- compute.brlen(Tree, method = "Grafen")
    Ainv <- inverseA(ultrametric_tree)$Ainv
    
    # Set up MCMC parameters (using a different prior for family-level analysis)
    prior <- list(
      G = list(G1 = list(V = 0.01, nu = 0.1)), 
      R = list(V = 1, fix = 1)
    )
    nitt <- 1050000
    burnin <- 50000
    thin <- 100
    
    # Convert response to character for threshold modeling
    analysis_df_filt$Values <- as.character(analysis_df_filt$Values)
    base_formula <- as.formula("Values ~ PC.1 + PC.2 + PC.3 + PC.4 + PC.5 + PC.6 + PC.7 + PC.8 + PC.9 + PC.10")
    
    # Run three MCMC chains (or replicate runs) for stability
    for (i3 in 1:3) {
      model_pglmm <- MCMCglmm(
        base_formula, 
        random = ~ Species, 
        family = "threshold", 
        ginverse = list(Species = Ainv), 
        data = analysis_df_filt, 
        prior = prior, 
        nitt = nitt, 
        burnin = burnin, 
        thin = thin, 
        verbose = TRUE
      )
      
      # Save outputs for the current run
      summary_pglm <- summary(model_pglmm)
      write.csv(summary_pglm$solutions, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, i3, ".csv"))
      write.csv(data.frame(model_pglmm$Sol), Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, i3, "_mcmc.csv"))
      effsize <- data.frame(effectiveSize(model_pglmm$Sol))
      colnames(effsize) <- Target_color
      write.csv(effsize, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, i3, "_effsize.csv"))
      save(model_pglmm, file = Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, i3, ".rda"))
    }
  }
}






#############################################
# MCMCglmm Analysis Without Phylogenetic Information
# (Using Taxonomic Information: Family and Genus as random effects)
#############################################

# Set analysis parameters
Target_color <- "RedType"   # Choose from: White, Yellow, RedType, Red, Blue, Purple
Prefix <- "ALL"             # Could also be set to a specific family name
Out_dir <- "Output/GLMM/"
Data_file <- "Data/Data_for_Analysis.csv"  # PCA variables are included
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Tree_file <- "Data/TimeTree_all_filt.nwk"
Tree <- read.tree(Tree_file)

# Prepare the data: match species names to tree labels and remove rows missing environmental info
Data$Species <- str_replace(Data$Species, " ", "_")
Data <- Data[!is.na(Data$GBIF_accept), ]
# Create aggregated red type variable
Data$RedType <- as.numeric(Data$Red == 1 | Data$Blue == 1 | Data$Purple == 1)

# Build an analysis data frame including species names, genus, family, the target variable, and PCA variables
analysis_df <- data.frame(
  Species = Data$Species,
  Genus = Data$genus,
  Family = Data$family,
  Values = Data[, Target_color],
  Data[colnames(Data)[grep("PC", colnames(Data))]]
)

# Set up MCMC parameters with two random effects: Family and Genus
prior <- list(
  G = list(
    G1 = list(V = 1, nu = 0.002), 
    G2 = list(V = 1, nu = 0.002)
  ), 
  R = list(V = 1, fix = 1)
)
nitt <- 1050000   # Total iterations
burnin <- 50000   # Burn-in period
thin <- 100       # Thinning interval

# Define the base model formula (using the first 10 principal components)
base_formula <- as.formula("Values ~ PC.1 + PC.2 + PC.3 + PC.4 + PC.5 + PC.6 + PC.7 + PC.8 + PC.9 + PC.10")

# Fit the model without including phylogenetic covariance (using Family and Genus as random effects)
model_pglmm <- MCMCglmm(
  base_formula, 
  random = ~ Family + Genus, 
  family = "threshold", 
  data = analysis_df, 
  prior = prior, 
  nitt = nitt, 
  burnin = burnin, 
  thin = thin, 
  verbose = TRUE
)

# Save the output and convergence diagnostics
summary_pglm <- summary(model_pglmm)
write.csv(summary_pglm$solutions, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, ".csv"))
write.csv(data.frame(model_pglmm$Sol), Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_mcmc.csv"))
effsize <- data.frame(effectiveSize(model_pglmm$Sol))
colnames(effsize) <- Target_color
write.csv(effsize, Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, "_effsize.csv"))
save(model_pglmm, file = Lc(Out_dir, "mcmcGLMM_", Prefix, Target_color, ".rda"))


#############################################
# GLM Analysis Without Random Effects
# (A standard logistic regression model)
#############################################

# Set the target color and prefix for the analysis
Target_color <- "Yellow"   # Choose the target phenotype (e.g., White, Yellow, RedType, etc.)
Prefix <- ""               # Could be set to a specific group (e.g., Asteraceae, etc.)
Out_dir <- "Output/GLM/"
Data_file <- "Data/Data_for_Analysis.csv"  # Data with PCA variables included
Data <- read.csv(Data_file, header = TRUE, row.names = 1)
Tree_file <- "Data/TimeTree_all_filt.nwk"
Tree <- read.tree(Tree_file)

# Format species names to match the tree labels
Data$Species <- str_replace(Data$Species, " ", "_")
# Remove rows missing environmental information
Data <- Data[!is.na(Data$GBIF_accept), ]
# Create the aggregated red type variable (if needed)
Data$RedType <- as.numeric(Data$Red == 1 | Data$Blue == 1 | Data$Purple == 1)

# Build the analysis data frame with species names, the target variable, and PCA variables
analysis_df <- data.frame(
  Species = Data$Species,
  Values = Data[, Target_color],
  Data[colnames(Data)[grep("PC", colnames(Data))]]
)

# Match the species between the tree and phenotype data
matched_species <- intersect(Tree$tip.label, analysis_df$Species)
Tree <- drop.tip(Tree, setdiff(Tree$tip.label, matched_species))
analysis_df_filt <- analysis_df[analysis_df$Species %in% matched_species, ]

# Convert the response variable to a factor (for logistic regression)
analysis_df_filt$Values <- as.factor(analysis_df_filt$Values)
# Define the GLM formula using the first 10 principal components as predictors
base_formula <- as.formula("Values ~ PC.1 + PC.2 + PC.3 + PC.4 + PC.5 + PC.6 + PC.7 + PC.8 + PC.9 + PC.10")

# Fit the logistic regression model (binomial family with logit link)
glm_model <- glm(formula = base_formula, data = analysis_df_filt, family = binomial("logit"))
result_glm <- summary(glm_model)
result_glm <- data.frame(result_glm$coefficients)

# Perform an ANOVA on the GLM model
result_anova <- Anova(glm_model, test = "F")
result_anova <- data.frame(result_anova)

# Save the coefficient estimates and ANOVA results to CSV files
write.csv(result_glm, Lc(Out_dir, Target_color, "_glm_coef.csv"))
write.csv(result_anova, Lc(Out_dir, Target_color, "_glm_anova.csv"))
