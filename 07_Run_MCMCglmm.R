#####
# Source custom functions
#####
source("Function_source.R")

#####
# User settings
# Change these manually for each run.
# Available colors: White, Yellow, Red, Blue, Purple, Green, MultiColor
#####
data_file <- "Data3/Data_for_withoutPhylogeny.csv"
phylo_file <- "Data3/Phylogeny_all_S1.rda"
Out_dir <- "Output_20260120_MCMCglmm_main/"

Prefix <- "MS_1_main"
Target_color <- "MultiColor"
Replicate_value <- "3"

#####
# Input data
#####
Data <- fread(data_file, header = TRUE, data.table = FALSE)

# Derived color classes used in the analysis
Data$RedBluePurple <- as.numeric(rowSums(Data[, c("Red", "Blue", "Purple")]) > 0)
Data$MultiColor <- as.numeric(rowSums(Data[, c("White", "Yellow", "RedBluePurple")]) > 1)

#####
# Phylogenetic data preprocessing
#####
load(phylo_file) # res_lcvp
res_lcvp$species.list$species <- str_replace(res_lcvp$species.list$species, " ", "_")
sp_list1 <- res_lcvp$species.list$species[res_lcvp$species.list$status == "prune"]

Tree <- drop.tip(res_lcvp$scenario.1, setdiff(res_lcvp$scenario.1$tip.label, sp_list1))
Tree$node.label <- NULL

Data$Species <- str_replace(Data$wfo_name, " ", "_")
Data <- Data[Data$Species %in% sp_list1, ]

#####
# Analysis dataset
#####
analysis_df <- data.frame(
  Species = Data$Species,
  Values = Data[, Target_color],
  Data[c("PC1", "PC2", "PC3", "PC4", "PC5")],
  Data[c("Con_OCEANIA", "Con_AFRICA", "Con_NORTHAMERICA", "Con_EUROPE", "Con_SOUTHAMERICA", "Con_ASIA")]
)

# Scale PC axes
pc_columns <- c("PC1", "PC2", "PC3", "PC4", "PC5")
analysis_df[, pc_columns] <- apply(analysis_df[, pc_columns], 2, scale)

#####
# Phylogenetic covariance matrix (Ainv)
#####
scale_tree_unit_diag <- function(tr) {
  V <- ape::vcv(tr)
  s <- mean(diag(V))
  tr$edge.length <- tr$edge.length / s
  tr
}

ultrametric_tree <- compute.brlen(Tree, method = "Grafen")
ultrametric_tree <- scale_tree_unit_diag(ultrametric_tree)
Ainv <- inverseA(ultrametric_tree)$Ainv

#####
# MCMCglmm settings
#####
prior <- list(
  G = list(G1 = list(V = 1, nu = 1)),
  R = list(V = 1, fix = 1)
)

nitt <- 1050000
burnin <- 50000
thin <- 100

analysis_df$Values <- as.factor(analysis_df$Values)
analysis_df[, c("Con_OCEANIA", "Con_AFRICA", "Con_NORTHAMERICA", "Con_EUROPE", "Con_SOUTHAMERICA", "Con_ASIA")] <-
  apply(
    analysis_df[, c("Con_OCEANIA", "Con_AFRICA", "Con_NORTHAMERICA", "Con_EUROPE", "Con_SOUTHAMERICA", "Con_ASIA")],
    2,
    as.factor
  )

base_formula <- as.formula(
  "Values ~ PC1 + PC2 + PC3 + PC4 + PC5 + Con_OCEANIA + Con_AFRICA + Con_NORTHAMERICA + Con_EUROPE + Con_SOUTHAMERICA + Con_ASIA"
)

#####
# Model fitting
#####
model_pglmm <- MCMCglmm(
  base_formula,
  random = ~ Species,
  family = "threshold",
  ginverse = list(Species = Ainv),
  data = analysis_df,
  prior = prior,
  nitt = nitt,
  burnin = burnin,
  thin = thin,
  verbose = TRUE
)

#####
# Output files
#####
summary_pglm <- summary(model_pglmm)
write.csv(summary_pglm$solutions, Lc(Out_dir, Prefix, Target_color, "_", Replicate_value, ".csv"))
write.csv(data.frame(model_pglmm$Sol), Lc(Out_dir, Prefix, Target_color, "_", Replicate_value, "_mcmc.csv"))

effsize <- data.frame(effectiveSize(model_pglmm$Sol))
colnames(effsize) <- Target_color
write.csv(effsize, Lc(Out_dir, Prefix, Target_color, "_", Replicate_value, "_effsize.csv"))

save(model_pglmm, file = Lc(Out_dir, Prefix, Target_color, "_", Replicate_value, ".rda"))
