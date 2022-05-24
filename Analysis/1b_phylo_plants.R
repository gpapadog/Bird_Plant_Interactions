# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# The directory where the original data are:
data_path <- 'Data/'
# The directory where the processed data should be saved:
save_path <- 'Data/'
# Whether the processed data should be saved or not:
save_files <- TRUE


# --------- BEGINNING -------- #

# Setting the working directory.
setwd(wd_path)

# Loading libraries.
library(ape)
library(V.PhyloMaker)

dta <- fread(paste0(data_path, 'ATLANTIC_frugivory.csv'))

# Correcting recorded coordinates. Some entries have recorded coordinates that
# are not correct as they are not in Brazil. These coordinates are set to NA.
wh <- which(dta$Longitude > 0)
dta$Longitude[wh] <- dta$Longitude[wh] * (- 1)
wh <- which(dta$Longitude > - 30)
dta$Longitude[wh] <- NA
dta$Latitude[wh] <- NA

# Correcting a wrong recorded interaction.
dta <- subset(dta, Plant_Species != 'Amaioua hybridus')
dta$Fruit_color[dta$Plant_Species == 'Amaranthus hybridus'] <- 'red'

# Subsetting the data to only bird-plant interactions
birds <- subset(dta, Frug_Class == 'Aves')
# All species and interactions for studies that followed at least one bird:
dta_subset <- subset(dta, `Study reference` %in% birds$`Study reference`)

# Keeping track of the unique plant species in the study:
uni_plants <- cbind(species = dta_subset$Plant_Species,
                    genus = dta_subset$Plant_genus,
                    family = dta_subset$Plant_family)
uni_plants <- as.data.frame(unique(uni_plants))


y <- phylo.maker(uni_plants, scenarios = "S3")

# Using the phylogenetic tree in the ape R package to get a phylogenetic
# correlation matrix.
Cv_phylo <- round(ape::vcv(y[[1]], corr = TRUE), 10)

# Getting rid of the underscore in the names of species.
new_names_Cv <- sapply(rownames(Cv_phylo), function(r) gsub('_', ' ', r))
rownames(Cv_phylo) <- new_names_Cv
colnames(Cv_phylo) <- new_names_Cv

# Re-ordering the correlation matrix in the way that is in our data.
use_order <- sapply(as.character(uni_plants$species),
                    function(r) which(rownames(Cv_phylo) == r))
Cv_phylo <- Cv_phylo[use_order, use_order]

if (save_files) {
  save(Cv_phylo, file = paste0(save_path, 'Cv_phylo.dat'))
}



