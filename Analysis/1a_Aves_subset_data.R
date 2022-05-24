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
library(data.table)
library(gplots)

# ---------- PART A: Loading in the data -------------- #


# Loading in the data set.
dta <- fread(paste0(data_path, 'ATLANTIC_frugivory.csv'))
# Studies with the same reference but different method are treated separately:
dta$Study_combined <- sapply(1 : nrow(dta), function(r) paste(dta$`Study reference`[r], '-', dta$Study_Method[r]))

# Correcting recorded coordinates. Some entries have recorded coordinates that
# are not correct as they are not in Brazil. These coordinates are set to NA.
wh <- which(dta$Longitude > 0)
dta$Longitude[wh] <- dta$Longitude[wh] * (- 1)
wh <- which(dta$Longitude > - 30)
dta$Longitude[wh] <- NA
dta$Latitude[wh] <- NA


# Correcting a wrong recorded interaction.
#
# After communications with the data manager, we realized that the species that
# was originally recorded as "Amaioua hybridus" does not exist, and it should
# have been recorded as "Amaranthus hubridus".
#
# We excluding the Amaioua hybridus interaction since an interactions between
# Amaranthus hybridus (the correct species) and the frugivore species already
# exists in the data.
#
# However, we are keeping the plant color (plant covariate).
#
subset(dta, Plant_Species %in% c('Amaioua hybridus', 'Amaranthus hybridus'))
dta <- subset(dta, Plant_Species != 'Amaioua hybridus')
dta$Fruit_color[dta$Plant_Species == 'Amaranthus hybridus'] <- 'red'



# ------ PART B: Subsetting the data to only bird-plant interactions -------- #

# We are interested in studying bird and plant interactions. The original data
# include interactions among five frugivore classes. Aves is the class of
# birds.
# We subset our data to studies that included bird-plant interactions. However,
# we want to keep in our data all plant species for studies that included
# birds, even if those plants were never recorded to interact with a bird. That
# is because we are interested in studying whether an unrecorded interaction is
# truly possible. A study that followed no birds could not record any
# bird-plant interaction so studies like that are excluded. However, a study
# that followed a bird and did not record an interaction with a certain plant
# is informative of whether that given plant interacts with the birds in the
# study.

# What interactions have been recorded that include bird species:
birds <- subset(dta, Frug_Class == 'Aves')
# Include all species and interactions for studies that followed at least one
# bird:
dta_subset <- subset(dta, `Study reference` %in% birds$`Study reference`)

# Keeping track of the unique bird and plant species in the study, and their
# number. Also keeping track of the unique studies.
uni_birds <- unique(birds$Frugivore_Species)
uni_plants <- unique(dta_subset$Plant_Species)
uni_studies <- unique(birds$Study_combined)
nB <- length(uni_birds)
nP <- length(uni_plants)
nS <- length(uni_studies)
cat(nB, nP, nS)
# 85 unique studies, but 2 of them are both species oriented and network
# studies, so that amounts to 87 networks.


# ------ PART C: Array of observed interactions -------- #

# Creating an array with dimensions that correspond to birds, plants and
# studies. Entries are equal to 0 if the species did not interact in the
# specific study, and equal to 1 if they were recorded to interact.
#
obs_A <- array(0, dim = c(nB, nP, nS))
dimnames(obs_A) <- list(uni_birds, uni_plants, uni_studies)
for (ss in 1 : nrow(birds)) {
  wh1 <- which(uni_birds == birds$Frugivore_Species[ss])
  wh2 <- which(uni_plants == birds$Plant_Species[ss])
  wh3 <- which(uni_studies == birds$Study_combined[ss])
  obs_A[wh1, wh2, wh3] <- 1
}
# Plotting a heatmap of the recorded interactions for a single study,
# and across all studies:
heatmap.2(obs_A[, , 1], dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
heatmap.2(apply(obs_A, c(1, 2), function(x) any(x == 1) * 1), dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
# Saving the binary array of recorded interactions.
if (save_files) {
  save(obs_A, file = paste0(save_path, 'obs_A.dat'))
}


# ------ PART D: Covariates and correlation matrices -------- #

# Here we re-format the original data to acquire matrices of covariate
# information for each set of species, and the phylogenetic correlation
# matrices.

# Starting for the first set of species.

# Including the bird covariate information in one data frame. This data frame
# includes phylogenetic information that will be used to create the species'
# correlation matrix.
#
all_X <- birds[, list(Frugivore_Species = Frugivore_Species,
                      Frug_Class = Frug_Class,
                      Frug_Order = Frug_Order, 
                      Frug_Family = Frug_Family,
                      Frug_Genus = Frug_Genus,
                      Frug_Group = Frug_Group,
                      Frug_Body_Mass = Frug_Body_Mass,
                      Frug_Mean_Gape_Size = Frug_Mean_Gape_Size,
                      Frugivory_score = Frugivory_score,
                      Frug_Migration_status = Frug_Migration_status,
                      Frug_IUCN = Frug_IUCN,
                      Frug_Population_Trend = Frug_Population_Trend)]
# Excluding replicates of birds:
all_X <- unique(all_X)
# Ensuring that the dimensions are correct:
cat(nB, dim(all_X)[1])
# Ensuring that the species in all_X are in the same order as in obs_A:
sum(all_X$Frugivore_Species != dimnames(obs_A)[[1]])  # Should be equal to 0.


# Creating the correlation matrix for the set of birds based on TAXONOMY.
# This correlation matrix is used in 1b_get_order.R to acquire the order
# of species in our data based on their taxonomic order.
#
# The phylogeny-based correlation matrix is created in a separate file.
# We use the phylogeny-based correlation matrix in the manuscript for the
# analysis of the interaction data. So this taxonomic correlation matrix is
# not used directly for our analysis.
#
# If the birds are in the same genus, family or order, the correlation is set
# to 0.75, 0.5, and 0.25, respectively.
#
Cu <- diag(nB)
dimnames(Cu) <- list(uni_birds, uni_birds)

for (i1 in 1 : (nB - 1)) {
  for (i2 in (i1 + 1) : nB) {
    
    if (all_X$Frug_Genus[i1] == all_X$Frug_Genus[i2]) {
      Cu[i1, i2] <- Cu[i2, i1] <- 0.75
    } else if (all_X$Frug_Family[i1] == all_X$Frug_Family[i2]) {
      Cu[i1, i2] <- Cu[i2, i1] <- 0.5
    } else if (all_X$Frug_Order[i1] == all_X$Frug_Order[i2]) {
      Cu[i1, i2] <- Cu[i2, i1] <- 0.25
    }
    
  }
}

# Visualizing it.
heatmap.2(Cu, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
# Saving it to our data information:
if (save_files) {
  save(Cu, file = paste0(save_path, 'Cu_tax.dat'))
}


# Excluding phylogenetic information from the matrix of covariates:
all_X[, c('Frug_Class', 'Frug_Order', 'Frug_Family', 'Frug_Genus') := NULL]

# Some information on the covariates:
#
# frugivory score” reliance on fruits. 1 (sporadic), 2 (moderate) and 3 (extensive frugivory)
#
# IUCN Red List Categories and Criteria are intended to be an easily and widely
# understood system for classifying species at high risk of global extinction.
# It divides species into nine categories: Not Evaluated, Data Deficient,
# Least Concern, Near Threatened, Vulnerable, Endangered, Critically Endangered,
# Extinct in the Wild and Extinct.

pB <- c(2, 3)  # Two continous and three binary covariates:
obs_X <- matrix(NA, nrow = nB, ncol = sum(pB))
rownames(obs_X) <- uni_birds
colnames(obs_X) <- c('logBodyMass', 'logGapeSize', 'isLarge', 'NotFruitDepend', 'isEndangered')
obs_X[, 1] <- log(all_X$Frug_Body_Mass)
obs_X[, 2] <- log(all_X$Frug_Mean_Gape_Size)
obs_X[, 3] <- (all_X$Frug_Group == 'Large birds') * 1
obs_X[, 4] <- (all_X$Frugivory_score == 1) * 1
obs_X[, 5] <- (all_X$Frug_IUCN != 'LC') * 1
if (save_files) {
  save(obs_X, file = paste0(save_path, 'obs_X.dat'))
}



# For the second set of species:
# 
# We perform a similar process to acquire their phylogenetic correlation matrix
# and their covariate information.


# All covariates and phylogenetic information in one data frame:
all_W <- dta_subset[, list(Plant_Species = Plant_Species,
                           Plant_family = Plant_family,
                           Plant_genus = Plant_genus,
                           Plant_distribution = Plant_distribution,
                           Plant_origin = Plant_origin,
                           Plant_Form = Plant_Form,
                           fruit_diameter = fruit_diameter,
                           fruit_length = fruit_length,
                           seed_diameter = seed_diameter,
                           seed_length = seed_length,
                           Fruit_color = Fruit_color,
                           Lipid_Score = Lipid_Score,
                           Plants_IUCN = Plants_IUCN)]
# Excluding species that are duplicated:
all_W <- unique(all_W)
# Making sure that the dimension and order of species is correct:
cat(nP, dim(all_W)[1])
sum(all_W$Plant_Species != dimnames(obs_A)[[2]])  # Should be 0.


# Creating the taxonomic correlation matrix for the second set of species.
# The correlation is 0.5 and 0.25 for the same genus and family, respectively.
# Again, the phylogenetic correlation matrix is created in a separate file.
#
Cv <- diag(nP)
dimnames(Cv) <- list(uni_plants, uni_plants)

for (i1 in 1 : (nP - 1)) {
  for (i2 in (i1 + 1) : nP) {
    
    if (all_W$Plant_genus[i1] == all_W$Plant_genus[i2]) {
      Cv[i1, i2] <- Cv[i2, i1] <- 0.5
    } else if (all_W$Plant_family[i1] == all_W$Plant_family[i2]) {
      Cv[i1, i2] <- Cv[i2, i1] <- 0.25
    }
    
  }
}

heatmap.2(Cv, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
if (save_files) {
  save(Cv, file = paste0(save_path, 'Cv_tax.dat'))
}



# Exclude phylogenetic information (we don't need it anymore).
all_W[, c('Plant_family', 'Plant_genus') := NULL]

pP <- c(4, 8)  # Number of continuous and binary covariates.
obs_W <- matrix(NA, nrow = nP, ncol = sum(pP))

# origin: native or not
# form: tree or not
# color: black, red, green
# lipid: 1 or other

obs_W[, 1] <- log(all_W$fruit_diameter)
obs_W[, 2] <- log(all_W$fruit_length)
obs_W[, 3] <- log(all_W$seed_diameter)
obs_W[, 4] <- log(all_W$seed_length)
obs_W[, 5] <- (all_W$Plant_origin == 'native') * 1
obs_W[, 6] <- (all_W$Plant_Form == 'tree') * 1
obs_W[, 7] <- (all_W$Fruit_color == 'black') * 1
obs_W[, 8] <- (all_W$Fruit_color == 'red') * 1
obs_W[, 9] <- (all_W$Fruit_color == 'orange') + (all_W$Fruit_color == 'yellow')
obs_W[, 10] <- (all_W$Fruit_color == 'green') * 1
obs_W[, 11] <- (all_W$Lipid_Score == 1) * 1
obs_W[, 12] <- (all_W$Plants_IUCN != 'LC') * 1
obs_W[all_W$Plants_IUCN == 'NE', 12] <- NA  # 'NE' means missing. Setting to NA

# Creating the column names.
colnames(obs_W) <- c('logFruitDiam', 'logFruitLen', 'LogSeedDiam', 'LogSeedLen',
                     'IsNative', 'IsTree', 'BlackFruit', 'RedFruit',
                     'YellowOrangeFruit', 'GreenFruit', 'LipidIs1',
                     'IsEndagered')

if (save_files) {
  save(obs_W, file = paste0(save_path, 'obs_W.dat'))
}




# ------ PART E: Focus and species occurrence across studies -------- #


# Creating two arrays that correspond to whether the studies would have
# recorded an observed interaction or not, and whether species co-occur in the
# study area. We will assume that these arrays are 0/1. For the focus array,
# this makes sense as we have access to whether the study is animal or plant
# focused or a network study. For the species occurrence, we know for a fact
# that a species with a recorded interaction occurs in the area. The only
# simplifying assumption is therefore when we set probability of occurrence to
# zero for species without any recorded interaction.
#

# First for study focus:
obs_F <- array(0, dim = c(nB, nP, nS))
dimnames(obs_F) <- list(uni_birds, uni_plants, uni_studies)

for (ss in 1 : nS) {
  wh_study <- uni_studies[ss]
  study_method <- unique(birds$Study_Method[birds$Study_combined == wh_study])
  if (length(study_method) > 1) {
    cat('Error. Each study should have a single method.')
  }
  # If this includes a network study, then create the most general focus matrix.
  if (study_method == 'Network study') {
    obs_F[, , ss] <- 1
  } else if (study_method == 'Animal-oriented') {
    these_birds <- unique(birds$Frugivore_Species[birds$Study_combined == wh_study])
    wh_birds <- which(uni_birds %in% these_birds)
    obs_F[wh_birds, , ss] <- 1
  } else if (study_method == 'Plant-oriented') {
    these_plants <- unique(birds$Plant_Species[birds$Study_combined == wh_study])
    wh_plants <- which(uni_plants %in% these_plants)
    obs_F[, wh_plants, ss] <- 1
  } else {
    cat('Error')
  }
}

if (save_files) {
  save(obs_F, file = paste0(save_path, 'obs_F.dat'))
}


# For occurrence of network study
obs_OB <- array(0, dim = c(nB, nS))
obs_OP <- array(0, dim = c(nP, nS))
dimnames(obs_OB) <- list(uni_birds, uni_studies)
dimnames(obs_OP) <- list(uni_plants, uni_studies)

for (ss in 1 : nS) {
  wh_study <- uni_studies[ss]
  these_birds <- unique(birds$Frugivore_Species[birds$Study_combined == wh_study])
  these_plants <- unique(birds$Plant_Species[birds$Study_combined == wh_study])
  wh_birds <- which(uni_birds %in% these_birds)
  wh_plants <- which(uni_plants %in% these_plants)
  obs_OB[wh_birds, ss] <- 1
  obs_OP[wh_plants, ss] <- 1
}

if (save_files) {
  save(obs_OB, file = paste0(save_path, 'obs_OB.dat'))
  save(obs_OP, file = paste0(save_path, 'obs_OP.dat'))
}




