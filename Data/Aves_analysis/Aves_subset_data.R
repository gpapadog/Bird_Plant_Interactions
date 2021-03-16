rm(list = ls())
dev.off()

setwd('~/Github/Birds_and_plants/Application/')

library(data.table)
library(rgdal)
library(ggmap)
library(ggplot2)
library(sf)
library(gplots)

save_files <- TRUE

# ---------- PART A: Loading in the data -------------- #

# Loading in the data set.

dta <- fread('Data/ATLANTIC_frugivory.csv')
wh <- which(dta$Longitude > 0)
dta$Longitude[wh] <- dta$Longitude[wh] * (- 1)
wh <- which(dta$Longitude > - 30)
dta$Longitude[wh] <- NA
dta$Latitude[wh] <- NA

# Excluding the Amaioua hybridus interaction since Amaioua hybridus should have been
# Amaranthus hybridus and that interaction already exists in the data.
# However, we are keeping the color.
subset(dta, Plant_Species %in% c('Amaioua hybridus', 'Amaranthus hybridus'))
dta <- subset(dta, Plant_Species != 'Amaioua hybridus')
dta$Fruit_color[dta$Plant_Species == 'Amaranthus hybridus'] <- 'red'


# ------ PART B: Subsetting the data to only bird-plant interactions -------- #

# Subsetting to birds:

birds <- subset(dta, Frug_Class == 'Aves')
dta_subset <- subset(dta, `Study reference` %in% birds$`Study reference`)

uni_birds <- unique(birds$Frugivore_Species)
uni_plants <- unique(dta_subset$Plant_Species)
nB <- length(uni_birds)
nP <- length(uni_plants)


# ------ PART C: Observed interaction and number of studies arrays -------- #

obs_A <- matrix(0, nrow = nB, ncol = nP, dimnames = list(uni_birds, uni_plants))
for (ss in 1 : nrow(birds)) {
  wh1 <- which(uni_birds == birds$Frugivore_Species[ss])
  wh2 <- which(uni_plants == birds$Plant_Species[ss])
  obs_A[wh1, wh2] <- 1
}
heatmap.2(obs_A, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
if (save_files) {
  save(obs_A, file = 'Data/Aves_analysis/obs_A.dat')
}


studies_birds <- lapply(uni_birds, function(x) unique(subset(dta_subset, Frugivore_Species == x)$`Study reference`))
studies_plants <- lapply(uni_plants, function(x) unique(subset(dta_subset, Plant_Species == x)$`Study reference`))
obs_n <- matrix(0, nrow = nB, ncol = nP, dimnames = list(uni_birds, uni_plants))
for (ii in 1 : nB) {
  for (jj in 1 : nP) {
    obs_n[ii, jj] <- length(intersect(studies_birds[[ii]], studies_plants[[jj]]))
  }
}
heatmap.2(obs_n, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
if (save_files) {
  save(obs_n, file = 'Data/Aves_analysis/obs_n.dat')
}


# ------ PART D: Covariates and correlation matrices -------- #

# All variables for the first set of species:

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
all_X <- unique(all_X)
nB
dim(all_X)

# Correlation matrix for the set of birds:

Cu <- diag(nB)
dimnames(Cu) <- list(uni_birds, uni_birds)

for (i1 in 1 : (nB - 1)) {
  for (i2 in (i1 + 1) : nB) {
    
    # Checking to make sure the order is the same.
    if (uni_birds[i1] != all_X$Frugivore_Species[i1]) print('Error')
    if (uni_birds[i2] != all_X$Frugivore_Species[i2]) print('Error')
    
    if (all_X$Frug_Genus[i1] == all_X$Frug_Genus[i2]) {
      Cu[i1, i2] <- Cu[i2, i1] <- 0.75
    } else if (all_X$Frug_Family[i1] == all_X$Frug_Family[i2]) {
      Cu[i1, i2] <- Cu[i2, i1] <- 0.5
    } else if (all_X$Frug_Order[i1] == all_X$Frug_Order[i2]) {
      Cu[i1, i2] <- Cu[i2, i1] <- 0.25
    }
    
  }
}

heatmap.2(Cu, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
if (save_files) {
  save(Cu, file = 'Data/Aves_analysis/Cu.dat')
}


# Covariate matrix of the first set of species:

all_X[, c('Frug_Class', 'Frug_Order', 'Frug_Family', 'Frug_Genus') := NULL]


# frugivory score” reliance on fruits. 1 (sporadic), 2 (moderate) and 3 (extensive frugivory)

# IUCN Red List Categories and Criteria are intended to be an easily and widely
# understood system for classifying species at high risk of global extinction.
# It divides species into nine categories: Not Evaluated, Data Deficient,
# Least Concern, Near Threatened, Vulnerable, Endangered, Critically Endangered,
# Extinct in the Wild and Extinct.

pB <- c(2, 3)
obs_X <- matrix(NA, nrow = nB, ncol = sum(pB))
rownames(obs_X) <- uni_birds
colnames(obs_X) <- c('logBodyMass', 'logGapeSize', 'isLarge', 'NotFruitDepend', 'isEndangered')
obs_X[, 1] <- log(all_X$Frug_Body_Mass)
obs_X[, 2] <- log(all_X$Frug_Mean_Gape_Size)
obs_X[, 3] <- (all_X$Frug_Group == 'Large birds') * 1
obs_X[, 4] <- (all_X$Frugivory_score == 1) * 1
obs_X[, 5] <- (all_X$Frug_IUCN != 'LC') * 1
if (save_files) {
  save(obs_X, file = 'Data/Aves_analysis/obs_X.dat')
}



# All variables for the second set of species:

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
all_W <- unique(all_W)
nP
dim(all_W)
apply(all_W, 2, function(x) mean(is.na(x)))


# Phylogenetic correlation for the second set of species:

Cv <- diag(nP)
dimnames(Cv) <- list(uni_plants, uni_plants)

for (i1 in 1 : (nP - 1)) {
  for (i2 in (i1 + 1) : nP) {
    
    # Checking to make sure the order is the same.
    if (uni_plants[i1] != all_W$Plant_Species[i1]) print('Error')
    if (uni_plants[i2] != all_W$Plant_Species[i2]) print('Error')
    
    if (all_W$Plant_genus[i1] == all_W$Plant_genus[i2]) {
      Cv[i1, i2] <- Cv[i2, i1] <- 0.5
    } else if (all_W$Plant_family[i1] == all_W$Plant_family[i2]) {
      Cv[i1, i2] <- Cv[i2, i1] <- 0.25
    }
    
  }
}

heatmap.2(Cv, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)
matrixcalc::is.positive.definite(Cv)
if (save_files) {
  save(Cv, file = 'Data/Aves_analysis/Cv.dat')
}



# Covariate information for the second set of species:
all_W[, c('Plant_family', 'Plant_genus') := NULL]

pP <- c(4, 8)
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
obs_W[all_W$Plants_IUCN == 'NE', 12] <- NA

colnames(obs_W) <- c('logFruitDiam', 'logFruitLen', 'LogSeedDiam', 'LogSeedLen',
                     'IsNative', 'IsTree', 'BlackFruit', 'RedFruit',
                     'YellowOrangeFruit', 'GreenFruit', 'LipidIs1',
                     'IsEndagered')

apply(obs_W, 2, function(x) mean(is.na(x)))
if (save_files) {
  save(obs_W, file = 'Data/Aves_analysis/obs_W.dat')
}



