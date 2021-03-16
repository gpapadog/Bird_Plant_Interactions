# DESCRIPTION: Getting the order of species by family, genus etc to help
# meaningful visualization of results:

# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# The directory where the original data are:
data_path <- 'Data/'
# The directory where the processed data are and should be saved:
save_path <- 'Data/'
# Whether the processed data should be saved or not:
save_files <- TRUE


# --------- BEGINNING -------- #

library(ggplot2)
library(reshape2)
library(data.table)

setwd(wd_path)

# Loading in the data set.
dta <- fread(paste0(data_path, 'ATLANTIC_frugivory.csv'))
# Changing Amaioua hybridus to the correct family.
dta$Plant_family[dta$Plant_Species == 'Amaioua hybridus'] <- 'Rubiaceae'


# Loading in the correlation matrices based on the original data:
load(paste0(save_path, 'Cu.dat'))
load(paste0(save_path, 'Cv.dat'))

# Sample sizes
nB <- nrow(Cu)
nP <- ncol(Cv)



# ------------ GETTING THE ORDERING OF BIRD SPECIES ------------- #

# For each species of bird, get the list of bird indices that are in the same
# taxonomic genus, family or order (using the values in the phylogenetic
# correlation matrix)
#
bird_families <- NULL
for (ii in 1 : nB) {
  bird_families[[ii]] <- list(eq0 = as.numeric(which(Cu[ii, ] == 0)),
                              eq25 = as.numeric(which(Cu[ii, ] == 0.25)),
                              eq50 = as.numeric(which(Cu[ii, ] == 0.5)),
                              eq75 = as.numeric(which(Cu[ii, ] == 0.75)),
                              eq100 = as.numeric(which(Cu[ii, ] == 1)))
}

# testing that the clustering is done properly by ensuring that if i has
# species j in their taxonomic cluster then j also has i in the same type
# of taxonomic cluster.
for (ii in 1 : nB) {
  for (gg in 1 : 4) {
    for (other in bird_families[[ii]][[gg]]) {
      wh <- which(sapply(bird_families[[other]], function(x) ii %in% x))
      if (any(wh != gg)) {
        print('something wrong')
      }
    }
  }
}


# Cluster indices for grouping birds by order, family and genus. Each
# row corresponds to the different level of taxonomy. Clustering is performed
# such that species in the same family but different genus are given adjacent
# genus cluster neumbers. This is illustrated when plotting the re-ordered
# correlation matrix below:
#
groupings <- matrix(NA, nrow = 3, ncol = nB)
# Indices of current frugivore group by order family and genus.
curr_group1 <- 1
curr_group2 <- 1
curr_group3 <- 1


for (ii in 1 : nB) {
  
  if (is.na(groupings[1, ii])) {
    groupings[1, as.numeric(unlist(bird_families[[ii]][- 1]))] <- curr_group1
    
    for (jj in which(groupings[1, ] == curr_group1)) {
      if (is.na(groupings[2, jj])) {
        groupings[2, as.numeric(unlist(bird_families[[jj]][- c(1, 2)]))] <- curr_group2
        
        for (tt in which(groupings[2, ] == curr_group2)) {
          if (is.na(groupings[3, tt])) {
            groupings[3, as.numeric(unlist(bird_families[[tt]][- c(1, 2, 3)]))] <- curr_group3
            curr_group3 <- curr_group3 + 1
          }
        }
        curr_group2 <- curr_group2 + 1
      }
    }
    curr_group1 <- curr_group1 + 1
  }
}


# Getting the bird order by ordering their cluster number for genus. In the
# re-ordered data, species will be taxonomically clustered.
bird_order <- order(groupings[3, ])
if (save_files) {
  save(bird_order, file = paste0(save_path, 'bird_order.dat'))
}

# Plotting the re-ordered correlation matrix to make sure that our ordering
# was performed correctly:
ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = as.character(value)), 
                     data = reshape2::melt(Cu[bird_order, bird_order])) +
  theme(axis.text = element_blank()) +
  scale_fill_manual(values = c('white', "#F79694", "#F54D49", "#C23C3A", '#752523'))


# Saving the species' taxonomic information in the re-ordered data:
bird_order_info <- data.frame(Species = rownames(Cu)[bird_order])
family_info <- dta[, list(Frug_Order = Frug_Order[1], 
                          Frug_Family = Frug_Family[1],
                          Frug_Genus = Frug_Genus[1]), by = Frugivore_Species]
bird_order_info <- merge(bird_order_info, family_info, by.x = 'Species', by.y = 'Frugivore_Species', sort = FALSE)
if (save_files) {
  save(bird_order_info, file = paste0(save_path, 'bird_order_info.dat'))
}



# ------------ GETTING THE ORDERING OF PLANT SPECIES ------------- #

# We perform almost identical steps for plant species, though we note that
# plants are organized in genera and families only.

plant_families <- NULL
for (ii in 1 : nP) {
  plant_families[[ii]] <- list(eq0 = as.numeric(which(Cv[ii, ] == 0)),
                               eq25 = as.numeric(which(Cv[ii, ] == 0.25)),
                               eq50 = as.numeric(which(Cv[ii, ] == 0.5)),
                               eq100 = as.numeric(which(Cv[ii, ] == 1)))
}


# testing that the clustering is done properly:
for (ii in 1 : nP) {
  for (gg in 1 : 3) {
    for (other in plant_families[[ii]][[gg]]) {
      wh <- which(sapply(plant_families[[other]], function(x) ii %in% x))
      if (any(wh != gg)) {
        print('something wrong')
      }
    }
  }
}



groupings <- matrix(NA, nrow = 2, ncol = nP)
curr_group1 <- 1
curr_group2 <- 1

for (ii in 1 : nP) {
  
  if (is.na(groupings[1, ii])) {
    groupings[1, as.numeric(unlist(plant_families[[ii]][- 1]))] <- curr_group1
    
    for (jj in which(groupings[1, ] == curr_group1)) {
      if (is.na(groupings[2, jj])) {
        groupings[2, as.numeric(unlist(plant_families[[jj]][- c(1, 2)]))] <- curr_group2
        curr_group2 <- curr_group2 + 1
      }
    }
    curr_group1 <- curr_group1 + 1
  }
}

# Saving the plant order by genera:
plant_order <- order(groupings[2, ])
if (save_files) {
  save(plant_order, file = paste0(save_path, 'plant_order.dat'))
}

# Plotting the correlation matrix for the re-ordered species to ensure that
# clustering was performed correctly:
ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = as.character(value)), 
                     data = reshape2::melt(Cv[plant_order, plant_order])) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c('white', "#F79694", "#C23C3A", '#752523'))


# Saving taxonomic information for the species in the new order:
plant_order_info <- data.frame(Species = rownames(Cv)[plant_order])
family_info <- dta[, list(Plant_genus = Plant_genus[1], Plant_family = Plant_family[1]), by = Plant_Species]
plant_order_info <- merge(plant_order_info, family_info, by.x = 'Species', by.y = 'Plant_Species', sort = FALSE)

plant_order_info$Plant_genus <- as.factor(plant_order_info$Plant_genus)
plant_order_info$Plant_family <- as.factor(plant_order_info$Plant_family)

if (save_files) {
  save(plant_order_info, file = paste0(save_path, 'plant_order_info.dat'))
}
