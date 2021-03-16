# Getting the species order by family, genus etc to help visualization of results:

# Plotting the results of the data analysis.

library(ggplot2)
library(reshape2)
library(data.table)
library(superheat)

save_files <- FALSE

setwd('~/Documents/Research/Birds_and_plants/')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/Cu.dat')
load('~/Github/Birds_and_plants/Application/Data/Aves_analysis/Cv.dat')

dta <- fread('~/Github/Birds_and_plants/Application/Data/ATLANTIC_frugivory.csv')
# Changing Amaioua hybridus to the correct family:
dta$Plant_family[dta$Plant_Species == 'Amaioua hybridus'] <- 'Rubiaceae'

nB <- nrow(Cu)
nP <- ncol(Cv)


# ------------ GETTING THE ORDERING OF BIRD SPECIES ------------- #

bird_families <- NULL
for (ii in 1 : nB) {
  bird_families[[ii]] <- list(eq0 = as.numeric(which(Cu[ii, ] == 0)),
                              eq25 = as.numeric(which(Cu[ii, ] == 0.25)),
                              eq50 = as.numeric(which(Cu[ii, ] == 0.5)),
                              eq75 = as.numeric(which(Cu[ii, ] == 0.75)),
                              eq100 = as.numeric(which(Cu[ii, ] == 1)))
}

# testing that the clustering is done properly:
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


groupings <- matrix(NA, nrow = 4, ncol = nB)
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


bird_order <- order(groupings[3, ])
if (save_files) {
  save(bird_order, file = '~/Github/Birds_and_plants/Application/Data/Aves_analysis/bird_order.dat')
}

ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = as.character(value)), 
                     data = reshape2::melt(Cu[bird_order, bird_order])) +
  theme(axis.text = element_blank()) +
  scale_fill_manual(values = c('white', "#F79694", "#F54D49", "#C23C3A", '#752523'))



bird_order_info <- data.frame(Species = rownames(Cu)[bird_order])
family_info <- dta[, list(Frug_Order = Frug_Order[1], 
                          Frug_Family = Frug_Family[1],
                          Frug_Genus = Frug_Genus[1]), by = Frugivore_Species]
bird_order_info <- merge(bird_order_info, family_info, by.x = 'Species', by.y = 'Frugivore_Species', sort = FALSE)
if (save_files) {
  save(bird_order_info, file = '~/Github/Birds_and_plants/Application/Data/Aves_analysis/bird_order_info.dat')
}


# ------------ GETTING THE ORDERING OF PLANT SPECIES ------------- #

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


plant_order <- order(groupings[2, ])

if (save_files) {
  save(plant_order, file = '~/Github/Birds_and_plants/Application/Data/Aves_analysis/plant_order.dat')
}

ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = as.character(value)), 
                     data = reshape2::melt(Cv[plant_order, plant_order])) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c('white', "#F79694", "#C23C3A", '#752523'))





plant_order_info <- data.frame(Species = rownames(Cv)[plant_order])
family_info <- dta[, list(Plant_genus = Plant_genus[1], Plant_family = Plant_family[1]), by = Plant_Species]
plant_order_info <- merge(plant_order_info, family_info, by.x = 'Species', by.y = 'Plant_Species', sort = FALSE)

plant_order_info$Plant_genus <- as.factor(plant_order_info$Plant_genus)
plant_order_info$Plant_family <- as.factor(plant_order_info$Plant_family)

if (save_files) {
  save(plant_order_info, file = '~/Github/Birds_and_plants/Application/Data/Aves_analysis/plant_order_info.dat')
}
