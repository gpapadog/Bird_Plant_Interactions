# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# The directory where the original data are:
data_path <- 'Data/'
# The directory where the downloaded bird phylogeny data are.
# This directory should include the file named "AllBirdsEricson1.tre"
# The bird phylogeny data are downloaded from BirdTree.org. The file exists under
# phylogeny subsets -> download full trees -> birdtree -> Stage 2 ->
# EricsonStage2_0001_1000.zip
phylo_path <- 'Data/'
# The directory where the processed data should be saved:
save_path <- 'Data/'
# Whether the processed data should be saved or not:
save_files <- TRUE


# --------- BEGINNING -------- #

# Setting the working directory.
setwd(wd_path)

# Loading libraries.
library(ape)


# --------------------
# PART 1: Loading the bird phylogenetic trees from the donwloaded file.
# --------------------

x <- read.tree(paste0(phylo_path, 'AllBirdsEricson1.tre'))
# There are 1000 phylogenetic trees. Getting the name of the species included
# in those trees:
x_names <- sapply(x, function(y) y$tip.label)
table(apply(x_names, 1, function(y) length(unique(y))))
length(unique(as.character(x_names)))
# From these results we see that the 1000 phylogenetic trees do not have the
# species in the same order, even though they include the same species overall.


# --------------------
# PART 2: Loading the interaction data.
# --------------------

dta <- fread(paste0(data_path, 'ATLANTIC_frugivory.csv'))
ave_species <- unique(dta$Frugivore_Species[dta$Frug_Class == 'Aves'])
ave_species <- as.character(ave_species)

# Get rid of space in the end of the word.
ave_species <- stringr::str_replace(ave_species, "\\s+$", "")
# Replace spaces with underscore to match the names of the phylogenetic trees.
ave_species <- gsub(" ", "_", ave_species)


# --------------------
# PART 3: Checking the overlap of species in our data and in the phylogenies.
# --------------------

# Checking how many ave species from our interaction data exist in the
# phylogenetic trees.
sum(ave_species %in% x_names[, 1])

# Which ave species are missing:
miss_ave <- ave_species[which(!(ave_species %in% x_names[, 1]))]

# Some of these species are also known under alternative names in ecology.
# Specify these alternative names:
alt_names <- rep(NA, length(miss_ave))
alt_names[2] <- 'Pteroglossus_bailloni'
alt_names[3] <- 'Pipile_jacutinga'
alt_names[6] <- 'Basileuterus_flaveolus'
alt_names[8] <- 'Parula_pitiayumi'

cbind(miss_ave, alt_names)

alt_names[!is.na(alt_names)] %in% x_names[, 1]  # The new names exist in phylogeny
alt_names[!is.na(alt_names)] %in% ave_species  # No overlap with other names

# We notice that the remaining 10 entries that are missing from the phylogenies
# correspond to Ave names that are not species, rather than the genus that a
# species belongs to. Therefore, these recorded interactions did not specify
# which species from that genus formed the interaction. For this reason, we
# have to drop these from our analysis, leading us to 232 bird species. The
# names of the species to keep are included in the file below.
birds_232 <- setdiff(ave_species, miss_ave[is.na(alt_names)])
# Replacing the underscore back to a space:
birds_232 <- gsub("_", " ", birds_232)
# Saving the file so that we can easily subset the data to the species.
if (save_files) {
  save(birds_232, file = paste0(save_path, 'birds_232.dat'))
}


# --------------------
# PART 4: Getting a phylogenetic correlation matrix for each phylogenetic tree
# --------------------

# Use all 1000 samples from the phylogenetic tree
use_post_samples <- 1 : 1000

num_species <- length(ave_species) - sum(is.na(alt_names))
s <- array(NA, dim = c(num_species, num_species, length(use_post_samples)))

# For each phylogenetic tree from the 1000 posterior samples, calculate the
# phylogenetic correlation matrix using the ape R package.
t1 <- Sys.time()
for (ii in 1 : length(use_post_samples)) {
  if (ii %% 20 == 0) print(ii)
  this_sample <- use_post_samples[ii]
  # Getting the correlation matrix:
  this_cov <- ape::vcv(x[[this_sample]], corr = TRUE)
  # Re-ordering to get species on the same order:
  this_cov <- this_cov[order(rownames(this_cov)), order(rownames(this_cov))]
  # Keep only the species that are in our analysis:
  keep_birds <- which(rownames(this_cov) %in% c(ave_species, alt_names[!is.na(alt_names)]))
  this_cov <- this_cov[keep_birds, keep_birds]
  
  s[, , ii] <- this_cov
  if (ii == 1) {
    dimnames(s)[[1]] <- rownames(this_cov)
    dimnames(s)[[2]] <- colnames(this_cov)
  }
}
t2 <- Sys.time()
t2 - t1


# --------------------
# PART 5: Combining the correlation matrix in one phylogenetic matrix
# --------------------

# Take the mean correlation matrix across posterior samples and use this as the
# phylogenetic correlation matrix.
Cu_phylo <- apply(s, c(1, 2), mean)
gplots::heatmap.2(Cu_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)

# For the alternative names we used, put in the original names:
for (ii in 1 : length(alt_names)) {
  if (!is.na(alt_names[ii])) {
    
    old_name <- miss_ave[ii]
    wh_row <- which(rownames(Cu_phylo) == alt_names[ii])
    rownames(Cu_phylo)[wh_row] <- old_name
    colnames(Cu_phylo)[wh_row] <- old_name
    
  }
}

# Re-order the rows and columns to match the order of species in our data:
use_order <- sapply(ave_species, function(r) which(rownames(Cu_phylo) == r))
use_order <- use_order[sapply(use_order, length) > 0]
use_order <- as.numeric(use_order)
Cu_phylo <- Cu_phylo[use_order, use_order]
gplots::heatmap.2(Cu_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE,
                  labRow = FALSE, labCol = FALSE)


save(Cu_phylo, file = paste0(save_path, 'Cu_phylo.dat'))


