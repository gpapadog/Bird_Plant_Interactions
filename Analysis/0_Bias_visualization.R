# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the analysis is performed:
wd_path <- 'Bird_Plant_Interactions/'
# The directory where the original data are:
data_path <- 'Bird_Plant_Interactions/Data/'


library(data.table)
library(ggmap)
library(ggplot2)


# ---------- PART A: Loading in the data -------------- #

# Loading in the data set.

dta <- fread(paste0(data_path, 'ATLANTIC_frugivory.csv'))

# Correcting recorded coordinates. Some entries have recorded coordinates that
# are not correct as they are not in Brazil. These coordinates are set to NA.
wh <- which(dta$Longitude > 0)
dta$Longitude[wh] <- dta$Longitude[wh] * (- 1)
wh <- which(dta$Longitude > - 30)
dta$Longitude[wh] <- NA
dta$Latitude[wh] <- NA


# ------ PART B: Subsetting the data to only bird-plant interactions -------- #

# Keeping interactions that only include birds:
birds <- subset(dta, Frug_Class == 'Aves')
# Getting the locations of these interactions:
unique_locs <- unique(data.frame(long = birds$Longitude, lat = birds$Latitude))
unique_locs <- na.omit(unique_locs)


# --------------- PART C: Geographic Bias --------------- #

# Getting the brazilian map:
map_br <- get_stamenmap(bbox = c(left = min(unique_locs$long) - 1,
                                 bottom = min(unique_locs$lat) - 1, 
                                 right = max(unique_locs$long) + 1,
                                 top = max(unique_locs$lat) + 1), 
                        zoom = 5, maptype = 'terrain',
                        color = "color", force = TRUE)


# Plotting the locations of the bird and plant interactions to
# visualize geographic bias:
ggmap(map_br, extent = "device", legend = "right") + 
  geom_point(data = unique_locs, aes(x = long, y = lat),
             color = '#E86861', size = 1, alpha = 0.5)



# --------------- PART D: Taxonomic Bias --------------- #

# Distribution of number of plants, birds across studies of different types.

dta_by_type <- data.table::as.data.table(birds)
# Counting the number of unique birds and plants by study and study type:
dta_by_type <- dta_by_type[, list(unique_birds = length(unique(Frugivore_Species)) * 1,
                                  unique_plants = length(unique(Plant_Species)) * 1),
                           by = c('Study reference', 'Study_Method')]

# Reverting back to a data frame:
dta_by_type <- as.data.frame(dta_by_type)
# Re-arranging the data frame by plants and birds:
dta_by_type <- rbind(data.frame(Study_Method = dta_by_type$Study_Method, class = 'Aves',
                                number = dta_by_type$unique_birds),
                     data.frame(Study_Method = dta_by_type$Study_Method, class = 'Plants',
                                number = dta_by_type$unique_plants))

# Turning the study types to factor:
dta_by_type$Study_Method <- factor(as.character(dta_by_type$Study_Method),
                                   levels = c('Animal-oriented', 'Plant-oriented',
                                              'Network study'))

# Plotting: Visualization of taxonomic bias.
ggplot(data = dta_by_type) +
  geom_boxplot(aes(y = number, group = Study_Method, x = Study_Method)) +
  facet_wrap(~ class) +
  theme(legend.position = 'bottom') +
  xlab('') + ylab('Number of species in a study') +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 20, size = 10, vjust = 0.8))


