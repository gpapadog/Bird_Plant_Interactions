rm(list = ls())
dev.off()

setwd('~/Documents/Research/Birds_and_plants/')

library(data.table)
library(rgdal)
library(ggmap)
library(ggplot2)
library(sf)

# ---------- PART A: Loading in the data -------------- #

# Loading in the data set.

dta <- fread('~/Github/Birds_and_plants/Application/Data/ATLANTIC_frugivory.csv')
wh <- which(dta$Longitude > 0)
dta$Longitude[wh] <- dta$Longitude[wh] * (- 1)
wh <- which(dta$Longitude > - 30)
dta$Longitude[wh] <- NA
dta$Latitude[wh] <- NA


# Loading in the map of Brazil:

# brazil <- rgdal::readOGR(dsn = "Data/Map_data/Level0/.")
# brazil@data$id <- brazil@data$ID_0
# brazil.points <- fortify(brazil, region = "id")
# brazil.df <- plyr::join(brazil.points, brazil@data, by = "id")
# all_brazil <- st_union(st_as_sf(brazil))



# ------ PART B: Subsetting the data to only bird-plant interactions -------- #

# Subsetting to birds:

birds <- subset(dta, Frug_Class == 'Aves')
apply(birds, 2, function(x) sum(is.na(x)))

dta_subset <- subset(dta, `Study reference` %in% birds$`Study reference`)

length(unique(dta_subset$Plant_Species))
length(unique(dta_subset$Frugivore_Species))

length(unique(birds$Plant_Species))
length(unique(birds$Frugivore_Species))

unique_locs <- unique(data.frame(long = birds$Longitude, lat = birds$Latitude))
unique_locs <- na.omit(unique_locs)


# Getting the brazilian map:

map_br <- get_stamenmap(bbox = c(left = min(unique_locs$long) - 1,
                                 bottom = min(unique_locs$lat) - 1, 
                                 right = max(unique_locs$long) + 1,
                                 top = max(unique_locs$lat) + 1), 
                        zoom = 5, maptype = 'terrain',
                        color = "color", force = TRUE)

ggmap(map_br, extent = "device", legend = "right") + 
  geom_point(data = unique_locs, aes(x = long, y = lat),
             color = '#E86861', size = 1, alpha = 0.5)

# 
# 
# ggplot() +
#   geom_sf(data = all_brazil, fill = 'transparent') +
#   geom_point(data = unique_locs, aes(x = long, y = lat),
#              color = '#E86861', size = 0.6, alpha = 1) 
# 



# Distribution of number of plants, birds across studies of different types.

dta_by_type <- data.table::as.data.table(birds)
dta_by_type <- dta_by_type[, list(unique_birds = length(unique(Frugivore_Species)) * 1,
                                  unique_plants = length(unique(Plant_Species)) * 1),
                           by = c('Study reference', 'Study_Method')]

ggplot(data = subset(dta_by_type, Study_Method != 'Network study')) +
  geom_histogram(aes(x = unique_birds / unique_plants, group = Study_Method, fill = Study_Method),
                 alpha = 1, bins = 5, color = 'black', position = position_dodge())


dta_by_type[, list(med_birds = median(unique_birds),
                   low_birds = quantile(unique_birds, 0.25),
                   low_birds = quantile(unique_birds, 0.75),
                   med_plants = median(unique_plants),
                   low_plants = quantile(unique_plants, 0.25),
                   high_plants = quantile(unique_plants, 0.75)),
            by = 'Study_Method']




dta_by_type2 <- as.data.frame(dta_by_type)
dta_by_type2 <- rbind(data.frame(Study_Method = dta_by_type2$Study_Method, class = 'Aves',
                                 number = dta_by_type2$unique_birds),
                      data.frame(Study_Method = dta_by_type2$Study_Method, class = 'Plants',
                                 number = dta_by_type2$unique_plants))


dta_by_type2$Study_Method <- factor(as.character(dta_by_type2$Study_Method),
                                    levels = c('Animal-oriented', 'Plant-oriented',
                                               'Network study'))


ggplot(data = dta_by_type2) +
  geom_boxplot(aes(y = number, group = Study_Method, x = Study_Method)) +
  facet_wrap(~ class) +
  theme(legend.position = 'bottom') +
  xlab('') + ylab('Number of species in a study') +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 20, size = 10, vjust = 0.8))



# ------ Looking at the missingness ----------- #

dta[, list(frug_spec = mean(is.na(Frugivore_Species)),
           plant_spec = mean(is.na(Plant_Species)),
           plant_fam = mean(is.na(Plant_family)),
           plant_epit = mean(is.na(Plant_specific.epiteth)),
           fruit_diam = mean(is.na(fruit_diameter)),
           fruit_len = mean(is.na(fruit_length)),
           seed_diam = mean(is.na(seed_diameter)),
           seed_len = mean(is.na(seed_length)),
           friut_color = mean(is.na(Fruit_color)),
           lipid = mean(is.na(Lipid_Score)),
           frug_epi = mean(is.na(Frug_Epitetus)),
           frug_mass = mean(is.na(Frug_Body_Mass)),
           frug_gape = mean(is.na(Frug_Mean_Gape_Size)),
           frug_Score = mean(is.na(Frugivory_score)),
           frug_mig = mean(is.na(Frug_Migration_status))),
    by = 'Frug_Class']


