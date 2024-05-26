library(tidyr)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(sp)
library(lubridate)
library(phenesse)
library(terra)
library(ggplot2)
library(raster)
library(readr)
library(performance)
library(lmerTest)
library(lme4)
library(readr)
library(AOI)
library(climateR)
library(stringr)
library(climateR)
library(AOI)
library(sjPlot)

#read_in_NPN_data
data <- read_csv("C:/Users/robgu/Downloads/SpeciesPhenologyPresent/carya_alba_present_flowering.csv")
data2 <- read_csv("C:/Users/robgu/Downloads/SpeciesPhenologyPresent/carya_alba_present_foliage.csv")

#combine_data_and_deal_with_dates
carya_alba_present <- rbind(data,data2)
carya_alba_present$ObsDate <- mdy(carya_alba_present$Observation_Date)
carya_alba_present$year <- year(carya_alba_present$ObsDate)
carya_alba_present$doy <- yday(carya_alba_present$ObsDate)

#group_by_site_state_and_year_get_earliest_day
carya_alba_present_2 <- carya_alba_present %>%
  group_by(Longitude, Latitude, state, year) %>%
  slice(which.min(doy))
#same_group_and_get_count
carya_alba_present_count <- carya_alba_present %>%
  group_by(Longitude, Latitude, state,year) %>%
  dplyr::summarize(count=n())

#merge_counts_and_min_day
carya_alba_present_3  <- merge(carya_alba_present_2, 
                               carya_alba_present_count, by=c("Longitude","Latitude", "state","year"))

#clean_up_outliers_and_filter_counts_less_than_3
carya_alba_present4 <- subset(carya_alba_present_3, Longitude>-100 & Longitude<0)
carya_alba_present6 <- subset(carya_alba_present4, doy<176)
carya_alba_present7 <- subset(carya_alba_present6, doy>15)
carya_alba_present8 <- subset(carya_alba_present7, count>2)

#remove_years_with_missing_Data_standardize_field_names_and_remove_unnecessary_fields
carya_alba_present8 <- subset(carya_alba_present8, year != "2010" & year != "2012" & year != "2013" &  year != "2014" &  year != "2016" &  year != "2020" &  year != "2023") 
carya_alba_present81  <- dplyr::rename(carya_alba_present8 , latitude = Latitude)
carya_alba_present82 <-  dplyr::rename(carya_alba_present81, longitude = Longitude)
carya_alba_present83  <- unite(carya_alba_present82, col='taxon', c('Genus', 'Species'), sep=' ')
carya_alba_present83$all_species_dist <- extract(distance_to_coastline_10, cbind(carya_alba_present83$longitude,carya_alba_present83$latitude))
carya_alba_present84  <- carya_alba_present83 %>% dplyr::select(latitude, longitude,period,taxon,doy,state,year,all_species_dist)

#write_present_phenometrics
write.csv(carya_alba_present84, file="carya_alba_present84.csv")
