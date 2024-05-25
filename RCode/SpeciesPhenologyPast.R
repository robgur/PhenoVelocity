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
library(lubridate)
library(distancetocoast)

#Assemble_past_phenology_data_species_by_species_from_Hough_1861
#example_assembly_for_one_species
#get unique latitude, longitude and year

#import_data
erythronium_americanum_flowers <- read_csv("C:/Users/robgu/Downloads/SpeciesPhenologyPast/erythronium_americanum_flowering_1850.csv",
                                           col_types = cols(`1851` = col_date(format = "%m/%d/%Y"),
                                                            `1852` = col_date(format = "%m/%d/%Y"),
                                                            `1853` = col_date(format = "%m/%d/%Y"),
                                                            `1854` = col_date(format = "%m/%d/%Y"),      
                                                            `1857` = col_date(format = "%m/%d/%Y"),
                                                            `1858` = col_date(format = "%m/%d/%Y"),
                                                            `1855` = col_date(format = "%m/%d/%Y"),
                                                            `1856` = col_date(format = "%m/%d/%Y"),
                                                            `1859` = col_date(format = "%m/%d/%Y")))

erythronium_americanum_foliage <- read_csv("C:/Users/robgu/Downloads/SpeciesPhenologyPast/erythronium_americanum_foliage_1850.csv",
                                    col_types = cols(`1851` = col_date(format = "%m/%d/%Y"),
                                                     `1852` = col_date(format = "%m/%d/%Y"),
                                                     `1853` = col_date(format = "%m/%d/%Y"),
                                                     `1854` = col_date(format = "%m/%d/%Y"),      
                                                     `1857` = col_date(format = "%m/%d/%Y"),
                                                     `1858` = col_date(format = "%m/%d/%Y"),
                                                     `1855` = col_date(format = "%m/%d/%Y"),
                                                     `1856` = col_date(format = "%m/%d/%Y"),
                                                     `1859` = col_date(format = "%m/%d/%Y")))


#pivot_longer_both_flowers_foliage
erythronium_americanum_flowers2  <- erythronium_americanum_flowers  %>% 
  pivot_longer(
    cols = starts_with("18"),
    names_to = "year",
    names_prefix = "yr",
    values_drop_na = TRUE
  )

erythronium_americanum_foliage2  <- erythronium_americanum_foliage  %>% 
  pivot_longer(
    cols = starts_with("18"),
    names_to = "year",
    names_prefix = "yr",
    values_drop_na = TRUE
  )

#add_state_variable
erythronium_americanum_flowers2$state <- rep("flowering", nrow(erythronium_americanum_flowers2))
erythronium_americanum_foliage2$state <- rep("foliage", nrow(erythronium_americanum_foliage2))

#bind_together_dataframes_remove_years_with_no_data
erythronium_americanum_all_1850s <- rbind(erythronium_americanum_flowers2,erythronium_americanum_foliage2)
erythronium_americanum_all_1850s3  <- subset(erythronium_americanum_all_1850s, year != "1851" &  year != "1856")
erythronium_americanum_all_1850s3  <- subset(erythronium_americanum_all_1850s, year != "1858")

#add_in_doy
erythronium_americanum_all_1850s3$doy <- yday(erythronium_americanum_all_1850s3$value)
erythronium_americanum_all_1850s3$all_species_dist <- extract(distance_to_coastline_10, cbind(erythronium_americanum_all_1850s3$longitude,erythronium_americanum_all_1850s3$latitude))

erythronium_americanum_all_1850s3_2 <- subset(erythronium_americanum_all_1850s3, longitude>-93)
erythronium_americanum_all_1850s3_3 <- subset(erythronium_americanum_all_1850s3_2, latitude>31 & latitude<46)

erythronium_americanum_all_1850s3_3$period <- rep("past", nrow(erythronium_americanum_all_1850s3_3))
erythronium_americanum_all_1850s3_3$taxon <- rep("erythronium_americanum", nrow(erythronium_americanum_all_1850s3_3))
erythronium_americanum_all_1850s3_3$inat <- rep("0", nrow(erythronium_americanum_all_1850s3_3))
erythronium_americanum_all_1850s3_done <- distinct(dplyr::select(erythronium_americanum_all_1850s3_3, latitude, longitude,year,period,taxon,doy,state,all_species_dist,inat))

write.csv(erythronium_americanum_all_1850s3_done, file="erythronium_americanum_all_1850s_done2.csv")


