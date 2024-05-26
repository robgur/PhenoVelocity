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
library(prism)
library(stringr)
library(fs)
library(distancetocoast)

#function_for_downloading_iNat_images
downlaod_images <- function(dat,size="medium",outpath="."){
  for (i in 1:dim(dat)[1]){
    iurl <- dat$image_url[i]
    iurl <- gsub("medium",size,iurl)
    iname <- paste(outpath,dat$id[i]," ",dat$scientific_name[i]," ",dat$observed_on[i],".jpg",sep="")
    tryCatch(
      {download.file(iurl,iname, mode = 'wb')},
      error = function(err) {print(paste("MY_ERROR:  ",err))}
    )
    Sys.sleep(2)
  }
}

#import_metadata_including_iNat_metadata_and_annotations
observations_390429 <- read_csv("C:/Users/robgu/Downloads/erythronium_americanum_metadata.csv")
erythronium_americanum_flowers <- read.csv("C:/Users/robgu/Downloads/erythronium_americanum.csv")
erythronium_americanum_flowers$file2 <- path_ext_remove(erythronium_americanum_flowers$file)
erythronium_americanum_flowers$file2  <- as.numeric(erythronium_americanum_flowers$file2)
erythronium_americanum_flowerj <- left_join(observations_390429,erythronium_americanum_flowers, by=c('ID'='file2'))
erythronium_americanum_flower_filter <- erythronium_americanum_flowerj %>%
  dplyr::select(Longitude, Latitude, everything()) %>%
  filter(leaves == 1) %>%
  tidyr::drop_na(Longitude, Latitude)
erythronium_americanum_flower_filter <- erythronium_americanum_flower_filter %>%   rename(latitude = Latitude, longitude=Longitude, observed_on = Observed_On, id=ID) 



#get_region_prject_to_albers_equal_area
test <- rnaturalearth::ne_states(country = c("United States of America"),returnclass = "sf") %>%
  filter(name %in% c("Pennsylvania", "Delaware", "Maryland", "District of Columbia", "Michigan", "Ohio", "Indiana", "Iowa", "Illinois", "Wisconsin", "Minnesota", "West Virginia", "Virginia", "North Carolina", "Tennessee", "Kentucky", "South Carolina", "Georgia", "Alabama", "Missouri", "Mississippi", "Florida"))
test2<-st_transform(test, sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
test3 <-st_union(test2)

#make_grid_over_region
grids <- st_make_grid(test3, cellsize = c(55000, 55000))
grids2 <- st_join(st_as_sf(grids), st_as_sf(test3))
grids2 = dplyr::mutate(st_sf(geometry = grids2), id_cells = 1:n())
grids3 <- grids2[test3, col = '#ff000088']

#intersect points to grids
erythronium_americanum_flower_sf <-erythronium_americanum_flower_filter  %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "+proj=longlat +datum=WGS84 +no_defs")
erythronium_americanum_flower_sf <- st_transform(erythronium_americanum_flower_sf,crs =st_crs(grids))
erythronium_americanum_flower_sf <- st_join(erythronium_americanum_flower_sf, grids3)
erythronium_americanum_flower_df <- st_drop_geometry(erythronium_americanum_flower_sf)


#get_year_filter_dates_outside_range
erythronium_americanum_flower_df$observed_on <-anytime::anydate(erythronium_americanum_flower_df$observed_on)
erythronium_americanum_flower_df$year <- year(erythronium_americanum_flower_df$observed_on)
erythronium_americanum_flower_df$doy <- yday(erythronium_americanum_flower_df$observed_on)
erythronium_americanum_flower_df2 <-erythronium_americanum_flower_df %>%
  filter(doy > 10 & doy < 225 )
erythronium_americanum_flower_df2 <- distinct(erythronium_americanum_flower_df2, doy,id_cells,id, .keep_all= TRUE)

#check_distinctdays_ceiling_recsperday_and_get_counts
erythronium_americanum_flower_ndistinct <-erythronium_americanum_flower_df2 %>%
  group_by(year, id_cells) %>%
  dplyr::summarise(ndistinct = n_distinct(doy))

erythronium_americanum_flower_df3 <- erythronium_americanum_flower_df2 %>%
  group_by(id_cells,year,doy) %>%  slice_sample(n = 5)

erythronium_americanum_flower_count <-erythronium_americanum_flower_df3 %>%
  group_by(year, id_cells) %>%
  dplyr::summarise(obsnum = n())  

erythronium_americanum_flower_count2 <-erythronium_americanum_flower_count %>% filter(obsnum >= 7)
erythronium_americanum_flower_count2 <- erythronium_americanum_flower_count2 %>% drop_na(id_cells)
##################################################################

######################################################
# Quantile_estimates_get_onset_mean_offset_although_target_here_just_onset
erythronium_americanum_flower_mean_quantile <-erythronium_americanum_flower_df3 %>%
  group_by(year, id_cells) %>%
  filter(!is.na(id_cells)) %>%
  filter(n() > 7)  %>%
  filter(n_distinct(doy)>4) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.5, bootstraps = 500)))
erythronium_americanum_flower_onset_quantile <-erythronium_americanum_flower_df3 %>%
  group_by(year, id_cells) %>%
  filter(!is.na(id_cells)) %>%
  filter(n() > 7)  %>%
  filter(n_distinct(doy)>4) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.05, bootstraps = 500)))
erythronium_americanum_flower_offset_quantile <-erythronium_americanum_flower_df3 %>%
  group_by(year, id_cells) %>%
  filter(!is.na(id_cells)) %>%
  filter(n() > 7)  %>%
  filter(n_distinct(doy)>4) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.95,bootstraps=500)))

#cleanup_phenometrics_and_assemble_into_formats_compatible_with_other_estimates
erythronium_americanum_flower_mean_quantile2 <- dplyr::select(erythronium_americanum_flower_mean_quantile, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
erythronium_americanum_flower_mean_quantile3 <- pivot_wider(erythronium_americanum_flower_mean_quantile2, names_from = column, values_from = mean)
erythronium_americanum_flower_onset_quantile2 <- dplyr::select(erythronium_americanum_flower_onset_quantile, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
erythronium_americanum_flower_onset_quantile3 <- pivot_wider(erythronium_americanum_flower_onset_quantile2, names_from = column, values_from = mean)
erythronium_americanum_flower_offset_quantile2 <- dplyr::select(erythronium_americanum_flower_offset_quantile, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
erythronium_americanum_flower_offset_quantile3 <- pivot_wider(erythronium_americanum_flower_offset_quantile2, names_from = column, values_from = mean)
erythronium_americanum_year2_phen_onoff_quan <- merge(erythronium_americanum_flower_onset_quantile3,erythronium_americanum_flower_offset_quantile3, by=c("year", "id_cells"))
erythronium_americanum_year2_phen_onoffmean_quan <- merge(erythronium_americanum_flower_mean_quantile3,erythronium_americanum_year2_phen_onoff_quan, by=c("year", "id_cells"))
erythronium_americanum_year2_phen_all_quan <- merge(erythronium_americanum_flower_count,erythronium_americanum_year2_phen_onoffmean_quan, by=c("year", "id_cells"))
erythronium_americanum_year2_phen_all2_quan <- merge(erythronium_americanum_flower_ndistinct,erythronium_americanum_year2_phen_all_quan, by=c("year", "id_cells"))
erythronium_americanum_year2_phen_all2_quan$species <-  rep("erythronium_americanum",nrow(erythronium_americanum_year2_phen_all2_quan))
erythronium_americanum_year2_phend_quan<-erythronium_americanum_year2_phen_all2_quan %>%
  dplyr::rename(mean=estimate,mean_low=low_ci, mean_high=high_ci, onset=estimate.x,onset_low=low_ci.x, onset_high=high_ci.x, offset =estimate.y,offset_low=low_ci.y, offset_high=high_ci.y)

#get_latitude_longitude_centroid_for_gridcells_and_join_back_to_results
grids4 <- st_transform(grids3,crs("EPSG:4326"))
grid_cent <- st_centroid(grids4)
grid_wcent <- dplyr::mutate(grids3,st_coordinates(grid_cent))
grid_wcent2 <-st_drop_geometry(grid_wcent)
test <- as.matrix(grid_wcent2)
cell_latitude_long <- cbind(test[,1],test[,2],test[,3])
cell_latlon2 <- as.data.frame(cell_latitude_long)
cell_latlon3 <- cell_latlon2 %>% dplyr::rename(id_cells=V1, longitude=V2, latitude=V3)
erythronium_americanum_year2_phend_quan2 <- left_join(erythronium_americanum_year2_phend_quan, cell_latlon3 , by="id_cells")

erythronium_americanum_year2_phend_quan2 <- dplyr::rename(erythronium_americanum_year2_phend_quan2, doy=onset)
erythronium_americanum_year2_phend_quan2$inat <- rep("1", nrow(erythronium_americanum_year2_phend_quan2))
erythronium_americanum_year2_phend_quan2$period <- rep("present",nrow(erythronium_americanum_year2_phend_quan2))
erythronium_americanum_year2_phend_quan2$taxon <- rep("erythronium_americanum", nrow(erythronium_americanum_year2_phend_quan2))
erythronium_americanum_year2_phend_quan2$state <- rep("foliage",nrow(erythronium_americanum_year2_phend_quan2))
erythronium_americanum_year2_phend_quan2$all_species_dist <- extract(distance_to_coastline_10, cbind(erythronium_americanum_year2_phend_quan2$longitude,erythronium_americanum_year2_phend_quan2$latitude))
erythronium_americanum_year2_phend_quan3 <- distinct(dplyr::select(erythronium_americanum_year2_phend_quan2, latitude, longitude,year,period,taxon,doy,state,all_species_dist,inat))

write.csv(erythronium_americanum_year2_phend_quan, "/Users/robgu/OneDrive/Desktop/image_data/phenology_results/quantile_erythronium_americanum.csv")

