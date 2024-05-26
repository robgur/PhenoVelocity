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
library(car)


#get_climate_for_one_species
setwd("../../Downloads")
all_species_pastpresent <- read_csv("C:/Users/robgu/OneDrive/Desktop/onset_pastpresent/all_species_pastpresent_n_outlier_5.csv")

########################################################################### past_data_assembly_using_Berkeley_Earth
#data_for_past_downloaded_from_Berkeley_Earth
#need_both_temperature_anomaly_and_climatology_for_this_dataset
tmp_ncdf <- brick("CONUS_TAVG_Gridded_0p25.nc",varname="temperature")
tmp_ncdf <- as(tmp_ncdf, "SpatRaster")
clim_ncdf <- brick("CONUS_TAVG_Gridded_0p25.nc",varname="climatology")
clim_ncdf <- as(clim_ncdf, "SpatRaster")

#Past_temperature
#get_right_years(1851-1959)_and_months(Jan-June)_and_create_stacks
sd3136 <- list()
 for (i in 1851:1859) {
   i <- as.character(i)
   stack3136_year  <- tmp_ncdf[i]
   stackyear35_abs <- sum(stack3136_year, clim_ncdf)
   stackyear3136_spring <- stackyear35_abs[[1:6]]
    sd3136[[i]] <-  stackyear3136_spring
     } 
spring3136SD <- rast(sd3136)
pastclim <- init(spring3136SD, "cell")
pastclim2 <- c(spring3136SD, pastclim)


#subset_past_records_and_convert_to_spatial_points
all_species_pastpresent_past <- subset(all_species_pastpresent, period=="past")
all_species_pastpresent_past2  <- distinct(dplyr::select(all_species_pastpresent_past, latitude, longitude,year))
all_species_pastpresent_past2_point <- st_as_sf(x = all_species_pastpresent_past2, 
                                            coords = c("longitude", "latitude"),
                                            crs = crs(pastclim))

#extract_cell_id
all_species_pastpresent_past2_point$lyr <- terra::extract(pastclim, vect(all_species_pastpresent_past2_point))
all_species_pastpresent_past2_point2 <- 
all_species_pastpresent_past2_point %>%
  mutate(longitude = st_coordinates(.)[,1],
         latitude = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  as.data.frame()
all_species_pastpresent_past2_point2$lyr1 <- all_species_pastpresent_past2_point2$lyr$lyr1
all_species_pastpresent_past2_point3 <- all_species_pastpresent_past2_point2 %>% dplyr::select(-lyr)
all_species_pastpresent_past3 <- left_join(all_species_pastpresent_past, all_species_pastpresent_past2_point3,by=c("latitude","longitude","year"))

#convert_clime_stack_to_dataframe_and_pivot_longer
pastclim2_df<- as.data.frame(pastclim2, xy=TRUE)
pastclim2_df2 <- pastclim2_df %>% 
  pivot_longer(cols = -c(x,y,lyr1), names_to = c("year", "month"), 
               names_pattern = "(\\d+)_(.*)")
pastclim_df3 <- pastclim2_df2 %>% 
  pivot_wider(names_from = month, values_from = value)
pastclim_df4 <- as.data.frame(pastclim_df3)
pastclim_df4$year <- as.numeric(pastclim_df4$year)

#make_sure_years_in_numeric_and_join_jan_june_temperature_to_dataframe
all_species_pastpresent_past3$year <- as.numeric(all_species_pastpresent_past3$year)
all_species_pastpresent_past_done <- left_join(all_species_pastpresent_past3, pastclim_df4, by=c("lyr1","year"))

################################################################### present_day_assembly_for_NPN
#present_npn
#processing_PRISM_data
prism_set_dl_dir("~/prismtmp")
sd_all <- get_prism_monthlys(type = "tmean", year = 2010:2023, mon=1:6, keepZip = FALSE)
sd_stack <- pd_stack(prism_archive_subset("tmean", "monthly", years = 2010:2023, mon = 1:6))
sd_stack2 <- as(sd_stack, "SpatRaster")

#reduce_area_closer_to_EasternUSA
e <-extent(-100, -60, 30, 48)
sd_stack3 <-crop(sd_stack2, e)

#get_PRISM_stack_with_data_for_right_years_and_months
st <- ym("2010-01")
en <- ym("2023-12")
dates <- st %m+% months(seq(0, round(interval(st, en) / months(1)), 1))
dates2 <- as.data.frame(dates) %>%
  mutate(monthyr = format_ISO8601(dates, precision = "ym")) 
dates2$monthyr2 <-  ym(dates2$monthyr) 
dates3 <- as.data.frame(dates2) %>%
  mutate(month = month(monthyr2)) %>%
  filter(month %in% c(1:6))
names(sd_stack2) <- dates3$monthyr
e <-extent(-100, -60, 30, 48)
sd_stack3 <-crop(sd_stack2, e)

#associate_cell_ids
presnpn_clim <- init(sd_stack3, "cell")
presnpn_clim2 <- c(sd_stack3,presnpn_clim)

#process_data_records_subset_to_presentNPN_and_convert_to_spatial_points_and_get_cell_id
all_species_pastpresent_presNPN <- subset(all_species_pastpresent, period=="present" & inat=="0")
all_species_pastpresent_presNPN2  <- distinct(dplyr::select(all_species_pastpresent_presNPN, latitude, longitude,year))
all_species_pastpresent_presNPN2_point <- st_as_sf(x = all_species_pastpresent_presNPN2, 
                                                           coords = c("longitude", "latitude"),
                                                           crs = crs(presnpn_clim2))
all_species_pastpresent_presNPN2_point$lyr <- terra::extract(presnpn_clim, vect(all_species_pastpresent_presNPN2_point))
all_species_pastpresent_presNPN2_point2 <- 
  all_species_pastpresent_presNPN2_point %>%
  mutate(longitude = st_coordinates(.)[,1],
         latitude = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  as.data.frame()
all_species_pastpresent_presNPN2_point2$lyr1 <- all_species_pastpresent_presNPN2_point2$lyr$lyr1
all_species_pastpresent_presNPN2_point3 <- all_species_pastpresent_presNPN2_point2 %>% dplyr::select(-lyr)
all_species_pastpresent_presNPN3 <- left_join(all_species_pastpresent_presNPN, all_species_pastpresent_presNPN2_point3, by=c("latitude","longitude","year"))
all_species_pastpresent_presNPN3$year <- as.character(all_species_pastpresent_presNPN3$year)

#convert_stack_to_data_frame_pivot_longer
presnpn_clim2_df<- as.data.frame(presnpn_clim2, xy=TRUE)
presnpn_clim2_df2 <- presnpn_clim2_df %>% 
  pivot_longer(cols = -c(x,y,lyr1), names_to = c("year", "month"), 
               names_pattern = "(\\d+)-(.*)")
presnpn_clim_df3 <- presnpn_clim2_df2 %>% 
  pivot_wider(names_from = month, values_from = value)
presnpn_clim_df4 <- as.data.frame(presnpn_clim_df3)

#join_temperature_data_to_NPNPresent_and_bind_together_past_and_presentNPN
all_species_pastpresent_presNPN_done <- left_join(all_species_pastpresent_presNPN3, presnpn_clim_df4, by=c("lyr1","year"))
names(all_species_pastpresent_presNPN_done)[15:20] <- sub("0", "", names(all_species_pastpresent_presNPN_done)[15:20])
all_species_pastpresent_n_outlier3 <- rbind(all_species_pastpresent_past_done, all_species_pastpresent_presNPN_done)

####################################################################### iNaturalist_temp_assembly
#iNaturalist_present
all_species_pastpresent_presiNat <- subset(all_species_pastpresent, period=="present" & inat=="1")

#get_region
test <- rnaturalearth::ne_states(country = c("United States of America"),returnclass = "sf") %>%
  filter(name %in% c("Pennsylvania", "Delaware", "Maryland", "District of Columbia", "Michigan", "Ohio", "Indiana", "Iowa", "Illinois", "Wisconsin", "Minnesota", "West Virginia", "Virginia", "North Carolina", "Tennessee", "Kentucky", "South Carolina", "Georgia", "Alabama", "Missouri", "Mississippi", "Florida"))

test2<-st_transform(test, sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
test3 <-st_union(test2)

#make_grid_over_region
grids <- st_make_grid(test3, cellsize = c(55000, 55000))
grids2 <- st_join(st_as_sf(grids), st_as_sf(test3))
grids2 = dplyr::mutate(st_sf(geometry = grids2), id_cells = 1:n())
grids3 <- grids2[test3, col = '#ff000088']

#aggregate_stack_to_grid
sd_stack4 <- aggregate(sd_stack3, fact=2, fun=mean)
sdnpn2 <- crop(sd_stack4, test3, mask= T)
sd_annual_bg <- extract(sdnpn2, grids3, fun="mean")
grids5sd <- dplyr::mutate(grids3,sd_annual_bg)

#housekeeping_to_get_latlon
grids4 <- st_transform(grids3,crs("EPSG:4326"))
grid_cent <- st_centroid(grids4)
grid_wcent <- dplyr::mutate(grids3,st_coordinates(grid_cent))
grid_wcent2 <-st_drop_geometry(grid_wcent)
test <- as.matrix(grid_wcent2)
cell_latitude_long <- cbind(test[,1],test[,2],test[,3])
cell_latlon2 <- as.data.frame(cell_latitude_long)
cell_latlon3 <- cell_latlon2 %>% dplyr::rename(id_cells=V1, longitude=V2, latitude=V3)

#housekeeping_on_latlon_joins
cell_latlon3$latitude <- round(cell_latlon3$latitude, digits=5)
cell_latlon3$longitude <- round(cell_latlon3$longitude, digits=5)
all_species_pastpresent_presiNat$latitude <- round(all_species_pastpresent_presiNat$latitude , digits=5)
all_species_pastpresent_presiNat$longitude <-  round(all_species_pastpresent_presiNat$longitude , digits=5)

#join_iNat_data_get_grid_cells
cell_latlon3$latitude2 <- as.character(cell_latlon3$latitude)
cell_latlon3$longitude2 <- as.character(cell_latlon3$longitude)
all_species_pastpresent_presiNat$latitude2 <- as.character(all_species_pastpresent_presiNat$latitude)
all_species_pastpresent_presiNat$longitude2 <- as.character(all_species_pastpresent_presiNat$longitude)
all_species_pastpresent2_presinat_idc <- left_join(all_species_pastpresent_presiNat,cell_latlon3, by=c("latitude","longitude"))

#do_filtering_remove_unneeded_gridcells
cell_filter <- unique(all_species_pastpresent2_presinat_idc$id_cells)
#grids6tempm <- grids5mean %>% dplyr::filter(id_cells %in% cell_filter)
grids6tempsd <- grids5sd %>% dplyr::filter(id_cells %in% cell_filter)
grids6tempsd2 <- grids6tempsd %>% st_drop_geometry()

#pivot_to_form
grids6tempsd2_2 <- grids6tempsd2 %>% 
  pivot_longer(cols = -c(id_cells,ID), names_to = c("year", "month"), 
               names_pattern = "(\\d+)-(.*)")
grids6tempsd2_3 <- grids6tempsd2_2 %>% 
  pivot_wider(names_from = month, values_from = value)
grids6tempsd2_4 <- as.data.frame(grids6tempsd2_3)

#join_everything_up
all_species_pastpresent2_presinat_idc$year <- as.character(all_species_pastpresent2_presinat_idc$year)
all_species_pastpresent_presiNat_done <- left_join(all_species_pastpresent2_presinat_idc, grids6tempsd2_4, by=c("id_cells","year"))
names(all_species_pastpresent_presiNat_done)[18:23] <- sub("0", "", names(all_species_pastpresent_presiNat_done)[18:23])

#all_together
all_species_pastpresent_n_outlier4 <- dplyr::select(all_species_pastpresent_n_outlier3, -lyr1,-x,-y)
all_species_pastpresent_presiNat_done2 <- dplyr::select(all_species_pastpresent_presiNat_done, -latitude2.x,-longitude2.x,-id_cells, -latitude2.y,-longitude2.y,-ID)     
all_species_pastpresent_n_outlier_d <- rbind(all_species_pastpresent_n_outlier4,all_species_pastpresent_presiNat_done2)
#write_results
write.csv(all_species_pastpresent_n_outlier_d, "all_species_pastpresent_n_outlier_d.csv" )

#############################################################################3
#temperature_and_slopes
all_species_pastpresent_n_slope <-  dplyr::select(all_species_pastpresent_n_outlier_d,-doy,-interval, -all_species_dist,-inat)
all_species_pastpresent_n_slope_longer <- all_species_pastpresent_n_slope  %>% 
  pivot_longer(cols = -c(...1,latitude,longitude,year,taxon,period,state), names_to = "month", values_to="temp")

all_species_pastpresent_n_slope_longer <- all_species_pastpresent_n_slope_longer %>% drop_na(temp)
all_species_pastpresent_n_slope_longer$month <- as.numeric(all_species_pastpresent_n_slope_longer$month)

#prepare_data_and_get_temp_sd_and_slope_for_Feb_June
all_species_pastpresent_n_slope2 <- subset(all_species_pastpresent_n_slope_longer, month!= "1")
#all_species_pastpresent_n_slope2_longer <- all_species_pastpresent_n_slope2  %>% 
#  pivot_longer(cols = -c(latitude,longitude,year,taxon,period,state), names_to = "month", values_to="temp")
all_species_pastpresent_n_slope2 <- all_species_pastpresent_n_slope2 %>% drop_na(temp)
all_species_pastpresent_n_slope2$month <- as.numeric(all_species_pastpresent_n_slope2$month)
temp_sd_out <- all_species_pastpresent_n_slope2 %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  summarize(springtemp.x = mean(temp), springsd2=sd(temp))
slope_out2 <- all_species_pastpresent_n_slope2 %>%
group_by(latitude,longitude,year,taxon,period,state) %>%
  do(model = lm(temp ~ month, data = .)) %>%
  mutate(slope2 = model$coefficients[2])

#join_data_back_to_main_data_frame
slope_out3 <- left_join(temp_sd_out,slope_out2, by=c("latitude","longitude","year","taxon","period","state"))

slope_out32<- slope_out3 %>% dplyr::select(-model)
         
all_species_done_2 <- left_join(all_species_pastpresent_n_outlier_d,slope_out32, by=c("latitude", "longitude","year","taxon","period","state"))

##################################################################3
#now_once_more_splitting_data_by_intervals
#create_lat_breaks
all_species_done_3 <- all_species_done_2 %>%
  mutate(interval = cut(`latitude`, 
                         breaks = c(-Inf, 36, 41, Inf),
                         labels = c("0", "1","2")))

#subset_by_breaks
all_species_done_3_3136 <- subset(all_species_done_3, interval=="0")
all_species_done_3_3641 <- subset(all_species_done_3, interval=="1")
all_species_done_3_4146 <- subset(all_species_done_3, interval=="2")

#subset_temp_values_for_each_break
#Jan_to_April_for_3136
#Feb_to_May_for_3641
#March_to_June_for_4146
all_species_done_3_3136_2 <- all_species_done_3_3136[,c(1:15,18:20)]
all_species_done_3_3641_2 <- all_species_done_3_3641[,c(1:11,13:16,18:20)]
all_species_done_3_4146_2 <- all_species_done_3_4146[,c(1:11,14:20)]

####################################################3
#run_temp_and_slopes_separately_for_each_interval_example_for_4146
all_species_done_3_4146_2 <-  dplyr::select(all_species_done_3_4146_2,-doy,-interval, -all_species_dist,-inat,-springtemp.x,-springsd2,-slope2)
all_species_done_3_4146_longer <-   all_species_done_3_4146_2 %>%
  pivot_longer(cols = -c(...1,latitude,longitude,year,taxon,period,state), names_to = "month", values_to="temp")
all_species_done_3_4146_longer <- all_species_done_3_4146_longer %>% drop_na(temp)
all_species_done_3_4146_longer$month <- as.numeric(all_species_done_3_4146_longer$month)

temp_sd_out_4146 <- all_species_done_3_4146_longer %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  summarize(springtemp.y = mean(temp), springsd1=sd(temp))

slope_out_4146 <-  all_species_done_3_4146_longer %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  do(model = lm(temp ~ month, data = .)) %>%
  mutate(slope1 = model$coefficients[2])

slope_out_4146_3 <- left_join(temp_sd_out_4146,slope_out_4146, by=c("latitude","longitude","year","taxon","period","state"))
slope_out_4146_32<- slope_out_4146_3 %>% dplyr::select(-model)

all_species_pastpresent_n_outlier_d3_4146 <- left_join(all_species_done_3_4146,slope_out_4146_32, by=c("latitude", "longitude","year","taxon","period","state"))

#################################################################
#band31-36
all_species_done_3_3136_2 <-  dplyr::select(all_species_done_3_3136_2,-doy,-interval, -all_species_dist,-inat,-springtemp.x,-springsd2,-slope2)
all_species_done_3_3136_longer <-   all_species_done_3_3136_2 %>%
  pivot_longer(cols = -c(...1,latitude,longitude,year,taxon,period,state), names_to = "month", values_to="temp")
all_species_done_3_3136_longer <- all_species_done_3_3136_longer %>% drop_na(temp)
all_species_done_3_3136_longer$month <- as.numeric(all_species_done_3_3136_longer$month)

temp_sd_out_3136 <- all_species_done_3_3136_longer %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  summarize(springtemp.y = mean(temp), springsd1=sd(temp))

slope_out_3136 <-  all_species_done_3_3136_longer %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  do(model = lm(temp ~ month, data = .)) %>%
  mutate(slope1 = model$coefficients[2])

slope_out_3136_3 <- left_join(temp_sd_out_3136,slope_out_3136, by=c("latitude","longitude","year","taxon","period","state"))
slope_out_3136_32<- slope_out_3136_3 %>% dplyr::select(-model)

all_species_pastpresent_n_outlier_d3_3136 <- left_join(all_species_done_3_3136,slope_out_3136_32, by=c("latitude", "longitude","year","taxon","period","state"))
######################################################################
#band36-41
all_species_done_3_3641_2 <-  dplyr::select(all_species_done_3_3641_2,-doy,-interval, -all_species_dist,-inat,-springtemp.x,-springsd2,-slope2)
all_species_done_3_3641_longer <-   all_species_done_3_3641_2 %>%
  pivot_longer(cols = -c(...1,latitude,longitude,year,taxon,period,state), names_to = "month", values_to="temp")
all_species_done_3_3641_longer <- all_species_done_3_3641_longer %>% drop_na(temp)
all_species_done_3_3641_longer$month <- as.numeric(all_species_done_3_3641_longer$month)

temp_sd_out_3641 <- all_species_done_3_3641_longer %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  summarize(springtemp.y = mean(temp), springsd1=sd(temp))

slope_out_3641 <-  all_species_done_3_3641_longer %>%
  group_by(latitude,longitude,year,taxon,period,state) %>%
  do(model = lm(temp ~ month, data = .)) %>%
  mutate(slope1 = model$coefficients[2])

slope_out_3641_3 <- left_join(temp_sd_out_3641,slope_out_3641, by=c("latitude","longitude","year","taxon","period","state"))
slope_out_3641_32<- slope_out_3641_3 %>% dplyr::select(-model)

all_species_pastpresent_n_outlier_d3_3641 <- left_join(all_species_done_3_3641,slope_out_3641_32, by=c("latitude", "longitude","year","taxon","period","state"))
#########################################################################################3
#final_Data_set
all_species_pastpresent_n_outlier_6 <- rbind(all_species_pastpresent_n_outlier_d3_3136, all_species_pastpresent_n_outlier_d3_3641, all_species_pastpresent_n_outlier_d3_4146)
write.csv(all_species_pastpresent_n_outlier_6, file="all_species_pastpresent_n_outlier_62.csv")
                                                              