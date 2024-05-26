library(lme4)
library(lmerTest)
#load_final_dataset
setwd("C:/Users/robgu/Downloads")
all_species_pastpresent_n_outlier6 <- read_csv("C:/Users/robgu/OneDrive/Desktop/onset_pastpresent/all_species_pastpresent_n_outlier_d62.csv")
all_species_pastpresent_no2 <- subset(all_species_pastpresent_n_outlier6, doy>13 & doy<215)
all_species_pastpresent_no2 <- all_species_pastpresent_no2 %>% drop_na() 
all_species_pastpresent_no2$latitude_sc <- scale(all_species_pastpresent_no2$latitude)
all_species_pastpresent_no2$springsd1_sc <- scale(all_species_pastpresent_no2$springsd1)
all_species_pastpresent_no2$springsd2_sc <- scale(all_species_pastpresent_no2$springsd2)
all_species_pastpresent_no2$springtemp1_sc <- scale(all_species_pastpresent_no2$springtemp.x)
all_species_pastpresent_no2$springtemp2_sc <- scale(all_species_pastpresent_no2$springtemp.y)
all_species_pastpresent_no2$disttc_sc <- scale(all_species_pastpresent_no2$all_species_dist)
all_species_pastpresent_no2$slope2_sc <- scale(all_species_pastpresent_no2$slope2)
all_species_pastpresent_no2$slope1_sc <- scale(all_species_pastpresent_no2$slope1)
all_species_pastpresent_no2$longitude_sc <- scale(all_species_pastpresent_no2$longitude)
all_species_pastpresent_no8$year2 <- as.factor(all_species_pastpresent_no8$year)


#subset_by_lat_band
all_species_pastpresent_no23641 <- subset(all_species_pastpresent_no2, latitude>36 & latitude<41)
all_species_pastpresent_no23136 <- subset(all_species_pastpresent_no2, latitude>31 & latitude<36)
all_species_pastpresent_no24146 <- subset(all_species_pastpresent_no2, latitude>41 & latitude<46)


#Clean_cases_where_sampling_is_too_low_per_year
all_species_pastpresent_no2_flow <- subset(all_species_pastpresent_no2, state=="flowering")
all_species_pastpresent_no2_fol <- subset(all_species_pastpresent_no2, state=="foliage")
all_species_pastpresent_no3_flow <- all_species_pastpresent_no2_flow %>% group_by(taxon,year,state) %>% filter(n()>4)
all_species_pastpresent_no3_fol <- all_species_pastpresent_no2_fol %>% group_by(taxon,year,state) %>% filter(n()>4)
all_species_pastpresent_no3 <- rbind(all_species_pastpresent_no3_flow, all_species_pastpresent_no3_fol)
all_species_pastpresent_no4 <- all_species_pastpresent_no3 %>% group_by(taxon,year) %>% filter(sum(state =="flowering")!= 0)
all_species_pastpresent_no5 <- all_species_pastpresent_no4 %>% group_by(taxon,year) %>% filter(sum(state =="foliage")!= 0)

#Create_50m_buffer_for_past_phenometrics_and_filter_present_falling_outside
all_species_pastpresent_no5present<- subset(all_species_pastpresent_no5, period=="present")
all_species_pastpresent_no5past<- subset(all_species_pastpresent_no5, period=="past")
all_species_pastpresent_no5present2 <- all_species_pastpresent_no5present   %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
all_species_pastpresent_no5past2  <-all_species_pastpresent_no5past   %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
all_species_pastpresent_no5present2_sf <- st_transform(all_species_pastpresent_no5present2 ,crs =st_crs(grids))
all_species_pastpresent_no5past2_sf <- st_transform(all_species_pastpresent_no5past2 ,crs =st_crs(grids))
all_species_pastpresent_no5past2_sf_buffer <- st_buffer(all_species_pastpresent_no5past2_sf$geometry, 50000)
all_species_pastpresent_no5past2_sf_buffer2 <- st_union(all_species_pastpresent_no5past2_sf_buffer, by_feature = FALSE)
all_species_pastpresent_no5present2_sf2 <- all_species_pastpresent_no5present2_sf[lengths(st_intersects(all_species_pastpresent_no5present2_sf,all_species_pastpresent_no5past2_sf_buffer2)) == 1,]
all_species_pastpresent_no5present2_sf3 <- st_transform(all_species_pastpresent_no5present2_sf2 ,crs ="EPSG:4326")
all_species_pastpresent_no5present2_sf4 <- tidyr::extract(all_species_pastpresent_no5present2_sf3, geometry, into = c('longitude', 'latitude'), '\\((.*),(.*)\\)', conv = T)
all_species_pastpresent_no6 <- rbind(all_species_pastpresent_no5past, all_species_pastpresent_no5present2_sf4)

#add_in_early_mid_late_timing_variable
timing2 <- read_csv("C:/Users/robgu/Downloads/timing7.csv")
all_species_pastpresent_no8 <- merge(all_species_pastpresent_no6,timing2, by=c("taxon","state"))

#import_tree_for_18species_and_process
tree18 <- "tree.tre"
tree18imp <- ape::read.tree(tree18)
tree18imp$tip.label <- str_replace(tree18imp$tip.label, "^\\w{1}", tolower)
tree18imp$tip.label <- gsub(" ","_",tree18imp$tip.label)
all_species_pastpresent_no8$sp <- all_species_pastpresent_no8$taxon

#modeling_lmmm
model_at4 <- lmer(doy ~ (springtemp2_sc  + slope2_sc  + period + state + timing + latitude_sc)^2  + longitude_sc + (1  + springtemp2_sc + latitude_sc | taxon) + (1 | year2), data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#backwards stepwise
step(model_at4)
#best model
model_at4 <- lmer(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1 + springtemp2_sc + latitude_sc | taxon) + (1 | year) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:timing + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc + state:timing + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#check VIFS
vif(model_at4)
#removed one interaction and re-run
model_at4 <- lmer(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1 + springtemp2_sc + latitude_sc | taxon) + (1 | year) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:timing + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc  + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#R2
r2_nakagawa(model_at4)
library(lme4)
library(lmerTest)
#load_final_dataset
setwd("C:/Users/robgu/Downloads")
all_species_pastpresent_n_outlier6 <- read_csv("C:/Users/robgu/OneDrive/Desktop/onset_pastpresent/all_species_pastpresent_n_outlier_d62.csv")
all_species_pastpresent_no2 <- subset(all_species_pastpresent_n_outlier6, doy>13 & doy<215)
all_species_pastpresent_no2 <- all_species_pastpresent_no2 %>% drop_na() 
all_species_pastpresent_no2$latitude_sc <- scale(all_species_pastpresent_no2$latitude)
all_species_pastpresent_no2$springsd1_sc <- scale(all_species_pastpresent_no2$springsd1)
all_species_pastpresent_no2$springsd2_sc <- scale(all_species_pastpresent_no2$springsd2)
all_species_pastpresent_no2$springtemp1_sc <- scale(all_species_pastpresent_no2$springtemp.x)
all_species_pastpresent_no2$springtemp2_sc <- scale(all_species_pastpresent_no2$springtemp.y)
all_species_pastpresent_no2$disttc_sc <- scale(all_species_pastpresent_no2$all_species_dist)
all_species_pastpresent_no2$slope2_sc <- scale(all_species_pastpresent_no2$slope2)
all_species_pastpresent_no2$slope1_sc <- scale(all_species_pastpresent_no2$slope1)
all_species_pastpresent_no2$longitude_sc <- scale(all_species_pastpresent_no2$longitude)
all_species_pastpresent_no8$year2 <- as.factor(all_species_pastpresent_no8$year)


#subset_by_lat_band
all_species_pastpresent_no23641 <- subset(all_species_pastpresent_no2, latitude>36 & latitude<41)
all_species_pastpresent_no23136 <- subset(all_species_pastpresent_no2, latitude>31 & latitude<36)
all_species_pastpresent_no24146 <- subset(all_species_pastpresent_no2, latitude>41 & latitude<46)


#Clean_cases_where_sampling_is_too_low_per_year
all_species_pastpresent_no2_flow <- subset(all_species_pastpresent_no2, state=="flowering")
all_species_pastpresent_no2_fol <- subset(all_species_pastpresent_no2, state=="foliage")
all_species_pastpresent_no3_flow <- all_species_pastpresent_no2_flow %>% group_by(taxon,year,state) %>% filter(n()>4)
all_species_pastpresent_no3_fol <- all_species_pastpresent_no2_fol %>% group_by(taxon,year,state) %>% filter(n()>4)
all_species_pastpresent_no3 <- rbind(all_species_pastpresent_no3_flow, all_species_pastpresent_no3_fol)
all_species_pastpresent_no4 <- all_species_pastpresent_no3 %>% group_by(taxon,year) %>% filter(sum(state =="flowering")!= 0)
all_species_pastpresent_no5 <- all_species_pastpresent_no4 %>% group_by(taxon,year) %>% filter(sum(state =="foliage")!= 0)

#Create_50m_buffer_for_past_phenometrics_and_filter_present_falling_outside
all_species_pastpresent_no5present<- subset(all_species_pastpresent_no5, period=="present")
all_species_pastpresent_no5past<- subset(all_species_pastpresent_no5, period=="past")
all_species_pastpresent_no5present2 <- all_species_pastpresent_no5present   %>%
   st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
all_species_pastpresent_no5past2  <-all_species_pastpresent_no5past   %>%
   st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
all_species_pastpresent_no5present2_sf <- st_transform(all_species_pastpresent_no5present2 ,crs =st_crs(grids))
all_species_pastpresent_no5past2_sf <- st_transform(all_species_pastpresent_no5past2 ,crs =st_crs(grids))
all_species_pastpresent_no5past2_sf_buffer <- st_buffer(all_species_pastpresent_no5past2_sf$geometry, 50000)
all_species_pastpresent_no5past2_sf_buffer2 <- st_union(all_species_pastpresent_no5past2_sf_buffer, by_feature = FALSE)
all_species_pastpresent_no5present2_sf2 <- all_species_pastpresent_no5present2_sf[lengths(st_intersects(all_species_pastpresent_no5present2_sf,all_species_pastpresent_no5past2_sf_buffer2)) == 1,]
all_species_pastpresent_no5present2_sf3 <- st_transform(all_species_pastpresent_no5present2_sf2 ,crs ="EPSG:4326")
all_species_pastpresent_no5present2_sf4 <- tidyr::extract(all_species_pastpresent_no5present2_sf3, geometry, into = c('longitude', 'latitude'), '\\((.*),(.*)\\)', conv = T)
all_species_pastpresent_no6 <- rbind(all_species_pastpresent_no5past, all_species_pastpresent_no5present2_sf4)

#add_in_early_mid_late_timing_variable
timing2 <- read_csv("C:/Users/robgu/Downloads/timing7.csv")
all_species_pastpresent_no8 <- merge(all_species_pastpresent_no6,timing2, by=c("taxon","state"))

#import_tree_for_18species_and_process
tree18 <- "tree.tre"
tree18imp <- ape::read.tree(tree18)
tree18imp$tip.label <- str_replace(tree18imp$tip.label, "^\\w{1}", tolower)
tree18imp$tip.label <- gsub(" ","_",tree18imp$tip.label)
all_species_pastpresent_no8$sp <- all_species_pastpresent_no8$taxon

#modeling_lmmm
model_at4 <- lmer(doy ~ (springtemp2_sc  + slope2_sc  + period + state + timing + latitude_sc)^2  + longitude_sc + (1  + springtemp2_sc + latitude_sc | taxon) + (1 | year2), data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#backwards stepwise
step(model_at4)
#best model
model_at4 <- lmer(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1 + springtemp2_sc + latitude_sc | taxon) + (1 | year) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:timing + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc + state:timing + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#check VIFS
vif(model_at4)
#removed one interaction and re-run
model_at4 <- lmer(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1 + springtemp2_sc + latitude_sc | taxon) + (1 | year) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:timing + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc  + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#R2
r2_nakagawa(model_at4)

#fit_pglmm_model
mod_at_phy <- phyr::pglmm(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1|sp__) + (springtemp2_sc | sp__) + (latitude_sc | sp__) + (1 | year2) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc  + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, cov_ranef = list(sp = tree18imp))

#check_model_in_DHARMa
simulationOutput <- simulateResiduals(fittedModel = model_at4, plot = F)
plot(simulationOutput)
