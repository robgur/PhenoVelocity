library(lme4)
library(lmerTest)
library(see)
library(ggeffects)
library(ggplot2)
library(sjPlot)
library(performance)
library(patchwork)
library(dplyr)
library(viridis)
library(rnaturalearth)
library(sf)
library(ncf)


#####################################################################
##################### Set global ggplot2 theme ######################
#####################################################################

library(ggplot2); theme_set(theme_classic(base_size = 16, base_family = "arial", base_line_size = 0))


## read in data

all_species_pastpresent_no8 <- read.csv("your_working_directory/all_species_pastpresent_no8.csv", stringsAsFactors = F)



#####################################################################
########################## Figure 1 #################################
#####################################################################


### for map figure ###

library(dplyr)

#get_region

test <- rnaturalearth::ne_states(country = c("United States of America"),returnclass = "sf")

sites<-c("New York", "Maine", "New Hampshire", "New Jersey", "Connecticut","Rhode Island", "Vermont","Pennsylvania", "Delaware", "Maryland", "District of Columbia", "Michigan", "Ohio", "Indiana", "Iowa", "Illinois", "Wisconsin", "Minnesota", "West Virginia", "Virginia", "North Carolina", "Tennessee", "Kentucky", "South Carolina", "Georgia", "Alabama", "Missouri", "Mississippi", "Florida", "Massachusetts", "Louisiana")

class(test$name)

test <- test %>% filter(name %in% sites)

test2<-st_transform(test, sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

test3 <-st_union(test2)

#make_grid_over_region

grids <- st_make_grid(test3, cellsize = c(55000, 55000))
grids2 <- st_join(st_as_sf(grids), st_as_sf(test3))
grids2 = dplyr::mutate(st_sf(geometry = grids2), id_cells = 1:n())
grids3 <- grids2[test3, col = '#ff000088']

plot(grids3)

st_write(grids3, "grids3.shp")

## intersect sample points with grid cells. This was done in QGIS v3.16 using the "join attributes by location" function, which sums the number of points within each polygon grid cell. The attribute table was then exported to a .csv file named samp_eff_id_cells_periods.csv and merged to the grids3 object in R for plotting.

pnts<-read.csv("samp_eff_id_cells_periods.csv",stringsAsFactors = F)

out<-merge(grids3,pnts,by.x = "id_cells", by.y = "id_cells")


pas<-out %>% filter(period == "past")
p<-out %>% filter(period == "present")


# build maps for past and present data sets of sampling of study area with latitude lines designated at seasonal windows

# Latitudinal intervals: 31° to 36°, 36° and 41° and 41° to 46° 


p1 <- ggplot() +
  geom_sf(data = test2, fill = "gray97", color = "gray20") +
  geom_sf(data = pas, aes(fill = log(sampling_effort))) +
  scale_fill_gradientn(colours = rev(inferno(6)),
                       limits=c(0,5)) +
  scale_y_continuous(breaks = c(25,31,36,41,46,50)) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "aliceblue"),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        panel.grid = element_line(color = "white", size = 0.8))


p2 <- ggplot() +
  geom_sf(data = test2, fill = "gray97", color = "gray20") +
  geom_sf(data = p, aes(fill = log(sampling_effort))) +
  scale_fill_gradientn(colours = rev(inferno(6)),
                       limits=c(0,5)) +
  scale_y_continuous(breaks = c(25,31,36,41,46,50)) +
  theme_bw() + 
  theme(#legend.position = "none",
    panel.background = element_rect(fill = "aliceblue"),
    #axis.ticks = element_blank(),
    #axis.text = element_blank(),
    panel.grid = element_line(color = "white", size = 0.8))

p1 + p2

## export as pdf and adjust formatting in Adobe Illustrator

#####################################################################
########################## Figure 3 #################################
#####################################################################

## final best LMM

model_at4 <- lmer(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1 + springtemp2_sc + latitude_sc | taxon) + (1 | year) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



p<-plot_model(model_at4, type = "pred", terms = c("latitude_sc","period","state","timing[early, mid,late]"),colors = c("cadetblue", "gold","orangered"),  title = "")

p

## create model objects with a consistent limits on the y-axis

p1<-p[[1]] + ylim(40,180)
p2<-p[[2]] + ylim(40,180)
p3<-p[[3]] + ylim(40,180)

## use patchwork to make a multipanel plot pulling them back together
p1+p2+p3

## export as pdf and adjust formatting in Adobe Illustrator

#####################################################################
########################## Figure 4 #################################
#####################################################################


p<-plot_model(model_at4, type = "pred", terms = c("springtemp2_sc","slope2_sc","period","timing[early, mid,late]"),colors = c("cadetblue", "gold","orangered"),  title = "")


## create model objects with a consistent limits on the y-axis

p1<-p[[1]] + ylim(60,160)
p2<-p[[2]] + ylim(60,160)
p3<-p[[3]] + ylim(60,160)

## using patchwork make multipanel plot pulling them back together
p1+p2+p3

## export as pdf and adjust formatting in Adobe Illustrator

#####################################################################
################### Supplementary Figure 1 ##########################
#####################################################################

data<-all_species_pastpresent_no8 

f<-data %>% group_by(year,state) %>%
  summarize(n = n())

f$year<-as.factor(year)

ggplot(f, aes(fill=state, y=n, x=as.factor(year))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## export as pdf and adjust formatting in Adobe Illustrator

#####################################################################
################### Supplementary Figure 2 ##########################
#####################################################################

## Run best LMM and extract residuals from model

model_at4 <- lmer(doy ~ springtemp2_sc + slope2_sc + period + state + timing + latitude_sc + (1 + springtemp2_sc + latitude_sc | taxon) + (1 | year) + springtemp2_sc:slope2_sc + springtemp2_sc:period + springtemp2_sc:latitude_sc + slope2_sc:timing + period:state + period:timing + period:latitude_sc + state:latitude_sc + timing:latitude_sc, data=all_species_pastpresent_no8, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


## check for residual spatial autocorrelation

resids<-residuals(model_at4)
out<-cbind(all_species_pastpresent_no8,resids)


fit1 <- spline.correlog(x = out2$longitude, y = out$latitude, z = out$resids, resamp = 10, na.rm = T)
summary(fit1)


plot(fit1, ylim = c(-0.2,0.2), xlim = c(0,8))

## export as pdf and adjust formatting in Adobe Illustrator


#####################################################################
################### Supplementary Figure 3 ##########################
#####################################################################

## plot best LMM residuals with latitude and longitude

plot(out2$longitude,out2$resids, cex.axis = 1.25)
plot(out2$latitude,out2$resids, cex.axis = 1.25)

## export as pdf and adjust formatting in Adobe Illustrator
