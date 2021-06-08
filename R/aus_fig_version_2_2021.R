# Night fires figure

# setup ========================================================================
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(ggsn)
library(stars)
library(ggpubr)
library(scales)
library(lubridate)
library(cowplot)
library(ggmap)
theme_set(theme_void())

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# thresholds ====================================================================
thresholds<-read_csv("data/updated-goes-af-vpd-thresholds.csv")
# thresholds[1:20,]

# landcover classification =====================================================
classification_table <- data.frame(
  name = c("Evergreen Needleleaf Forests", "Evergreen Broadleaf Forests",
           "Deciduous Needleleaf Forests","Deciduous Broadleaf Forests",
           "Mixed Forests","Closed Shrublands",
           "Open Shrublands","Woody Savannas",
           "Savannas","Grasslands",
           "Permanent Wetlands","Croplands",
           "Urban and Built-up Lands","Cropland/Natural  Vegetation  Mosaics",
           "Permanent Snow and Ice","Barren",
           "Water Bodies", "Unclassified"),
  value = c(1:17,255)
)


lut_lc <-  c("1"="Evergreen Needleleaf Forests", "2"= "Evergreen Broadleaf Forests",
             "3"= "Deciduous Needleleaf Forests","4"= "Deciduous Broadleaf Forests",
             "5"= "Mixed Forests","6"= "Closed Shrublands",
             "7"= "Open Shrublands","8"= "Woody Savannas",
             "9"=  "Savannas","10"= "Grasslands",
             "11"= "Permanent Wetlands","12"= "Croplands",
             "13"=  "Urban and Built-up Lands",
             "14"= "Cropland/Natural  Vegetation  Mosaics",
             "15"=  "Permanent Snow and Ice","16"= "Barren",
             "17"=  "Water Bodies", "255"= "Unclassified")

lut_colors <- c("Evergreen Needleleaf Forests" = "#006400",
                "Evergreen Broadleaf Forests" = "#228B22",
                "Deciduous Needleleaf Forests"= "#458B00",
                "Deciduous Broadleaf Forests" = "#008B45", 
                "Mixed Forests"= "#3CB371",
                "Closed Shrublands" = "#6E8B3D",
                "Open Shrublands" = "#9ACD32", 
                "Woody Savannas" = "#6B8E23",
                "Savannas" = "#8B8B00",
                "Grasslands" = "#CDC673", 
                "Permanent Wetlands" = "#00868B", 
                "Croplands" = "#EE9572",
                "Urban and Built-up Lands" = "#B3B3B3", 
                "Cropland/Natural  Vegetation  Mosaics" = "#EE8262",
                "Permanent Snow and Ice" = "#FFFFFF",
                "Barren" = "#DEB887",
                "Water Bodies" = "#87CEEB",
                "Unclassified" = "#BEBEBE"
                )

df_colors <- read_csv("data/landcover-colors.csv")

lut_colors <- df_colors %>%
  dplyr::select(color) %>%
  pull 
names(lut_colors) <- df_colors$lc_name
lut_colors[12] <- "#87CEEB"
names(lut_colors)[12] <- "Water Bodies"
lut_colors[13] <- "#BEBEBE"
names(lut_colors)[13] <- "Unclassified"
lut_colors[14] <- "#B3B3B3"
names(lut_colors)[14] <- "Urban and Built-up Lands"
lut_colors[15] <- "#DEB887"
names(lut_colors)[15] <- "Barren"

daynight_cols <- c("#2166AC","#B2182B") # red is #B2182B

snowy_lc_file <- "data/MCD12Q1_tifs/h29v12.tif"
braz_lc_file <- "data/MCD12Q1_tifs/h12v09.tif"
tubbs_lc_file <- "data/MCD12Q1_tifs/h08v05.tif"

modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

# snowy complex ================================================================

snowy <- st_read("data/snowy_complex.gpkg") %>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))%>%
  mutate(acq_date = as.Date(acq_date)) %>%
  filter(acq_date > as.Date("2019-12-28"),
         acq_date < as.Date("2020-01-6"))
daterange_s <- snowy$acq_date %>% as.Date() %>% range()
bbox <- st_bbox(snowy)%>% as.numeric()

snowy_modis <- read_csv("data/snowy_modis_afd.csv")%>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = modis_crs)%>%
  dplyr::mutate(daynight = factor(daynight, levels = c("night", "day")))

snowy_modis_detections <- snowy_modis %>%
    mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
    group_by(acq_datetime) %>%
    summarise(value = n(),
              daynight= first(daynight)) %>%
    ungroup() %>%
    mutate(variable = "MODIS Detections",
           lty = 0,
           label_y = max(value)*0.85) %>%
    dplyr::select(date = acq_datetime, variable, value,lty, label_y, daynight) %>%
    st_set_geometry(NULL) %>%
    filter(date>daterange_s[1], date < daterange_s[2])
  
snowy_n <- snowy_modis %>%
  st_set_geometry(NULL)%>%
  mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
  group_by(acq_datetime) %>%
  summarise(value = n(),
            daynight = first(daynight)) %>%
  ungroup() %>%
  mutate(variable = "Detections",
         lty = 0,
         lty_main = 0) %>%
  mutate(label_y = max(value)*0.85) %>%
  dplyr::select(datetime = acq_datetime, variable, value,lty, label_y, daynight, lty_main) 

snowy_clim <- read_csv("data/snowyave.csv",col_names = F) %>%
  dplyr::rename(day=X1, hour=X2, value = X3)%>% 
  mutate(year = ifelse(day>300,2019,2020),
         month = ifelse(day>300,"12","01"),
         value = value/10) %>%
  mutate(day = replace(day, day > 300, day-334),
         day = str_pad(day, 2, "left", "0"),
         date = paste0(year, "-",month,"-",day) %>% as.Date(),
         datetime = lubridate::ymd_hms(paste0(date," ", hour, ":00:00")),
         variable = "VPD (kPa)",
         lty=2,
         lty_main = 1) %>%
  mutate(label_y = max(value)*0.85) %>%
  filter(date > as.Date("2019-12-28"),
         date < as.Date("2020-01-6"))%>%
  dplyr::mutate(latitude = first(snowy$latitude),
                longitude = first(snowy$longitude),
                acq_datetime = datetime,
                acq_hour = hour(acq_datetime),
                acq_min = minute(acq_datetime),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  dplyr::select(datetime, variable, value,lty, label_y,daynight, lty_main) %>%
  rbind(snowy_n)

lc_s <- snowy_modis %>%
  mutate(lc=raster::extract(x=raster(snowy_lc_file), y=.)) %>%
  mutate(lc_name = lut_lc[lc]) %>%
  group_by(lc_name) %>%
  summarise(percent_lc = round((n()/nrow(.))*100)) %>%
  ungroup() %>%
  st_set_geometry(NULL)

stats_s <- snowy_modis %>%
  st_set_geometry(NULL) %>%
  group_by(daynight) %>%
  summarise(percent_dn = round((n()/nrow(.))*100,2)) %>%
  ungroup() 
  


snowy_area <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Australia") %>%
  st_crop(snowy %>% st_buffer(1))

australia <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Australia") 


if(!file.exists("data/fishnet.RDS")){
  fishnet <- st_bbox(snowy)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01)%>%
    as("Spatial")%>%
    raster::extract(x=raster(snowy_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/fishnet.RDS")
  }else{
  fishnet <- readRDS("data/fishnet.RDS")
}

if(!file.exists("data/fishnet_lc.RDS")){
  fishnet_lc <- st_bbox(snowy) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc, "data/fishnet_lc.RDS")
}else{
    fishnet_lc<-readRDS("data/fishnet_lc.RDS")
  }


locator_box <- st_bbox(snowy) %>% st_as_sfc()


ls_snow<-list.files("data/landsat_snowy/may7", full.names = TRUE) %>%
  as_tibble()%>%
  filter(str_detect(value,c("B[234]"))) %>%
  pull(value) %>%
  terra::rast() %>%
  terra::aggregate(10)

names(ls_snow)<- c("b", "g", "r")
snowy_area_rast<- st_crs(ls_snow)

# removing extreme values
# ls_snow[ls_snow<5000] <-NA
# ls_snow[ls_snow>20000] <-NA
library(MODIS)

filenames<-MODIS::getSds("data/modis_snowy/MOD09A1.A2020049.h29v12.006.2020058053841.hdf")

# rgb 143
# r=1, g=4, b=3

hdf<-terra::rast("data/modis_snowy/MOD09A1.A2020049.h29v12.006.2020058053841.hdf")

library(stars)

st_hdf<-st_as_stars(hdf)

hdf_t <- st_transform(st_hdf, st_crs(snowy_area))

plot(hdf_t)

names(hdf) <- filenames$SDSnames
snowy_area_rast<- st_crs(hdf)
r_bbox <- snowy%>%
  st_transform(crs=snowy_area_rast)%>%
  st_bbox()%>% as.numeric()


RStoolbox::ggRGB(as(hdf#%>%terra::stretch(0,255)
                    , "Raster"),
                 r=1,g=4,b=3,
                 stretch = "lin") +
  # geom_sf(data = st_transform(snowy_area, crs=snowy_area_rast),
  #         fill = "transparent") +
  geom_sf(data = st_transform(snowy_modis, crs = snowy_area_rast), 
          aes(color=daynight), show.legend = "point")  +
  scale_fill_manual(values = lut_colors,
                    name = "Land Cover")+
  scale_alpha_manual(values = c(0.5,0.15))+
  scale_color_manual(values = daynight_cols)+
  ggsn::scalebar(data = snowy, location = "topleft", model = "WGS84",
                 dist = 30, dist_unit = "km",transform = TRUE,st.dist = 0.05) +
  guides(fill = FALSE)+
  xlim(c(r_bbox[c(1,3)])) +
  ylim(c(r_bbox[c(2,4)])) +
  ggtitle(paste("A. Snowy Complex. December 29, 2019 - January 5, 2020")) +
  theme(legend.position ="none",
        # legend.justification = c(0.95,0),
        # legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA))
                          

ls_snow_rgb <- ls_snow %>%
  as.data.frame(xy=TRUE)%>%
  na.omit 



ggplot(ls_snow_rgb,aes(x=b)) +
  geom_histogram(stat="count")

testplot<-ggplot(ls_snow_rgb,aes(x=x,y=y,fill=dnbr))+
  geom_tile() +
  ggsave("testplot.png")

RStoolbox::ggRGB(as(ls_snow, "Raster"))

# snowy plots ==================================================================
main_plot_s <- ggplot() +
  geom_sf(data = snowy_area, fill = "transparent")+
  geom_sf(data = snowy_modis, 
          aes(color=daynight), show.legend = "point") +
  scale_fill_manual(values = lut_colors,
                    name = "Land Cover")+
  scale_alpha_manual(values = c(0.5,0.15))+
  scale_color_manual(values = daynight_cols)+
  ggsn::scalebar(data = snowy, location = "topleft", model = "WGS84",
                 dist = 30, dist_unit = "km",transform = TRUE,st.dist = 0.05) +
  guides(fill = FALSE)+
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggtitle(paste("A. Snowy Complex. December 29, 2019 - January 5, 2020")) +
  theme(legend.position ="none",
        # legend.justification = c(0.95,0),
        # legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA))

locator_plot_s <- ggplot(australia) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=st_centroid(locator_box), color = "red", size=2)+
  ylim(c(-43, -11))+
  xlim(c(114,153))

inset_s <- ggplot(snowy_clim %>% filter(variable == "VPD (kPa)"), aes(x=datetime, y=value)) +
  geom_line() +
  geom_bar(data = snowy_modis_detections,stat = "identity", width = 20000,
           aes(x=date, y=value, fill=daynight))+
  geom_abline(aes(slope = 0,intercept = 0.794), lty=2)+
  scale_linetype_manual(values = c(0,2))+
  scale_y_continuous(position = "right", labels = label_number_si())+
  scale_x_datetime(date_breaks = "1 week", date_labels = "%b %d")+
  facet_wrap(~variable, scales = "free_y", 
             nrow = 2#, strip.position = "left"
             ) +
  scale_fill_manual(values = daynight_cols)+
  xlab("Date") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

snowy_cow <- ggdraw() +
  draw_plot(main_plot_s) +
  draw_plot(locator_plot_s, 0.8,.68,.25,.25)

# some fires in brazil ==========================================================
# pull in global forest cover and/or mod17 to make sure it's not an ag fire

brazil_viirs <- st_read("data/brazil_fire.gpkg")%>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))%>%
  mutate(acq_date = as.Date(acq_date)) %>%
  filter(acq_date > as.Date("2019-07-28"),
         acq_date < as.Date("2019-08-18"))
daterange_b <- brazil_viirs$acq_date %>% as.Date() %>% range()

brazil_modis <- readr::read_csv("data/brazil_modis_afd.csv")%>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = modis_crs)%>%
  dplyr::mutate(daynight = factor(daynight, levels = c("night", "day")))

brazil_viirs_detections <- brazil_viirs%>%
  st_set_geometry(NULL) %>%
  mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
  group_by(acq_datetime) %>%
  summarise(value = n()) %>%
  ungroup() %>%
  mutate(variable = "VIIRS Detections",
         lty = 0,
         label_y = max(value)*0.85) %>%
  dplyr::select(datetime = acq_datetime, variable, value,lty, label_y) 

brazil_modis_detections <- brazil_modis %>%
  mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
  group_by(acq_datetime) %>%
  summarise(value = n(),
            daynight= first(daynight)) %>%
  ungroup() %>%
  mutate(variable = "MODIS Detections",
         lty = 0,
         label_y = max(value)*0.85) %>%
  dplyr::select(date = acq_datetime, variable, value,lty, label_y, daynight) %>%
  st_set_geometry(NULL) %>%
  filter(date>daterange_b[1], date < daterange_b[2])

brazil_goes  <- read_csv("data/brazil_goes_counts.csv") %>%
  dplyr::mutate(variable = "GOES Detections",
                lty = 0,
                label_y = max(n)*0.85)%>%
  dplyr::select(datetime = rounded_datetime, variable, value = n, lty, label_y)

braz_clim <- read_csv("data/brazil2019vpd.csv") %>%
  dplyr::rename(value = `vpd(hPa)`)%>% 
  mutate(value = value/10,
         year = 2019,
         month = ifelse(day>212,"08","07")) %>%
  mutate(day = ifelse(day > 212, day-212, day-181))%>%
  mutate(day = str_pad(day, 2, "left", "0"),
         date = paste0(year, "-",month,"-",day) %>% as.Date(),
         datetime = lubridate::ymd_hms(paste0(date," ", hour, ":00:00")),
         variable = "VPD (kPa)",
         lty = 2,
         label_y = max(value)*0.85) %>%
  filter(date > as.Date("2019-07-28"),
         date < as.Date("2019-08-18"))%>%
  dplyr::select(datetime, variable, value,lty, label_y) %>%
  rbind(brazil_goes)

lc_b <- brazil_modis %>%
  mutate(lc=raster::extract(x=raster(braz_lc_file), y=.)) %>%
  mutate(lc_name = lut_lc[lc]) %>%
  group_by(lc_name) %>%
  summarise(percent_lc = round((n()/nrow(.))*100)) %>%
  ungroup() %>%
  st_set_geometry(NULL)

stats_b <- brazil_modis %>%
  st_set_geometry(NULL) %>%
  group_by(daynight) %>%
  summarise(percent_dn = round((n()/nrow(.))*100,2)) %>%
  ungroup() 

braz_area <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Brazil") %>%
  st_crop(brazil_viirs %>% st_buffer(1))

brazil <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Brazil") 

bbox <- st_bbox(brazil_viirs)%>% as.numeric()

locator_box <- st_bbox(brazil_viirs) %>% st_as_sfc()
if(!file.exists("data/braz_fishnet.RDS")){
  fishnet <- st_bbox(braz)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.005)%>%
    as("Spatial")%>%
    raster::extract(x=raster(braz_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/braz_fishnet.RDS")
}else{
  fishnet <- readRDS("data/braz_fishnet.RDS")
}

if(!file.exists("data/braz_fishnet_lc.RDS")){
  fishnet_lc <- st_bbox(braz) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.005) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc, "data/braz_fishnet_lc.RDS")
}else{
  fishnet_lc<-readRDS("data/braz_fishnet_lc.RDS")
}
brazil_fired<- st_read("data/brazil_fired.gpkg")

# brazil plots =================================================================
main_plot_b <- ggplot() +
  geom_sf(data = fishnet_lc, 
          aes(fill=classes),
          color = "transparent", alpha = 0.5)+
  geom_sf(data = braz_area, fill = "transparent")+
  geom_sf(data = brazil_fired, fill = "transparent") +
  geom_sf(data = brazil_modis, 
          aes(color=daynight), 
          show.legend = "point") +
  scale_fill_manual(values = lut_colors,
                    name = "Land Cover")+
  scale_alpha_manual(values = c(0.75,0.5))+
  scale_color_manual(values = daynight_cols)+
  guides(fill = FALSE, color=FALSE)+
  ggsn::scalebar(data = brazil_viirs, location = "topleft", model = "WGS84",st.dist = 0.03,
                 dist = 5, dist_unit = "km",transform = TRUE) +
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggtitle(paste("C. Fires in Brazil. July 14 - September 30, 2019")) +
  theme(legend.position ="none",
        legend.justification = c(0,0),
        legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA))


locator_plot_b <- ggplot(brazil) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=locator_box, fill="transparent", color = "red", lwd=2)

inset_b <- ggplot(braz_clim,
                  aes(x=datetime, y=value)) +
  geom_line(color = "grey40") +
  geom_bar(data = brazil_modis_detections,stat="identity",width=20000,
           aes(y = value, fill=daynight, x=date), 
           color = "transparent")+
  geom_abline(aes(slope = 0,intercept = 01.33, lty=as.factor(lty)))+
  scale_linetype_manual(values = c(0,2))+
  facet_wrap(~variable, scales = "free_y", nrow = 3) +
  scale_y_continuous(position = "right", labels = label_number_si())+
  scale_x_datetime(date_breaks = "1 week", date_labels = "%b %d")+
  theme_bw() +
  xlab("Date")+
  scale_fill_manual(values = daynight_cols)+
  theme(legend.position = "none",
        legend.justification = c(0,1),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()); inset_b

braz_cow <- ggdraw() +
  draw_plot(main_plot_b) +
  draw_plot(locator_plot_b, .72,.65,.25,.25) 

# tubbs ==========================================================

tubbs_viirs <- st_read("data/tubbs.gpkg") %>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))
daterange_t <- tubbs$acq_date %>% as.Date() %>% range()

tubbs_modis <- read_csv("data/tubbs_modis_afd.csv")%>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = modis_crs)%>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))

tubbs_modis_detections<- tubbs_modis %>%
  mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
  group_by(acq_datetime) %>%
  summarise(value = n(),
            daynight= first(daynight)) %>%
  ungroup() %>%
  mutate(name = "MODIS Detections",
         lty = 0,
         label_y = max(value)*0.85) %>%
  dplyr::select(datetime_utc = acq_datetime, name, value,lty, label_y, daynight) %>%
  st_set_geometry(NULL) %>%
  filter(datetime_utc>daterange_t[1], datetime_utc < daterange_t[2])

# tubbs_p <- st_read("~/data/fire/mtbs/mtbs_perimeter_data/") %>%
#   filter(Fire_Name == "TUBBS")

tubbs_fired_p <- st_read("data/tubbs_fired.gpkg")

lct <- tubbs %>%
  mutate(lc=raster::extract(x=raster(tubbs_lc_file), y=.)) %>%
  mutate(lc_name = lut_lc[lc]) %>%
  group_by(lc_name) %>%
  summarise(percent_lc = round((n()/nrow(.))*100)) %>%
  ungroup() %>%
  st_set_geometry(NULL)

statst <- tubbs_modis %>%
  st_set_geometry(NULL) %>%
  group_by(daynight) %>%
  summarise(percent_dn = round((n()/nrow(.))*100,2)) %>%
  ungroup() 


tubbs_area <- st_read("~/data/background/CUS/CUS.shp") %>%
  st_transform(crs=st_crs(tubbs)) %>%
  st_crop(tubbs %>% st_buffer(1))

usa <- st_read("~/data/background/CUS/CUS.shp") %>%
  st_transform(crs=st_crs(tubbs)) %>%
  filter(STUSPS == "CA")

bbox <- st_bbox(tubbs)%>% as.numeric()

locator_box <- st_bbox(tubbs_area) %>% st_as_sfc()

if(!file.exists("data/fishnett.RDS")){
  fishnet <- st_bbox(tubbs)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01)%>%
    as("Spatial")%>%
    raster::extract(x=raster(tubbs_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/fishnett.RDS")
}else{
  fishnet <- readRDS("data/fishnett.RDS")
}

if(!file.exists("data/fishnet_lct.RDS")){
  fishnet_lc <- st_bbox(tubbs) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc, "data/fishnet_lct.RDS")
}else{
  fishnet_lc<-readRDS("data/fishnet_lct.RDS")
}

# tubbs plots ==================================================================
main_plot_t <- ggplot() +
  geom_sf(data = fishnet_lc, 
          aes(fill=classes),
          color = "transparent", alpha = 0.5)+
  geom_sf(data = tubbs_area, fill = "transparent")+
  # geom_sf(data = tubbs_p, fill = "transparent",lwd=1) +
  geom_sf(data = tubbs_fired_p, fill = "transparent",lwd=1) +
  geom_sf(data = tubbs_modis, 
          aes(color=daynight, alpha = daynight), size=4,show.legend = "point") +
  scale_alpha_manual(values = c(0.9,0.5))+
  scale_color_manual(values = daynight_cols)+
  scale_fill_manual(values = lut_colors)+
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggsn::scalebar(data = tubbs, location = "topleft", model = "WGS84",st.dist = 0.03,
                 dist = 2, dist_unit = "km",transform = TRUE) +
  guides(fill = FALSE) +
  ggtitle("B. Tubbs Fire. Oct 9-15, 2017") +
  theme(legend.position = "none",
        legend.justification = c(1,0),
        legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA))

inset_t <- read_csv("data/tubbs-hourly.csv") %>% 
  mutate(lty = ifelse(name == "VPD (kPa)",2,0),
         name = str_replace_all(name, "Active Fire detections", "Detections")) %>%
  group_by(name) %>%
  mutate(label_y = max(value)*0.85,
         daynight = NA,
         lty_main = 1) %>%
  ungroup()%>%
  dplyr::select(-fire_name, -id)%>%
  rbind(tubbs_modis_detections %>% mutate(lty_main = 0)) %>%
  filter(datetime_utc > as.Date("2017-10-08") &
           datetime_utc < as.Date("2017-10-17"))%>%
  ggplot(aes(x=datetime_utc, y=value)) +
  geom_line(color = "grey40", aes(lty = as.factor(lty_main))) +
  geom_bar(data = tubbs_modis_detections,
           stat="identity",width = 10000,
           aes(y = value, fill=daynight), 
           color = "transparent")+
    xlab("Date") +
    scale_x_datetime(date_breaks = "1 week", date_labels = "%b %d")+
    geom_hline(aes(yintercept = 0.794, lty=as.factor(lty)))+
    scale_y_continuous(position = "right", 
                       labels = label_number_si())+
    scale_linetype_manual(values=c(0,1,2)) +
    scale_fill_manual(values =daynight_cols)+
    facet_wrap(~name, ncol = 1, scales = "free_y", 
               nrow = 3, strip.position = "left") +
    theme_bw() +
    theme(legend.position = "none",
          # strip.background = element_blank(),
          # strip.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())

locator_plot_t <- ggplot(usa) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=st_centroid(locator_box), color = "red", size=2)

tubbs_cow<-ggdraw() +
  draw_plot(main_plot_t) +
  draw_plot(locator_plot_t, .75,.10,.2,.2)

# all together =================================================================
leg <- get_legend(ggplot()+
                    geom_bar(data = snowy_n,stat = "identity", width = 20000,
                             aes(x=datetime, y=value, fill=daynight))+
                    scale_fill_manual(values=daynight_cols,
                                      name= "Modis Active\nFire Detections")
                  )
insets_ls <- ggarrange(leg, inset_s,
                       nrow=2, ncol=1,heights = c(1,2.2),
                       labels = c("","F. Snowy"),
                       label.x = 0.15,
                       label.y = 0.95)

insets_tb <- ggarrange(inset_t, inset_b, 
                    nrow=1, ncol=2,
                    labels = c( "D. Tubbs", "E. Brazil"),
                    label.x = c(.85, .39), 
                    label.y = .95, hjust="right", vjust="top")

finalfig <- ggdraw(xlim = c(0,7.5), ylim = c(0, 12.75)) +
  draw_plot(snowy_cow, x = 0, y = 8.75, width = 7.5, height = 4) +
  # draw_plot(leg, x=6.14,y=8.1, width=2,height=2) +
  draw_plot(tubbs_cow, x = 0, y = 4, width = 3, height = 4.75) +
  draw_plot(braz_cow, x = 3, y=4, width = 4.5, height = 4.75) +
  draw_plot(insets_tb, x=0,y=0, width=5, height=4) +
  draw_plot(insets_ls, x=5, y=0, width=2.5, height=4)+
  # draw_label(label = "1 kPa", x = 2.18, y=0.7, size=12,fontface = "bold") +
  ggsave("images/panel_fig.png", height = 12.75, width = 7.5)
