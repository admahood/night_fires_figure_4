# Night fires figure

# setup ========================================================================
library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(raster)
library(ggsn)
library(stars)
library(ggpubr)
library(scales)
library(basemaps)
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
           "Urban and Built-up","Cropland/Natural  Vegetation  Mosaics",
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
             "13"=  "Urban and Built-up",
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
                "Urban and Built-up" = "grey10",
                "Cropland/Natural  Vegetation  Mosaics" = "#EE8262",
                "Permanent Snow and Ice" = "#FFFFFF",
                "Barren" = "#DEB887",
                "Water Bodies" = "#87CEEB",
                "Unclassified" = "#BEBEBE"
                )

lut_lc_simple <-  c("Evergreen Needleleaf Forests" = "Forests",
                    "Evergreen Broadleaf Forests" = "Forests",
                    "Deciduous Needleleaf Forests"= "Forests",
                    "Deciduous Broadleaf Forests" = "Forests",
                    "Mixed Forests"= "Forests",
                    "Closed Shrublands" = "Grasslands, Shrublands, Savannas",
                    "Open Shrublands" = "Grasslands, Shrublands, Savannas",
                    "Woody Savannas" = "Grasslands, Shrublands, Savannas",
                    "Savannas" = "Grasslands, Shrublands, Savannas",
                    "Grasslands" = "Grasslands, Shrublands, Savannas",
                    "Permanent Wetlands" = "Croplands, Other",
                    "Croplands" = "Croplands, Other",
                    "Urban and Built-up Lands" = "Urban and Built-up",
                    "Cropland/Natural  Vegetation  Mosaics" = "Croplands, Other",
                    "Permanent Snow and Ice" = "Croplands, Other",
                    "Barren" = "Croplands, Other",
                    "Water Bodies" = "Water Bodies",
                    "Unclassified" = "Croplands, Other"
)

lut_colors_simple <- c("Forests" = "#228B22",
                "Grasslands, Shrublands, Savannas" = "#CDC673",
                "Croplands, Other" = "grey",
                "Urban and Built-up" = "grey10",
                "Water Bodies" = "#87CEEB"
)

df_colors <- read_csv("data/landcover-colors-2021.csv") %>%
  mutate(lc_name= replace(lc_name, lc_name == "Urban and Built-up Lands",
                          "Urban and Built-up"))

lut_colors <- df_colors %>%
  dplyr::select(color) %>%
  pull
names(lut_colors) <- df_colors$lc_name
lut_colors[17] <- "#87CEEB"
names(lut_colors)[17] <- "Water Bodies"
lut_colors[18] <- "#BEBEBE"
names(lut_colors)[18] <- "Unclassified"




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

# crsss<- basemaps::basemap_stars(snowy_area) %>% st_crs

bbox_s <- snowy %>%
  # st_transform(crs = crsss)%>%
  # st_buffer(dist=10000) %>%
  st_bbox()%>%
  as.numeric()


australia <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Australia")


if(!file.exists("data/fishnet.RDS")){
  fishnet <- st_bbox(snowy)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01)%>%
    as("Spatial")%>%
    exactextractr::exact_extract(x=raster(snowy_lc_file), y=., 
                   fun = "mode") 
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

# snowy plots ==================================================================
# basemap_gglayer(snowy_area)->pp

main_plot_s <- ggplot() +
  geom_sf(data = fishnet_lc, 
          aes(fill=lut_lc_simple[classes]),
          alpha=0.5,
          color = "transparent")+
  scale_fill_manual(values = lut_colors_simple, name = "Landcover Classes")+
  geom_sf(data = snowy_modis,# %>% st_transform(crs=crsss), 
          aes(color=daynight), show.legend = "point") +
  scale_color_manual(values = daynight_cols)+
  # ggsn::scalebar(data = snowy_modis,#%>% st_transform(crs=crsss),
  #                location = "topleft", model = "WGS84",
  #                dist = 30, dist_unit = "km",transform = TRUE,st.dist = 0.05) +
  xlim(c(bbox_s[c(1,3)])) +
  ylim(c(bbox_s[c(2,4)])) +
  ggtitle(paste("a. Snowy Complex. December 29, 2019 - January 5, 2020")) +
  theme(legend.position ="none",
        plot.title = element_text(face = "bold"),
        panel.spacing = unit(0, "cm"),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(color="black", fill=NA))

locator_plot_s <- ggplot(australia) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=st_centroid(locator_box), color = "red", size=2)+
  ylim(c(-43, -11))+
  xlim(c(114,153))

inset_s <- ggplot(snowy_clim %>% filter(variable == "VPD (kPa)"), aes(x=datetime, y=value)) +
  geom_line(color = "grey30") +
  geom_bar(data = snowy_modis_detections,stat = "identity", width = 20000,
           aes(x=date, y=value, fill=daynight))+
  geom_abline(aes(slope = 0,intercept = 0.794), lty=2)+
  scale_linetype_manual(values = c(0,2))+
  scale_y_continuous(position = "right", labels = label_number_si())+
  scale_x_datetime(date_breaks = "1 week", date_labels = "%b %d")+
  facet_wrap(~variable, scales = "free_y", 
             nrow = 2, strip.position = "left"
             ) +
  scale_fill_manual(values = daynight_cols)+
  xlab("Date") +
  theme_bw() +
  theme(legend.position = "none",
        # strip.background = element_blank(),
        # strip.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

snowy_cow <- ggdraw() +
  draw_plot(main_plot_s) +
  draw_plot(locator_plot_s, x=0.823,y=.62,width=.2,height=.2)
#
# tubbs ==========================================================

tubbs_viirs <- st_read("data/tubbs.gpkg") %>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))
daterange_t <- tubbs_viirs$acq_date %>% as.Date() %>% range()

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

tubbs_fired_p <- st_read("data/tubbs_fired.gpkg")

lct <- tubbs_viirs %>%
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
  st_transform(crs=st_crs(tubbs_viirs)) %>%
  st_crop(tubbs_viirs %>% st_buffer(1))

usa <- st_read("~/data/background/CUS/CUS.shp") %>%
  st_transform(crs=st_crs(tubbs_viirs)) %>%
  filter(STUSPS == "CA")

bbox <- st_bbox(tubbs_viirs)%>% as.numeric()

locator_box <- st_bbox(tubbs_area) %>% st_as_sfc()

if(!file.exists("data/fishnett.RDS")){
  fishnet_t <- st_bbox(tubbs_viirs)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01)%>%
    as("Spatial")%>%
    raster::extract(x=raster(tubbs_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/fishnett.RDS")
}else{
  fishnet_t <- readRDS("data/fishnett.RDS")
}

if(!file.exists("data/fishnet_lct.RDS")){
  fishnet_lc_t <- st_bbox(tubbs_viirs) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc_t, "data/fishnet_lct.RDS")
}else{
  fishnet_lc_t<-readRDS("data/fishnet_lct.RDS")
}

# tubbs plots ==================================================================
main_plot_t <- ggplot() +
  geom_sf(data = fishnet_lc_t, 
          aes(fill=lut_lc_simple[classes]),
          alpha = 0.5, 
          color = "transparent")+
  geom_sf(data = tubbs_area, fill = "transparent")+
  geom_sf(data = tubbs_fired_p, fill = "transparent",lwd=1) +
  geom_sf(data = tubbs_modis, 
          aes(color=daynight), size=6,show.legend = "point") +
  scale_alpha_manual(values = c(0.9,0.5))+
  scale_color_manual(values = daynight_cols)+
  scale_fill_manual(values = lut_colors_simple)+
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggsn::scalebar(data = tubbs_viirs, location = "topleft", model = "WGS84",st.dist = 0.03,
                 dist = 2, dist_unit = "km",transform = TRUE) +
  guides(fill = FALSE) +
  ggtitle("b. Tubbs Fire. Oct 9-15, 2017") +
  theme(legend.position = "none",
        plot.title = element_text(face="bold"),
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
    geom_line(color = "grey30", aes(lty = as.factor(lty_main))) +
    geom_bar(data = tubbs_modis_detections,
             stat="identity",width = 25000,
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


# full plot without brazil =====================================================
leg_dn <- get_legend(ggplot()+
                    geom_bar(data = snowy_n,stat = "identity", width = 20000,
                             aes(x=datetime, y=value, fill=daynight))+
                    scale_fill_manual(values=daynight_cols,
                                      name= "Modis Active\nFire Detections"))

leg_lc <- get_legend(fishnet_lc %>%
                       st_set_geometry(NULL) %>%
                       group_by(classes) %>%
                       summarise(x = rnorm(1)) %>%
                       ungroup() %>%
                       na.omit() %>%
                       mutate(simple_classes = lut_lc_simple[classes]) %>%
                       ggplot() +
                        geom_raster(aes(fill=simple_classes, x=x,y=x), alpha=0.5) +
                     scale_fill_manual(values = lut_colors_simple)+
                     guides(fill=guide_legend(ncol=1, title = "Landcover Classes")))

insets_ls <- ggarrange(leg_dn, inset_s,
                       nrow=2, ncol=1,heights = c(1,4),
                       # labels = c("","d. Snowy"),
                       label.x = 0.08,
                       label.y = 0.95)

insets_ts <- ggarrange(inset_t, insets_ls, 
                       nrow=1, ncol=2,
                       labels = c( "c. Tubbs", ""),
                       label.x = c(.85, .39), 
                       label.y = .97, hjust="right", vjust="top")

finalfig <- ggdraw(xlim = c(0, 10), ylim = c(0, 9.25)) +
  draw_plot(snowy_cow, x = 0, y = 4.7, width =7.75, height = 4.75) +
  draw_plot(inset_s, x = 7.75, y= 4.79, width = 2.25, height = 3.75) +
  
  draw_plot(tubbs_cow, x = 0, y = 0.265, width = 3, height = 4.75) +
  draw_plot(inset_t, x = 3, y = 0, width = 2.25, height = 4.85)+
  draw_plot(leg_dn, x=5.4, y=2)+
  draw_plot(leg_lc,x=5.5, y=3, width = 2, height = 2) +
  ggsave("images/panel_fig_no_braz.png", height = 9.25, width =10)
