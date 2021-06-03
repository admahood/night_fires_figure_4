# extract GOES16 data to fired polygons

libs<- c("tidyverse", "sf", "lubridate", "data.table", "vroom",
         "doParallel", "foreach")
invisible(sapply(libs, library, character.only=TRUE, quietly=TRUE))

# metadata ===========================================
modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

brazil_dates<- c(as.Date("2019-07-28"),as.Date("2019-08-18"))

# data import ==================================================================

# choose appropriate goes files for dl

raw_goes_aws <- system(paste("aws s3 ls", 
                             "s3://earthlab-mkoontz/goes16",
                             "--recursive"), 
                   intern = TRUE) 

goes_files <- raw_goes_aws %>%
  as_tibble() %>%
  tidyr::separate(value, into=c("s3date", "s3time", "size", "filename"), sep = "\\s+") %>%
  dplyr::select(filename) %>%
  tidyr::separate(filename, sep = "_", into = c("x1", "x2", "x3", "x4", "x5",
                                         "starttime", "endtime", "creationtime"),
                  remove= FALSE) %>%
  dplyr::select(starttime, filename) %>%
  dplyr::mutate(date = as.Date(str_sub(starttime,2,8), "%Y%j"),
                datetime = ymd_hms(paste(date,
                                     str_sub(starttime, 9,14)))) %>%
  dplyr::filter(date > as.Date("2002-01-01"))


# getting brazil files

brazil_files <- goes_files %>%
  dplyr::filter(date > brazil_dates[1] & date < brazil_dates[2]) %>%
  pull(filename)

# the slowest possible way
for(i in brazil_files){
  system(paste("aws s3 cp", 
               file.path("s3://earthlab-mkoontz", i),
               file.path("data/goes16", "brazil",i),
               "--only-show-errors"))
}

braz <- st_read("data/brazil_fire.gpkg")%>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))%>%
  mutate(acq_date = as.Date(acq_date)) %>%
  filter(acq_date > as.Date("2019-07-28"),
         acq_date < as.Date("2019-08-18")) %>%
  st_transform(crs = modis_crs) %>%
  st_coordinates %>%
  as_tibble() %>%
  summarise(xmin = min(X),
            ymin = min(Y),
            xmax = max(X),
            ymax = max(Y))

braz_goes<-list.files("data/goes16/brazil/goes16",full.names=TRUE, pattern = ".csv") %>%
  vroom() 

braz_goes_counts <- braz_goes %>%
  filter(sinu_x>=braz$xmin[1], 
         sinu_x<=braz$xmax[1],
         sinu_y>=braz$ymin[1],
         sinu_y<=braz$ymax[1]) %>%
  group_by(rounded_datetime) %>%
  summarise(n = n()) %>%
  ungroup()

# ggplot(braz_goes %>%
#          filter(sinu_x>=braz$xmin[1], 
#                 sinu_x<=braz$xmax[1],
#                 sinu_y>=braz$ymin[1],
#                 sinu_y<=braz$ymax[1]) %>%
#          st_as_sf(coords = c("sinu_x", "sinu_y"), crs = modis_crs)%>%
#          st_buffer(dist = 1500, endCapStyle = "SQUARE")) +
#   geom_sf(alpha = 0.25, fill="blue", color = "transparent")

all_times <- seq(braz_goes_counts$rounded_datetime %>% min,
                 braz_goes_counts$rounded_datetime %>% max, by=3600)

tibble(rounded_datetime = all_times) %>%
  left_join(braz_goes_counts)%>%
  replace_na(list(n=0)) %>%
  write_csv("data/brazil_goes_counts.csv")
