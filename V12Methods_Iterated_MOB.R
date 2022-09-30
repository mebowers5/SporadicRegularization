############################################################################################################
############## Create random tracks and test biases in the distribution summaries that result ##############
############## from crawl movement models and pathroutr re-routing in the FACT Network #####################
############################################################################################################
##############################################

#Remove everything from environment
rm(list=ls())

# Clear cache
gc()

# Increase R memory limit
memory.limit(size = 1.0000E+12)

##MOB: Don't need to worry about exporting the parallel package to the cluster
library(sf)
library(dplyr)

#Load required libraries
library(momentuHMM)
library(crawl)
library(lubridate)

#Create function to make line from each pair of coordinates
make_line <- function(start_x, start_y, end_x, end_y) {
  sf::st_linestring(matrix(c(start_x, end_x, start_y, end_y), 2, 2))
}

# Set working directory to file destination of iteration objects
# setwd("D:/Data/IterObjects")
setwd("/Volumes/Bowers_PhD/IterObjects/")

# Objects required for iterations to run
Bathy500m <- sf::st_read("Bathy500m.shp")

# Call in receiver stations
stations <- sf::st_read("stations.shp") 

### Check whether many retryFits are needed###
anims <- 1
yr <- 3
npts <- anims*yr*365

# Run test lines with one individual
test_lines <- Bathy500m %>% 
  sf::st_sample(geom, size=npts, type ="random", exact = TRUE) %>%  #Create random points inside 500m bathymetric isoline polygon
  as.data.frame %>%  
  mutate(AnimalID = as.factor(sample(1:anims, npts, replace = TRUE)),  #Randomly assign an Animal ID
         year = as.numeric(sample(1:yr, npts, replace = TRUE)),  #Randomly assign a year
         lat_direct = as.numeric(sample(1:2, npts, replace = TRUE)),   #Randomly assign a direction (north (1) or south (2))
         geom = geometry, 
         POINT_X = sf::st_coordinates(geom)[,1], POINT_Y = sf::st_coordinates(geom)[,2]) %>%  #Obtain x and y coordinates so that you can force northward and southward movement
  group_by(AnimalID, year) %>%
  arrange(AnimalID, year, lat_direct, ifelse(lat_direct == 2, desc(POINT_Y), POINT_Y), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(uid = 1:nrow(.)) %>%   #Assign sequential numbers to rows to preserve directionality
  dplyr::select(uid, AnimalID, year, geom, POINT_Y, POINT_X) %>%
  arrange(uid) %>%
  group_by(AnimalID) %>%
  mutate(t = row_number(),
         time = as.Date(t, origin = "2000-01-01"),   #Make up dates for the row numbers by Iteration and Animal ID
         start_x = POINT_X, start_y = POINT_Y, end_x = lead(POINT_X), end_y = lead(POINT_Y)) %>%  #Prep data for making lines - ensures proper order
  st_as_sf

#Make lines and unite geometries by uid's so that each line segment maintains time information
test_lines <- test_lines %>%
  filter(!is.na(end_y)) %>%
  tidyr::nest() %>% 
  mutate(
    data = purrr::map(data, 
                      ~ mutate(.x,
                               x = purrr::pmap(.l = list(start_x, start_y, end_x, end_y),
                                               .f = make_line))))

# Maintain detailed time info in lines to derive acoustic telemetry detection data later
test_lines <- test_lines %>%
  tidyr::unnest(data) %>%
  mutate(x = st_sfc(x)) %>% 
  st_as_sf(sf_column_name = 'x', crs = 26917)

### Find intersections of tracks and receiver stations
# These intersections represent data that a researcher might receive from an animal tagged with an acoustic transmitter that traveled along the track that we created

# Derive acoustic detection data through the intersection of lines and receiver ranges
# Cannot just take the intersection because the lines self intersect
# Take the symmetrical difference between the lines and self intersection to get the true derived detection data
test_lines <- sf::st_intersection(test_lines, stations$geometry) %>%  # x provides intersection geometry - linestrings that cross station polygons
  st_centroid() # Take the centroid of the linestrings

test_lines <- st_snap(test_lines, st_centroid(stations$geometry), tolerance = 650) # Snap the linestring centroids to the station point geometry - geometry 'x' is altered

### Run "derived acoustic telemetry data through continuous time correlated random walk model
#Load data----

#In real data, problem caused by the zero movement between locations at the same acoustic receiver, close together in time - solved by averaging daily locations
#Average daily locations
test_lines <- test_lines %>%
  arrange(AnimalID, time) %>% 
  rename(ID = AnimalID) %>% 
  group_by(ID, time) %>% 
  summarize(mean_pos = st_combine(x)) %>% 
  st_centroid() %>% 
  mutate(time = as.POSIXct(time + lubridate::hours(12)), 
         x = sf::st_coordinates(mean_pos)[,1], 
         y = sf::st_coordinates(mean_pos)[,2],
         locType = "o") %>% # Obtain x and y coordinates so that you can merge with predicted locations later
  as.data.frame

# Create a blank data frame
all_aic <- data.frame(matrix(ncol = 4, nrow = 0))

# Rename columns
colsnames <- c("aic_0", "aic_1", "aic_3", "aic_5")

colnames(all_aic) <- colsnames

# Number of replicates desired to test retryFits paramater values 
run <- 1:10

# Temporary memory is causing an issue. Run inside another replicate function to "bank" information along the way.
# aics <- sapply(1:2 , FUN = function(I) {
base::replicate(10, {
  
  # Clear cache
  gc()
  
  # Run replicates of crawlWrap function with multiple retryFits parameters
  df <- sapply(run, FUN = function(I) {
    crw_aic_0 <- momentuHMM::crawlWrap(obsData = test_lines, timeStep = "15 min", fixPar = c(NA, NA), theta = c(0,0),
                                       attempts = 15, retryFits = 0)
    
    crw_aic_1 <- momentuHMM::crawlWrap(obsData = test_lines, timeStep = "15 min", fixPar = c(NA, NA), theta = c(0,0),
                                       attempts = 15, retryFits = 1)
    
    crw_aic_3 <- momentuHMM::crawlWrap(obsData = test_lines, timeStep = "15 min", fixPar = c(NA, NA), theta = c(0,0),
                                       attempts = 15, retryFits = 3)
    
    crw_aic_5 <- momentuHMM::crawlWrap(obsData = test_lines, timeStep = "15 min", fixPar = c(NA, NA), theta = c(0,0),
                                       attempts = 15, retryFits = 5)
    
    aic_0 = as.numeric(crw_aic_0[["crwFits"]][['1']][["aic"]]) 
    aic_1 = as.numeric(crw_aic_1[["crwFits"]][['1']][["aic"]]) 
    aic_3 = as.numeric(crw_aic_3[["crwFits"]][['1']][["aic"]])
    aic_5 = as.numeric(crw_aic_5[["crwFits"]][['1']][["aic"]])
    
    vals <- list(aic_0, aic_1, aic_3, aic_5)
    
    return(unlist(vals))
    
  }
  )
  
  # Create data frame from results
  df <- as.data.frame(t(df)) %>% # Produces column for each run that contains AIC values with rows for each retryFits condition - need to transpose rows and columns
    rename(aic_0 = V1, aic_1 = V2, aic_3 = V3, aic_5 = V4) %>%
    mutate(trials = seq_along(aic_0)) %>%
    dplyr::select(trials, aic_0, aic_1, aic_3, aic_5)
  
  # Bind columns to create factored aic values
  df <- cbind(df[1], stack(df[2:5])) %>%
    mutate(ind = as.factor(ind)) %>%
    as.data.frame
  
  # Read in working aics file
  all_aic <- read.csv(file = "aics.csv") %>% dplyr::select(-X)
  
  # Append results
  all_aic <- rbind(all_aic, df)
  
  write.csv(all_aic, file = "aics.csv")
}
)

aics <- read.csv(file = "aics.csv") %>% dplyr::select(-X) 

# Run a power analysis to make sure enough samples were run
# Need help? https://www.statmethods.net/stats/power.html
pwr::pwr.anova.test(k = 4, n = NULL, f = .1, sig.level = 0.05, power = 0.8) # k = np. of groups, n = lowest common sample size between groups, sig.level = the difference you want to be able to predict, power = 0.8 is standard but can leave NULL to calculate
# power = 0.09872555
# Can determine effect size of f = 0.3788014 with 80% power

# Compare averages from each AIC value to one another (ANOVA? followed by Tukey test?)
aic.aov <- anova(lm(values ~ ind, data = aics))
aic.aov

# Analysis of Variance Table
# 
# Response: values
# Df Sum Sq Mean Sq F value Pr(>F)
# ind        2   2227 1113.60   1.955 0.1509
# Residuals 57  32469  569.63  

### End retryFits testing


{   
#Remove everything from environment
rm(list=ls())

# Clear cache
gc()

# Increase R memory limit
memory.limit(size = 1.0000E+12)

# If desired, change the temporary directory for R
# write("TMP = 'C:\Program Files\R\R-4.1.2\rtemp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

library(dplyr)
library(sf)
library(pathroutr)

#########################################################
# Setup parallel
#########################################################

# Find out how many cores you have
library(parallel)

## MOB: Avoid using "logical cores". Apparently this can actually wind up
##    slowing things down.
numCores <- detectCores(logical = F)
numCores

# Create cluster
# BB: Should I just use all cores now that we are setting logical to FALSE? Are the logical cores the ones required to run background processing?
cl <- makeCluster(numCores-2, outfile = "cluster_out")

# Load libraries
clusterEvalQ(cl, {
  
  ##MOB: Don't need to worry about exporting the parallel package to the cluster
  library(sf)
  library(dplyr)
  
  ##MOB: Don't load/export Tidyverse -- there's a lot of overhead associated
  ##    with it. Best to just load/export the packages you need -
  # BB: I do need tidyverse for pathroutr - I'm not sure what all inside tidyverse is used in the background code for pathroutr
  library(momentuHMM)
  library(crawl)
  library(lubridate)
  library(pathroutr)
  library(sfnetworks)
  library(tidyverse)
  library(stplanr)
  
})

# Set working directory to file destination of iteration objects
# setwd("D:/Data/IterObjects")
setwd("/Volumes/Bowers_PhD/IterObjects/")

# Objects required for iterations to run
Bathy500m <- sf::st_read("Bathy500m.shp")

# we need to get our land polygon into a proper format; essentially, we want it to be
# a series of POLYGONs (and not MULTIPOLYGONs or GEOMETRYCOLLECTION).
land_barrier <- sf::st_read("atlcoast.shp") %>%
  sf::st_transform(26917) %>%
  #sf::st_crop(bb) %>%
  sf::st_collection_extract('POLYGON') %>%
  sf::st_cast('POLYGON')

# Call in receiever stations
stations <- sf::st_read("stations.shp") 

# Call in grid to analyze counts
sg50km <- sf::st_read("sg500m.shp") 

#Make a new table with all gid's to create a basis by which tables should be merged
gid <- seq(1:max(sg50km$gid))  #Create a range of every grid cell
all <- as.data.frame(gid)  #Create every combination of iteration and gid

# Create a spatial network graph ("visibility graph") outside of the iterations
# prt_visgraph will build our visual graph network from our land barrier object and 
# return a SpatialLinesNetwork / sfNetwork that has no edges that cross land
vis_graph <- pathroutr::prt_visgraph(land_barrier)

# Set working directory back to where other data may be pulled from and where items should be saved
setwd("D:/Data/")
# setwd("/Volumes/LaCie/Data")

# Number of iterations desired 
runs <- 1:10

# Adjustments in simulated data
anims <- 35     # Number of animals (35 for my data)
yr <- 3      # Number of years (3 for my data)
npts <- anims*yr*365  # Number of random points used to make tracks for total number of animals - gives you the potential for a daily position per animal per day

#Create function to make line from each pair of coordinates
make_line <- function(start_x, start_y, end_x, end_y) {
  sf::st_linestring(matrix(c(start_x, end_x, start_y, end_y), 2, 2))
}

# Mirror global environment objects on all cores
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))

################## Start loops ##########################

base::replicate(75, {
  all_difs <- parallel::parLapply(cl, runs, fun = function(I) {
    
### Populate random points to create tracks that are coerced to move north and south
    # If the variability itself is a problem, try set.seed
    # If that doesn't work try tryCatch(wrap everything that you want to "catch")
#Create random tracks
lines <- Bathy500m %>% 
  sf::st_sample(geom, size=npts, type ="random", exact = TRUE) %>%  #Create random points inside 500m bathymetric isoline polygon
  as.data.frame %>%  
  mutate(AnimalID = as.factor(sample(1:anims, npts, replace = TRUE)),  #Randomly assign an Animal ID
         year = as.numeric(sample(1:yr, npts, replace = TRUE)),  #Randomly assign a year
         lat_direct = as.numeric(sample(1:2, npts, replace = TRUE)),   #Randomly assign a direction (north (1) or south (2))
         geom = geometry, 
         POINT_X = sf::st_coordinates(geom)[,1], POINT_Y = sf::st_coordinates(geom)[,2]) %>%  #Obtain x and y coordinates so that you can force northward and southward movement
  group_by(AnimalID, year) %>%
  arrange(AnimalID, year, lat_direct, ifelse(lat_direct == 2, desc(POINT_Y), POINT_Y), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(uid = 1:nrow(.)) %>%   #Assign sequential numbers to rows to preserve directionality
  dplyr::select(uid, AnimalID, year, geom, POINT_Y, POINT_X) %>%
  arrange(uid) %>%
  group_by(AnimalID) %>%
  mutate(t = row_number(),
         time = as.Date(t, origin = "2000-01-01"),   #Make up dates for the row numbers by Iteration and Animal ID
         start_x = POINT_X, start_y = POINT_Y, end_x = lead(POINT_X), end_y = lead(POINT_Y)) %>%  #Prep data for making lines - ensures proper order
  st_as_sf

#Check that points are within the polygon - commented out for iterations
#plot(Bathy500m$geometry, col="white")
#plot(rand$geom, size = 0.5, col="blue", add=TRUE)

#Make lines and unite geometries by uid's so that each line segment maintains time information
lines <- lines %>%
  filter(!is.na(end_y)) %>%
  tidyr::nest() %>% 
  mutate(
    data = purrr::map(data, 
                      ~ mutate(.x,
                               x = purrr::pmap(.l = list(start_x, start_y, end_x, end_y),
                                               .f = make_line))))
  
#Unite geometries by AnimalID's so that you have complete tracks/lines summarized by ID for the grid cell count later
id_lines <- lines %>% 
  mutate(
    data = purrr::map(data,
                      ~ mutate(.x, x = st_sfc(x))),
    x = purrr::map(data, ~ st_union(st_set_geometry(.x, 'x'))), ##Preserves order
    x = purrr::map(x, ~ st_cast(.x, 'MULTILINESTRING')),
    ID = AnimalID  # BB added
  ) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(x) %>% 
  st_as_sf(sf_column_name = 'x', crs = 26917)

# Maintain detailed time info in lines to derive acoustic telemetry detection data later
lines <- lines %>%
  tidyr::unnest(data) %>%
  mutate(x = st_sfc(x)) %>% 
  st_as_sf(sf_column_name = 'x', crs = 26917)

### Find intersections of tracks and receiver stations
# These intersections represent data that a researcher might receive from an animal tagged with an acoustic transmitter that traveled along the track that we created

# Derive acoustic detection data through the intersection of lines and receiver ranges
# Cannot just take the intersection because the lines self intersect
# Take the symmetrical difference between the lines and self intersection to get the true derived detection data
lines <- sf::st_intersection(lines, stations$geometry) %>%  # x provides intersection geometry - linestrings that cross station polygons
 st_centroid() # Take the centroid of the linestrings
  
lines <-  st_snap(lines, st_centroid(stations$geometry), tolerance = 650) # Snap the linestring centroids to the station point geometry - geometry 'x' is altered

### Run "derived acoustic telemetry data through continuous time correlated random walk model
#Load data----

#In real data, problem caused by the zero movement between locations at the same acoustic receiver, close together in time - solved by averaging daily locations
#Average daily locations
lines <- lines %>%
  arrange(AnimalID, time) %>% 
  rename(ID = AnimalID) %>% 
  group_by(ID, time) %>% 
  summarize(mean_pos = st_combine(x)) %>% 
  st_centroid() %>% 
  mutate(time = as.POSIXct(time + lubridate::hours(12)), 
        x = sf::st_coordinates(mean_pos)[,1], 
        y = sf::st_coordinates(mean_pos)[,2],
        locType = "o") %>% # Obtain x and y coordinates so that you can merge with predicted locations later
  as.data.frame

# Clear cache
gc()

#Impute missing locations with continuous time movement model
crw_mod <- momentuHMM::crawlWrap(obsData = lines, timeStep = "15 min", fixPar = c(NA, NA), theta = c(0,0),
                      attempts = 5, retryFits = 0)
# plot(crw_mod)

# Clear cache
gc()

#Continue with movement model data products
#Get the data.frame of predicted locations
lines <- data.frame(crw_mod$crwPredict)
 
#Turn crwPredict product into data frame
lines <- lines %>%
  dplyr::select(ID, locType, time, mu.x, mu.y, speed) %>%
  mutate(locType = ifelse(is.na(locType), "p", locType)) %>%  #Specify predicted locs
  as.data.frame

cat("Reading in track data")

# Read in the track data; here, I'm filtering to just include the predicted locations
# and making sure the time is a proper POSIX. Neither of those two steps are
# necessary
lines <- lines %>%
  #dplyr::filter(locType == "p") %>% 
  dplyr::select(ID, time, mu.x, mu.y) %>%
  dplyr::mutate(time = lubridate::ymd_hms(time), ID = as.factor(ID))

cat("Converting the track data to CRS 26917")

# convert the track_data to sf and set the CRS; the bb step is just a way to limit
# the size of the land polygon and save some computation time when creating vis_graph
lines <- lines %>% sf::st_as_sf(coords = c("mu.x","mu.y"), crs = 26917)
#bb <- sf::st_as_sfc(sf::st_bbox(track_path))

cat("Nesting the tracks")
# there are multiple paths identified by ID; we'll group and nest for a proper
# tidyverse/list-column workflow
lines <- lines %>%
  group_by(ID) %>%
  tidyr::nest()

# Clear cache
gc()

cat("Trimming points from the land_barrier")

# the track cannot start or end within the land barrier; prt_trim() trims those out
lines <- lines %>%
  rowwise() %>%
  mutate(trim_data = list(pathroutr::prt_trim(data, land_barrier)))

cat("Creating re-routed points")

# here, we create our re-routed points; the return is a two column data frame with the
# index location in the original point data and the new geometry. The user can handle
# updating of those original point data or pass the result on to prt_update_points()
lines <- lines %>% dplyr::rowwise() %>%
  mutate(rrt_pts = list(prt_reroute(trim_data, land_barrier, vis_graph)))

# Clear cache
gc() 

cat("Updating paths with re-routed points")

# NOTE: previous versions of prt_update_points() had the argument order reversed from
# what it now requires. The updated geometry points are passed first (here, `rrt_pts`)
# and, then, the original data to be updated. This order should allow for easy piping
# from prt_reroute()
lines <- lines %>% dplyr::rowwise() %>%
  mutate(path_pts = list(prt_update_points(rrt_pts, trim_data)),
         path_lines = list(path_pts %>% summarise(do_union = FALSE) %>% sf::st_cast('LINESTRING')))  # do_union MUST be FALSE!

# Clear cache
gc() 

cat("Compiling all points and lines into single object")

# we need to rbind all of our lines and points into single objects that can be plotted
# Need to maintain "ID" column here - THIS WORKS!
lines$geom <- do.call(rbind, lines$path_lines)
lines$geom <- sf::st_set_crs(lines$geom, 26917)

lines <- lines %>% dplyr::select(ID, geom) 

cat("Convert rowwise df to sf object")

### Create grid that maintains row for each grid cell
#Convert from rowwise_df to sf object
lines <- lines %>%
  rowwise() %>%
  mutate(geom = sf::st_geometry(geom)) %>%
  dplyr::select(ID, geom) %>%
  ungroup() %>%
  st_as_sf(., sf_column_name = "geom")

cat("Counting animals per grid cell in derived data")

### Count distinct AnimalIDs in each grid cell
# Count IDs per grid cell for derived data
der_count <- sf::st_join(sg50km, lines, join = st_intersects)

der_count <- aggregate(ID ~ gid, data = der_count, FUN = length) %>%
  rename(dr_count = ID)

cat("Counting animals per grid cell in complete data")

# Count IDs per grid cell for complete data
comp_count <- sf::st_join(sg50km, id_lines, join = st_intersects)

comp_count <- aggregate(ID ~ gid, data = comp_count, FUN = length) %>%
  rename(cmp_count = ID) 

cat("Join count data to all grid cells file")

#Now that we have every combination of gid and Iteration in one variable, merge one of the count files to it - doesn't really matter which on but we'll use the complete data one here
# Note: Is the above comment wrong or is the code wrong? I think the comment is wrong - leftover from when iterations were assigned instead of looped.
tc <- left_join(all, comp_count, by = "gid") %>%
  mutate_all(~replace(., is.na(.), 0))  # Replace NA's with zeros so we can check count

# Check for code accuracy - Sum should equal # of observations in tc_comp
#sum(tc_comp$cmp_count > 0)

cat("Join derived count data to complete count data")

#Join derived count to data frame
tc <- left_join(tc, der_count, by = "gid") %>%
  mutate_all(~replace(., is.na(.), 0))  # Replace NA's with zeros so we can check count

cat("Find the difference in counts")

#Calculate differences in each grid cell by iteration and gid
tc <- tc %>%
  dplyr::select(gid, cmp_count, dr_count) %>%
  mutate(cmp_count = as.double(cmp_count),
  dr_count = as.double(dr_count),
  dif = cmp_count - dr_count) %>%
  dplyr::select(gid, cmp_count, dr_count, dif) %>%
  group_by(gid) %>%
  as.data.frame

return(tc$dif)

# Clear cache
gc()

}
)

################# After loops finish ##################
all_difs # Returns list return object from all runs per each node

#Turn list into data frame
difs2add_df <- t(sapply(all_difs, '[', 1:max(sapply(all_difs, length)))) %>%
  as.data.frame

difs_df <- read.csv(file = "difs_df.csv") %>%
  select(-X)

# Append results to others
difs_df <- rbind(difs_df, difs2add_df)

# Save the results
write.csv(difs_df, file = "difs_df.csv")

# Clear cache
gc()

})
}
# Stop the cluster
stopCluster(cl)

# Take average of differences by cell over iterations
#Calculate column means
difs_matrix <- read.csv(file = "Results_difs_df.csv") %>%
  select(-X) %>%
  as.matrix

avg <- apply(difs_matrix, 2, mean)

avg <- avg %>%
  as.data.frame %>%
  mutate(gid = seq_along(avg)) %>%
  rename(mean = ".")

# Take standard error of differences by cell over iterations
se <- apply(difs_matrix, 2, plotrix::std.error)

se <- se %>%
  as.data.frame %>%
  mutate(gid = seq_along(se)) %>%
  rename(se = ".")

# Merge averages and SE
final <- merge(avg, se, by = "gid") 

# Save averages and SE
write.csv(final, file = "avg330iters_se.csv")

#Merge spatial grid with averages and standard errors so it can be mapped
sg50km_final <- merge(sg50km, avg, all=TRUE)

sg50km_final <- merge(sg50km_final, se, all=TRUE) %>%
  dplyr::mutate(g_fit = ifelse(abs(mean) <= 1 & se <= 0.12608, 1, 0)) # Label good fit for all grid cells with an absolute value of the mean <= 1 and a SE <= 0.12608, which is the maximum SE for all excellent fits

head(sg50km_final)

# Save shapefile
sf::st_write(sg50km_final, "sg50km_final.shp")

# Plot results
# library(basemaps)
# get_maptypes()
# set_defaults(map_service = "esri", map_type = "world_ocean_base")
# bm <- basemap_gglayer(sg50km_final)
bm <- basemaps::basemap_ggplot(ext = st_bbox(sg50km_final), map_service = "esri", map_type = "world_ocean_base")

library("rnaturalearth")
library("rnaturalearthdata")

bounds <- st_bbox(sg50km_final)

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(., 26917)

world <- world$geometry %>%  
  st_make_valid() %>%
  st_crop(., xmin =430000, ymin = 2689000, xmax = 1502000, ymax = 46780000)
class(world)

sg50km_final <- st_transform(sg50km_final, 26917)
st_bbox(sg50km_final)

mean_plot <- ggplot() +
  basemapR::base_map(bounds, basemap = 'google-terrain') + 
  geom_sf(data = sg50km_final, aes(fill = mean)) + 
  scale_fill_viridis_c(option = "C") +
  theme_bw() +
  ggsn::north(sg50km_final, symbol = 12, location = "bottomright", scale = 0.15) +
  ggsn::scalebar(sg50km_final, location = "bottomright", dist = 100, dist_unit = "km", st.size = 3, transform = FALSE) # anchor = c(y = 50000, x = 72000),

# Possible basemap types in ggspatial::annotation_map_tile
# rosm::osm.types()
# [1] "osm"                    "opencycle"              "hotstyle"               "loviniahike"            "loviniacycle"           "hikebike"              
# [7] "hillshade"              "osmgrayscale"           "stamenbw"               "stamenwatercolor"       "osmtransport"           "thunderforestlandscape"
# [13] "thunderforestoutdoors"  "cartodark"              "cartolight" 

# library(OpenStreetMap)
# base <- openmap(upperLeft = c(bounds[1], bounds[2]), lowerRight = c(bounds[3], bounds[4]), type = "esri")
# base <- openmap(upperLeft = c(437896.5, 2689791.0), lowerRight = c(1502333.6, 4678301.8), type = "osm", minNumTiles = 9L)
# base <- openmap(upperLeft = c(42, -70), lowerRight = c(25, -81), type = "osm")
# base <- openmap(upperLeft = c(42, -70), lowerRight = c(25, -81), type = "osm")

mean_plot <- ggplot(sg50km_final) +
  # ggspatial::annotation_map_tile("hikebike", zoomin = -1) + 
  geom_sf(data = sg50km_final, aes(fill = mean)) + 
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme_classic()
  # ggsn::north(sg50km_final, symbol = 12, location = "bottomright", scale = 0.15) +
  # ggsn::scalebar(sg50km_final, location = "bottomright", dist = 100, dist_unit = "km", st.size = 3, transform = FALSE) # anchor = c(y = 50000, x = 72000),
mean_plot

dev.off()
jpeg(file = "methods_mean09012022.jpeg", units="in", width=5.38, height=3.90, res=300)

se_plot <- ggplot(sg50km_final) +
  geom_sf(aes(fill = se)) + 
  scale_fill_viridis_c(option = "D", direction = -1) +
  theme_classic()
se_plot

dev.off()
jpeg(file = "methods_se09012022.jpeg", units="in", width=5.38, height=3.90, res=300)

sg50km_df <- sg50km %>%
  as.data.frame() %>%
  dplyr::select(gid, p_a)

sg50km_plot <- merge(sg50km_df, sg50km_final) %>%
  mutate(g_fit = as.factor(g_fit),
         ex_fit = as.factor(ex_fit),
         plot_fit = as.factor(if_else(ex_fit == 1, 2, if_else(g_fit == 1, 1, 0)))) %>%
  st_as_sf

sf::st_write(sg50km_plot, "sg50km_plot.shp")

mapView(sg50km_plot, zcol = "plot_fit")


fill_col = c("yellow", "blue", "violet")
outline_col = c("grey", "black")

fit_plot <- ggplot(sg50km_plot) +
  # geom_sf(aes(fill = p_a), color = alpha("white", 0.5)) +
  geom_sf(aes(fill = plot_fit, color = p_a), stroke = 2) + 
  # scale_fill_viridis_d(option = "D", direction = -1) +
  theme_classic() +
  scale_fill_manual(values = fill_col) + 
  scale_color_manual(values = outline_col) +
  guides(fill = guide_legend(override.aes = list(color = fill_col)),
         color = guide_legend(override.aes = list(shape = 21)))

fit_plot

dev.off()
jpeg(file = "methods_fit09012022.jpeg", units="in", width=5.38, height=3.90, res=300)


# Count number of receivers in grid cells
# Tell whether receivers are in grid cells or not
rcvr_PA <- sf::st_intersects(stations, sg50km) %>%
  as.data.frame %>%
  rename(gid = col.id)

rcvr_PA <- merge(sg50km, rcvr_PA, by = "gid") %>%
  mutate(rcvr = "TRUE") %>%
  as.data.frame %>%
  dplyr::select(-geometry) 

rcvr_PA = rcvr_PA[!duplicated(rcvr_PA$gid),]

rcvr_PA <- merge(sg50km, rcvr_PA, by = "gid", all.x = TRUE) %>%
  mutate(rcvr = ifelse(is.na(.$rcvr), FALSE, TRUE))

sf::st_write(rcvr_PA, "rcvr_PA.shp")

### How many grid cells fit certain good fit definitions?
se <- function(a) sd(a) / sqrt(length(a))

wilcox <- function(x) {
  wilcox.test(x ~ 1, alternative = "two.sided", mu = 0)
}

setwd("/Volumes/Bowers_PhD/Data/")
df <- read.csv(file = "Results_difs_df.csv") %>%
  summarise(across(everything(), list(mean = mean, SE = se))) %>%
  tidyr::gather(., GID, values) %>%
  mutate(gid = gsub("_.*", "", GID), fx = gsub(".*_", "", GID)) %>%
  filter(!gid == "X") %>%
  dplyr::select(-GID) %>%
  tidyr::pivot_wider(., id_cols = gid, names_from = fx, values_from = values)

# Mean less than 1.75 (5% of total) and SE less than 0.1
# good_fit <- df %>%
#   mutate(g_fit = if_else(abs(mean) <= 1.75 & SE <= 0.1, 1, 0)) %>%
#   summarise(sum(g_fit))
# good_fit #74

# Mean less than 1.75 (5% of total) and SE less than 0.05
# good_fit <- df %>%
#   mutate(g_fit = if_else(abs(mean) <= 1.75 & SE <= 0.05, 1, 0)) %>%
#   summarise(sum(g_fit))
# good_fit #62

# Mean less than 1 and SE less than 0.12608
good_fit <- df %>%
  mutate(g_fit = if_else(abs(mean) <= 1 & SE <= 0.12608, 1, 0)) %>%
  summarise(sum(g_fit))
good_fit #66

# Mean less than 1 and SE less than 0.05
# good_fit <- df %>%
#   mutate(g_fit = if_else(abs(mean) <= 1 & SE <= 0.05, 1, 0)) %>%
#   summarise(sum(g_fit))
# good_fit #62

# Check for more restrictive fit
res_fit <- df %>%
  mutate(g_fit = if_else(abs(mean) <= 1 & SE <= 0.12608, 1, 0)) %>%
  dplyr::filter(g_fit == 1)

# Wilcoxon sign rank sum test
wilc <- read.csv(file = "Results_difs_df.csv") %>% 
  dplyr::select(-X) %>%
  select_if(is.numeric) %>%
  sapply(wilcox) %>%
  t() %>%
  data.frame() %>%
  dplyr::select(p.value) %>%
  mutate(Is_equal_to_zero = p.value >= 0.05) %>%
  dplyr::count(Is_equal_to_zero)

# Add excellent fit labels to shapefile
# Determine excellent fits - those grid cells whose means are not significantly different from zero
ex_fit <- read.csv(file = "Results_difs_df.csv") %>% 
  dplyr::select(-X) %>%
  select_if(is.numeric) %>%
  sapply(wilcox) %>%
  t() %>%
  data.frame() %>%
  dplyr::select(p.value) %>%
  mutate(ex_fit = if_else(p.value >= 0.05 | p.value == "NaN", 1, 0),
         gid = seq_along(ex_fit)) %>%
  dplyr::select(-p.value)
  

# Merge excellent fit information with rest of good fit information
sg50km_final <- merge(sg50km_final, ex_fit, all=TRUE) 

head(sg50km_final)

# Save Rdata file
save(sg50km_final, file = "sg50km_final.Rdata")

# Save shapefile
sf::st_write(sg50km_final, "sg50km_final.shp", 
             delete_dsn = TRUE) # Required if overwriting shapefile
