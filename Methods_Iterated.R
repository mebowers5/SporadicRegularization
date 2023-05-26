############################################################################################################
############## Create random tracks and test biases in the distribution summaries that result ##############
############## from crawl movement models and pathroutr re-routing off the US East Coast ###################
############################################################################################################
##############################################

#Remove everything from environment
rm(list=ls())

# Clear cache
gc()

# Increase R memory limit
# library(usethis) 
# usethis::edit_r_environ()

# If desired, change the temporary directory for R
# write("TMP = 'C:\Program Files\R\R-4.2.3\rtemp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

library(dplyr)
library(sf)
library(sp)

#########################################################
# Setup parallel processing
#########################################################

# Find out how many cores you have
library(parallel)

## Avoid using "logical cores". Apparently this can actually wind up
## slowing things down.
numCores <- detectCores(logical = F)
numCores

# Create cluster
cl <- makeCluster(numCores-4, outfile = "cluster_out")

# Load libraries
clusterEvalQ(cl, {
  
  library(sf)
  library(sp)
  library(glatos)
  library(dplyr)
  library(momentuHMM)
  library(crawl)
  library(lubridate)
  library(pathroutr)
  library(sfnetworks)
  library(tidyverse)
  library(stplanr)
  library(future)
  library(doFuture)
  
})

# Set working directory to file destination of iteration objects
setwd("E:/IterObjects")
# setwd("/Volumes/Bowers_PhD/IterObjects/")

# Objects required for iterations to run
study_site <- sf::st_read("fo_Bathy500m.shp")

# Convert sf study site into sp object
study_site <- study_site %>% as(., 'Spatial') 

crs <- "+init=epsg:3857"

# Set CRS in study site
sp::proj4string(study_site) <- sp::CRS(crs)

# Call in release site(s)
rel_site <- sf::st_read("m_fo_rel_site.shp")

# we need to get our land polygon into a proper format; essentially, we want it to be
# a series of POLYGONs (and not MULTIPOLYGONs or GEOMETRYCOLLECTION).
land_barrier <- sf::st_read("fo_land_barrier.shp") %>%
  sf::st_transform(3857) %>%
  #sf::st_crop(bb) %>%
  sf::st_collection_extract('POLYGON') %>%
  sf::st_cast('POLYGON')

# Call in receiver stations
stations <- sf::st_read("fo_stations.shp") 

# Set working directory back to where other data may be pulled from and where items should be saved
setwd("E:/Data/")
# setwd("/Volumes/Bowers_PhD/Data/reg_spoR_Data")

# Call in grid to analyze counts
grid_res <- function(km, study_site) {
  grid_spacing <- km * 1000 # Desired kilometers times 1000 m per 1 km
  
  HSgrid <- sf::st_make_grid(study_site, square = T, cellsize = c(grid_spacing, grid_spacing)) # Create a grid inside the coastal 500 m isobath polygon
  
  HSgrid <- sf::st_intersection(sf::st_as_sf(study_site), HSgrid) %>%
    sf::st_as_sf() %>%
    mutate(gid = seq_along(geometry)) %>%
    sf::st_make_valid()
  
  return(HSgrid)
}

HS_10km_grid <- grid_res(10, study_site)
HS_25km_grid <- grid_res(25, study_site)
HS_50km_grid <- grid_res(50, study_site)
HS_100km_grid <- grid_res(100, study_site)

# Create an empty data frame and save to .csv file
# csv_cols <- c("iteration", "gid", "sim_count", "mod_count", "dif")
# 
# csv_output <- data.frame(matrix(nrow = 1, ncol = length(csv_cols))) %>%
#   rlang::set_names(csv_cols) %>%
#   mutate(iteration = as.numeric(0),
#          gid = as.numeric(0),
#          sim_count = as.numeric(0),
#          mod_count = as.numeric(0),
#          dif = as.numeric(0))
# 
# csv_filename <- "difs_df_100km.csv"
# 
# write.csv(csv_output, file = csv_filename)

# Create a spatial network graph ("visibility graph") outside of the iterations
# prt_visgraph will build our visual graph network from our land barrier object and 
# return a SpatialLinesNetwork / sfNetwork that has no edges that cross land
vis_graph <- pathroutr::prt_visgraph(land_barrier)

# Parameters in simulated data
total_iterations <- 200 # Total number of iterations desired
reps <- 200 # Number of replicates that will result in desired iterations given reps
runs <- 1:(total_iterations/reps) # Interval at which outputs should be written to file
anims <- 30   # Number of animals (35 for example data)
yr <- 1       # Number of years (3 for example data)
theta <- c(0, 1.74)  # Turning angle
vmin <- 0.98  # Lower bound of mean recorded velocity for animal in m/s (ex: blacktip shark) (average - standard deviation)
vmax <- 1.58  # Upper bound of mean recorded velocity for animal in m/s (ex: blacktip shark) (average + standard deviation)
rel_site <- rel_site  # Intersection of bounding box of release site(s) and study site 
n_days <- 365*yr # Number of days animal is tracked
initHeading <- 0  # Initial heading of animal

# Create list of different grid cell resolutions
HSgrid <- list(HS_100km_grid, HS_50km_grid, HS_25km_grid, HS_10km_grid)

# # Create empty list object and save it
# df <- list()

# Assign the filename
filename <- "results_post_crawlupdate.Rdata"

# # Save the empty list object to the filename
# save(df, file = filename)

# Remove unnecessary objects before sending mirroring environment to cluster
rm(list = c("HS_100km_grid", "HS_50km_grid", "HS_25km_grid", "HS_10km_grid", "df", 
            "fail_test", "failed_trks", "mod_trks", "sim_trks", "test", "grid_res"))

# Mirror global environment objects on all cores
parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))

################## Start iterations ##########################
system.time(
  # pryr::mem_used(
  base::replicate(n = reps, expr = {
    
    # Clear cache
    gc()
    
    iterated_results <- parallel::parLapply(cl, runs, fun = function(I) {
      
      # Simulate tracks with wrapper function using cluster
      df <- simul_trks(anims, study_site, theta, vmin, vmax, rel_site, n_days, initHeading)
      
      # Model tracks from simulated tracks and compare counts across multiple grid cell resolutions
      df <- comp_trks(df, stations, land_barrier, vis_graph, HSgrid, multi.grid = TRUE)
      
      return(df)
    }
    )
    
    load(file = filename)
    
    df <- append(df, iterated_results)
    
    save(df, file = filename)
    
    # Clear cache
    gc()
    
    return(df)
  }
  )
# )
)

# Stop the cluster
parallel::stopCluster(cl)
