iters <- function(reps, cl, runs, anims, yr, npts, vis_graph, land_barrier, HSgrid, 
                  timeStep, fixPar, theta = c(0,0), attempts, retryFits, all) {
  base::replicate(n = reps, expr = {
    all_difs <- parallel::parLapply(cl, runs, fun = function(I) {
      
      ### Populate random points to create tracks that are coerced to move north and south
      # If the variability itself is a problem, try set.seed
      # If that doesn't work try tryCatch(wrap everything that you want to "catch")
      #Create random tracks
      sim_trks <- study_site %>% 
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
      #plot(study_site$geometry, col="white")
      #plot(rand$geom, size = 0.5, col="blue", add=TRUE)
      
      #Make lines and unite geometries by uid's so that each line segment maintains time information
      sim_trks <- sim_trks %>%
        filter(!is.na(end_y)) %>%
        tidyr::nest() %>% 
        mutate(
          data = purrr::map(data, 
                            ~ mutate(.x,
                                     x = purrr::pmap(.l = list(start_x, start_y, end_x, end_y),
                                                     .f = make_line))))
      
      # Maintain detailed time info in lines to derive acoustic telemetry detection data later
      mod_trks <- sim_trks %>%
        tidyr::unnest(data) %>%
        mutate(x = st_sfc(x)) %>% 
        st_as_sf(sf_column_name = 'x', crs = 3857)
      
      #Unite geometries by AnimalID's so that you have complete tracks/lines summarized by ID for the grid cell count later
      sim_trks <- sim_trks %>% 
        mutate(
          data = purrr::map(data,
                            ~ mutate(.x, x = st_sfc(x))),
          x = purrr::map(data, ~ st_union(st_set_geometry(.x, 'x'))), ##Preserves order
          x = purrr::map(x, ~ st_cast(.x, 'MULTILINESTRING')),
          ID = AnimalID  
        ) %>% 
        dplyr::select(-data) %>% 
        tidyr::unnest(x) %>% 
        st_as_sf(sf_column_name = 'x', crs = 3857)
      
      ### Find intersections of tracks and receiver stations
      # These intersections represent data that a researcher might receive from an animal tagged with an acoustic transmitter that traveled along the track that we created
      
      # Derive acoustic detection data through the intersection of lines and receiver ranges
      # Cannot just take the intersection because the lines self intersect
      # Take the symmetrical difference between the lines and self intersection to get the true derived detection data
      mod_trks <- sf::st_intersection(mod_trks, stations$geometry) %>%  # x provides intersection geometry - linestrings that cross station polygons
        st_centroid() # Take the centroid of the linestrings
      
      mod_trks <- st_snap(mod_trks, st_centroid(stations$geometry), tolerance = 650) # Snap the linestring centroids to the station point geometry - geometry 'x' is altered
      
      ### Run "derived acoustic telemetry data through continuous time correlated random walk model
      #Load data----
      
      #In real data, problem caused by the zero movement between locations at the same acoustic receiver, close together in time - solved by averaging daily locations
      #Average daily locations
      mod_trks <- mod_trks %>%
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
      crw_mod <- momentuHMM::crawlWrap(obsData = mod_trks, timeStep = "15 min", fixPar = c(NA, NA), theta = c(0,0),
                                       attempts = 5, retryFits = 0)
      # plot(crw_mod)
      
      # Clear cache
      gc()
      
      #Continue with movement model data products
      #Get the data.frame of predicted locations
      mod_trks <- data.frame(crw_mod$crwPredict)
      
      #Turn crwPredict product into data frame
      mod_trks <- mod_trks %>%
        dplyr::select(ID, locType, time, mu.x, mu.y, speed) %>%
        mutate(locType = ifelse(is.na(locType), "p", locType)) %>%  #Specify predicted locs
        as.data.frame
      
      cat("Reading in track data")
      
      # Read in the track data; here, I'm filtering to just include the predicted locations
      # and making sure the time is a proper POSIX. Neither of those two steps are
      # necessary
      mod_trks <- mod_trks %>%
        #dplyr::filter(locType == "p") %>% 
        dplyr::select(ID, time, mu.x, mu.y) %>%
        dplyr::mutate(time = lubridate::ymd_hms(time), ID = as.factor(ID))
      
      cat("Converting the track data to CRS 3857")
      
      # convert the track_data to sf and set the CRS; the bb step is just a way to limit
      # the size of the land polygon and save some computation time when creating vis_graph
      mod_trks <- mod_trks %>% sf::st_as_sf(coords = c("mu.x","mu.y"), crs = 3857)
      #bb <- sf::st_as_sfc(sf::st_bbox(track_path))
      
      cat("Nesting the tracks")
      # there are multiple paths identified by ID; we'll group and nest for a proper
      # tidyverse/list-column workflow
      mod_trks <- mod_trks %>%
        group_by(ID) %>%
        tidyr::nest()
      
      # Clear cache
      gc()
      
      cat("Trimming points from the land_barrier")
      
      # the track cannot start or end within the land barrier; prt_trim() trims those out
      mod_trks <- mod_trks %>%
        rowwise() %>%
        mutate(trim_data = list(pathroutr::prt_trim(data, land_barrier)))
      
      cat("Creating re-routed points")
      
      # here, we create our re-routed points; the return is a two column data frame with the
      # index location in the original point data and the new geometry. The user can handle
      # updating of those original point data or pass the result on to prt_update_points()
      mod_trks <- mod_trks %>% dplyr::rowwise() %>%
        mutate(rrt_pts = list(prt_reroute(trim_data, land_barrier, vis_graph)))
      
      # Clear cache
      gc() 
      
      cat("Updating paths with re-routed points")
      
      # NOTE: previous versions of prt_update_points() had the argument order reversed from
      # what it now requires. The updated geometry points are passed first (here, `rrt_pts`)
      # and, then, the original data to be updated. This order should allow for easy piping
      # from prt_reroute()
      mod_trks <- mod_trks %>% dplyr::rowwise() %>%
        mutate(path_pts = list(prt_update_points(rrt_pts, trim_data)),
               path_lines = list(path_pts %>% summarise(do_union = FALSE) %>% sf::st_cast('LINESTRING')))  # do_union MUST be FALSE!
      
      # Clear cache
      gc() 
      
      cat("Compiling all points and lines into single object")
      
      # we need to rbind all of our lines and points into single objects that can be plotted
      # Need to maintain "ID" column here - THIS WORKS!
      mod_trks$geom <- do.call(rbind, mod_trks$path_lines)
      mod_trks$geom <- sf::st_set_crs(mod_trks$geom, 3857)
      
      mod_trks <- mod_trks %>% dplyr::select(ID, geom) 
      
      cat("Convert rowwise df to sf object")
      
      ### Create grid that maintains row for each grid cell
      #Convert from rowwise_df to sf object
      mod_trks <- mod_trks %>%
        rowwise() %>%
        mutate(geom = sf::st_geometry(geom)) %>%
        dplyr::select(ID, geom) %>%
        ungroup() %>%
        st_as_sf(., sf_column_name = "geom")
      
      cat("Counting animals per grid cell in derived data")
      
      ### Count distinct AnimalIDs in each grid cell
      # Count IDs per grid cell for derived data
      mod_count <- sf::st_join(HSgrid, mod_trks, join = st_intersects)%>%
        distinct(gid, ID, geometry)
      
      mod_count <- aggregate(ID ~ gid, data = mod_count, FUN = length) %>%
        rename(mod_count = ID)
      
      cat("Counting animals per grid cell in complete data")
      
      # Count IDs per grid cell for complete data
      sim_count <- sf::st_join(HSgrid, sim_trks, join = st_intersects)%>%
        distinct(gid, ID, geometry)
      
      sim_count <- aggregate(ID ~ gid, data = sim_count, FUN = length) %>%
        rename(sim_count = ID) 
      
      cat("Join count data to all grid cells file")
      
      #Now that we have every combination of gid and Iteration in one variable, merge one of the count files to it - doesn't really matter which on but we'll use the complete data one here
      # Note: Is the above comment wrong or is the code wrong? I think the comment is wrong - leftover from when iterations were assigned instead of looped.
      tc <- left_join(all, sim_count, by = "gid") %>%
        mutate_all(~replace(., is.na(.), 0))  # Replace NA's with zeros so we can check count
      
      # Check for code accuracy - Sum should equal # of observations in tc_comp
      #sum(tc_comp$sim_count > 0)
      
      cat("Join derived count data to complete count data")
      
      #Join derived count to data frame
      tc <- left_join(tc, mod_count, by = "gid") %>%
        mutate_all(~replace(., is.na(.), 0))  # Replace NA's with zeros so we can check count
      
      cat("Find the difference in counts")
      
      #Calculate differences in each grid cell by iteration and gid
      tc <- tc %>%
        dplyr::select(gid, sim_count, mod_count) %>%
        mutate(sim_count = as.double(sim_count),
               mod_count = as.double(mod_count),
               dif = sim_count - mod_count) %>%
        arrange(gid) %>%
        dplyr::select(gid, sim_count, mod_count, dif) %>%
        as.data.frame
      
      return(tc)
      
      # Clear cache
      gc()
      
    }
    )
    ################# After loops finish ##################
    all_difs # Returns list return object from all runs per each node
    
    difs_df <- read.csv(file = "difs_df_25km.csv", check.names = FALSE) %>%
      dplyr::select(-1)
    
    #Turn list into data frame
    difs2add_df <- sapply(all_difs, '[', 1:max(sapply(all_difs, length))) %>%
      as.data.frame %>%
      t() %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "iteration") %>%
      tidyr::unnest(cols = c("gid", "sim_count", "mod_count", "dif")) %>%
      mutate(iteration = as.numeric(gsub("V", "", .$iteration)) + max(difs_df$iteration))
    
    # Append results to others
    difs_df <- rbind(difs_df, difs2add_df)
    
    # Save the results
    write.csv(difs_df, file = "difs_df_25km.csv")
    
    # Clear cache
    gc()
    
  }
  )
}