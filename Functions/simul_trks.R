# Simulate tracks
simul_trks <- function(anims, study_site, theta, vmin, vmax, rel_site, n_days, initHeading){
  library(sp)
  
  sim <- base::replicate(n = anims, expr = {
    
    # Clear cache
    gc()
    
    initPos <- sf::st_sample(rel_site$geometry, size=1, type ="random", exact = TRUE) %>%
      as.data.frame() %>%
      mutate(lon = sf::st_coordinates(geometry)[,1],
             lat = sf::st_coordinates(geometry)[,2]) 
    
    stepLen <- as.numeric(sample(vmin:vmax, 1)*60*60*24/24) # Sample between minimum and maximum velocity of blacktip sharks to set step length per hour for each individual
    
    simu <- glatos::crw_in_polygon(study_site, 
                           theta = theta, 
                           stepLen = stepLen,
                           initPos = c(initPos$lon, initPos$lat), 
                           nsteps = n_days*24, initHeading = initHeading) %>% 
      sf::st_as_sf()
    
    return(simu)
  }
  )
  
  simu <- sim %>% as.data.frame() %>% 
    setNames(gsub("geometry.", perl = TRUE, "", names(.))) %>%
    setNames(gsub("geometry", perl = TRUE, "0", names(.))) %>%
    tidyr::pivot_longer(., cols = everything(), names_to = "AnimalID", values_to = "geom") %>% 
    mutate(ID = as.numeric(AnimalID) + 1) %>%
    group_by(ID) %>%
    mutate(uid = seq_along(AnimalID)) #Assign sequential numbers to rows to preserve directionality
  
  simu <- simu %>%  
    mutate(POINT_X = sf::st_coordinates(geom)[,1], POINT_Y = sf::st_coordinates(geom)[,2]) %>%  #Obtain x and y coordinates so that you can force northward and southward movement
    dplyr::select(uid, 
                  ID, 
                  geom, POINT_Y, POINT_X) %>%
    group_by(ID) %>%
    arrange(uid) %>%
    mutate(time = seq.POSIXt(from = as.POSIXct("2000-01-01 01:00:00"), by = "1 hour", along.with = uid),   #Make up dates for the row numbers by Iteration and Animal ID
           start_x = POINT_X, start_y = POINT_Y, end_x = lead(POINT_X), end_y = lead(POINT_Y)) %>% #Prep data for making lines - ensures proper order
    sf::st_as_sf()
  
  #Make lines and unite geometries by uid's so that each line segment maintains time information
  simu <- simu %>%
    filter(!is.na(end_y)) %>%
    tidyr::nest() %>% 
    mutate(
      data = purrr::map(data, 
                        ~ mutate(.x,
                                 x = purrr::pmap(.l = list(start_x, start_y, end_x, end_y),
                                                 .f = make_line)))) 
  
  return(simu)
}
