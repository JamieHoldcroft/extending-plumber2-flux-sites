#####################################################################
# add_bom_lai_to_met.R
#
# Anna's Step3_Process_Copernicus_LAI_for_all_sites.R, run on the
# BOM 5 km monthly LAI (lai_bom_<year>_Jul2025.nc, 2000-2024) instead
# of Copernicus. Adds a variable LAI_BOM to every met file in
# Met_with_LAI (the folder the MODIS script fills), like Step3 adds
# LAI_Copernicus.
#
# What it does, in order:
#  1. Opens every met file in Met_with_LAI and reads each site's lat/lon.
#  2. Reads the BOM LAI at each site, one year at a time
#     (5 km pixel the site is in; if that pixel is empty, e.g. a
#     coastal site, the nearest non-empty pixel, as Anna's
#     sample_raster_NA does).
#     -> replaces Anna's lines 75-137 (Copernicus 3x3 averaging), as Anna asked.
#  3. Fills the 5 months the BOM files have no data for
#     (Jan 2000 and Sep-Dec 2024) with that month's mean over all years.
#  4. From here on Anna's code (lines 145-342): smooth with a spline,
#     make a 12-month climatology, add a 12-month rolling anomaly,
#     fill years outside 2000-2024 with the climatology, repeat each
#     monthly value for every half-hour of the month, write LAI_BOM.
#
# Every change from Anna's script is marked "CHANGED:".
#####################################################################

library(ncdf4)
library(zoo)

#clear R environment
rm(list=ls(all=TRUE))

#CHANGED: TRUE = linear interpolation between mid-month values (Anna's request),
#so LAI changes gradually instead of jumping on the 1st of each month.
#FALSE = Anna's original (each monthly value repeated over the whole month).
interpolate_lai <- TRUE

#CHANGED: TRUE = replace LAI_BOM in files that already have it (needed to switch
#interpolate_lai). FALSE = skip those files.
overwrite_lai_bom <- TRUE

#Set path
path <- "C:/Users/jamie/extending-plumber2-flux-sites/processed/"

#Which fluxnet data to use?
flx_path <- paste0(path, "/FLUXNET2015_sites/")   #CHANGED: Jamie's path


#Met files with processed LAI time series (made by add_lai_to_met.R, the MODIS step)
outdir <- paste0(flx_path, "/Nc_files/Met_with_LAI/")

#CHANGED: this script can run before or after add_lai_to_met.R (MODIS). It makes
#Met_with_LAI if needed and copies in any Met file that isn't there yet, never
#overwriting a file that already has LAI in it.
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
file.copy(from = list.files(paste0(flx_path, "/Nc_files/Met/"), full.names=TRUE),
          to = outdir, overwrite=FALSE)

#Get sites
site_files <- list.files(outdir, full.names=TRUE)

#Open file handles
site_nc <- lapply(site_files, nc_open, write=TRUE)

#Get site codes
site_codes <- sapply(site_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)

cat("Found", length(site_nc), "met files\n")


### Get LAI data ###

#Set LAI path   CHANGED: BOM files
lai_path <- "C:/Users/jamie/extending-plumber2-flux-sites/BOM_LAI"

#Get files
lai_files <- list.files(lai_path, pattern="lai_bom_.*_Jul2025.nc", full.names=TRUE)

#Should span 2000-2024   CHANGED: was 19 files 1999-2017
if(length(lai_files) != 25) stop("incorrect LAI files found")


#Set BOM start and end year   CHANGED: was copernicus_startyr 1999 / endyr 2017
#(kept Anna's variable names so the rest of her code is unchanged)
copernicus_startyr <- 2000
copernicus_endyr   <- 2024


#Get site coordinates
lat <- lapply(site_nc, ncvar_get, varid="latitude")
lon <- lapply(site_nc, ncvar_get, varid="longitude")

#CHANGED: plain numbers (ncvar_get returns 1x1 arrays here)
lat <- lapply(lat, as.numeric)
lon <- lapply(lon, as.numeric)


##############################################################
### CHANGED: replaces Anna's lines 75-137                   ###
### Read the BOM pixel at each site, one year at a time      ###
##############################################################

#BOM grid (same in every file: lat -10 to -44, lon 112 to 154, 0.05 degrees)
nc_lai   <- nc_open(lai_files[1])
lai_lats <- ncvar_get(nc_lai, "latitude")
lai_lons <- ncvar_get(nc_lai, "longitude")
nc_close(nc_lai)

#Pixel each site falls in (nearest pixel centre, same as raster::extract)
site_i <- sapply(lon, function(x) which.min(abs(lai_lons - x)))
site_j <- sapply(lat, function(x) which.min(abs(lai_lats - x)))

#Initialise
site_tseries <- lapply(1:length(site_nc), function(x) vector())

pb <- txtProgressBar(min=0, max=length(lai_files), style=3)

for (y in 1:length(lai_files)) {

  nc_lai <- nc_open(lai_files[y])

  for (s in 1:length(site_nc)) {

    #All 12 months for the site's pixel (lai is lon x lat x time in R)
    vals <- ncvar_get(nc_lai, "lai", start=c(site_i[s], site_j[s], 1), count=c(1, 1, -1))

    #If the pixel is empty in every month (sea or masked), use the nearest
    #non-empty pixel within 3 degrees (as Anna's sample_raster_NA does)
    if (all(is.na(vals))) {

      i_win <- which(abs(lai_lons - lon[[s]]) <= 3)
      j_win <- which(abs(lai_lats - lat[[s]]) <= 3)
      block <- ncvar_get(nc_lai, "lai", start=c(min(i_win), min(j_win), 1),
                         count=c(length(i_win), length(j_win), -1))

      #Use a month that has data (some months are empty everywhere)
      has_data <- which(apply(block, 3, function(m) any(!is.na(m))))[1]
      ok       <- which(!is.na(block[, , has_data]), arr.ind=TRUE)
      dist     <- (lai_lons[i_win][ok[, 1]] - lon[[s]])^2 + (lai_lats[j_win][ok[, 2]] - lat[[s]])^2
      best     <- ok[which.min(dist), ]
      vals     <- block[best[1], best[2], ]

      if (y == 1) cat("\n", site_codes[s], ": pixel empty, using nearest pixel at",
                      lai_lats[j_win][best[2]], lai_lons[i_win][best[1]], "\n")
    }

    #Append year values
    site_tseries[[s]] <- append(site_tseries[[s]], vals)
  }

  nc_close(nc_lai)
  setTxtProgressBar(pb, y)
}
close(pb)


#CHANGED: BOM files have no data for Jan 2000 and Sep-Dec 2024 (empty
#everywhere). Fill those months with that calendar month's mean over all
#years, otherwise smooth.spline below fails for every site.
for (s in 1:length(site_tseries)) {

  month_of_step <- rep(1:12, length.out=length(site_tseries[[s]]))
  missing       <- which(is.na(site_tseries[[s]]))

  for (m in missing) {
    site_tseries[[s]][m] <- mean(site_tseries[[s]][month_of_step == month_of_step[m]], na.rm=TRUE)
  }
}


#Check that no missing LAI values (Anna's line 135)
if (any(sapply(site_tseries, function(x) any(is.na(x))))) {
  stop("Missing LAI values present")
}


##########################
### Loop through sites ###
##########################

summary_table <- data.frame()

for (s in 1:length(site_nc)) {

  print(paste0("Processing site ", s, "/", length(site_nc), " (", site_codes[s], ")"))

  #CHANGED: skip if LAI_BOM already added (re-running would crash in ncvar_add)
  lai_bom_exists <- "LAI_BOM" %in% names(site_nc[[s]]$var)
  if (lai_bom_exists & !overwrite_lai_bom) {
    print("  LAI_BOM already in file, skipping")
    nc_close(site_nc[[s]])
    next
  }

  #Set time stamps (do inside loop as gets adjusted below)
  lai_time <- seq.Date(from=as.Date(paste0(copernicus_startyr, "-01-01")), by="month",   #CHANGED: was "1999-01-01"
                       length.out=length(site_tseries[[s]]))


  #Smooth LAI time series with spline (and cap negative values)
  smooth_lai_ts = tryCatch(smooth.spline(lai_time, site_tseries[[s]])$y,
                           error= function(e) NA)
  if (all(is.na(smooth_lai_ts))) {
    warning(paste0("could not process site", site_codes[s], ", missing values in LAI"))
    next
  }

  smooth_lai_ts[smooth_lai_ts < 0] <- 0



  ##########################
  ### Create climatology ###
  ##########################


  #Copernicus has 12 time steps per year
  no_tsteps <- length(which(grepl(copernicus_startyr, lai_time)))

  #Initialise
  copernicus_clim <- vector(length=no_tsteps)

  for (c in 1:no_tsteps) {

    #Indices for whole years
    inds <- seq(c, by=no_tsteps, length.out=length(lai_time)/no_tsteps)

    #Calculate average for time step
    copernicus_clim[c] <- mean(smooth_lai_ts[inds])

  }




  ###################################
  ### Calculate running anomalies ###
  ###################################

  #Initialise
  lai_clim_anomalies <- rep(NA, length(lai_time))

  #Repeat climatology for whole time series
  copernicus_clim_all <- rep_len(copernicus_clim, length(lai_time))


  #Calculate running mean anomaly (+/- 6 months either side of each time step)
  anomaly <- rollmean(smooth_lai_ts - copernicus_clim_all, k=12, fill=NA)

  #Add rolling mean anomaly to climatology
  lai_clim_anomalies <- copernicus_clim_all + anomaly


  #Check if remaining NA values from missing time steps, gapfill if found
  if (any(is.na(lai_clim_anomalies))) {

    #Find missing values
    missing <- which(is.na(lai_clim_anomalies))

    #Repeat climatology for all years and gapfill time series
    clim_all_yrs <- rep(copernicus_clim, floor(length(lai_time)/no_tsteps))
    lai_clim_anomalies[missing] <- clim_all_yrs[missing]

  }


  ###################################
  ### Match with site time series ###
  ###################################


  #Each year has 12 time steps

  #Get timing info for site
  site_start_time <- ncatt_get(site_nc[[s]], "time")$units
  site_time       <- ncvar_get(site_nc[[s]], "time")
  site_tstep_size <- 86400 / (site_time[2] - site_time[1])

  #Extract year
  startyr    <- as.numeric(substr(site_start_time, start=15, stop=18))
  obs_length <- length(site_time)
  nyr        <- round(obs_length/(site_tstep_size*365))
  endyr      <- startyr + nyr - 1



  #Add climatological values to smoothed lai time series

  #Overwrite original time series with new extended data

  #If start year earlier than Copernicus
  if (startyr < copernicus_startyr) {

    #LAI data
    lai_clim_anomalies <- append(rep(copernicus_clim, copernicus_startyr - startyr),
                                 lai_clim_anomalies)

    #Time vector
    lai_time <- append(seq.Date(as.Date(paste0(startyr, "-01-01")), by="month",
                       length.out=no_tsteps * (copernicus_startyr - startyr)),
                       lai_time)

  }

  #If end year later than Copernicus
  if (endyr > copernicus_endyr) {

    #LAI data
    #CHANGED: Anna's line was
    #  append(smooth_lai_ts, rep(lai_clim_anomalies, endyr - copernicus_endyr))
    #which drops the anomaly step and pads with the whole series from the
    #start (2000...) instead of the climatology. Now: climatology, as the
    #comment "Add climatological values" and the start-year branch do.
    lai_clim_anomalies <- append(lai_clim_anomalies,
                                 rep(copernicus_clim, endyr - copernicus_endyr))


    #Time vector
    lai_time <- append(lai_time,
                       seq.Date(as.Date(paste0(copernicus_endyr+1, "-01-01")), by="month",
                       length.out=no_tsteps * (endyr - copernicus_endyr)))

  }



  #Find modis time step corresponding to site start time
  start_ind <- which(lai_time == paste0(startyr, "-01-01"))
  end_ind   <- tail(which(grepl(endyr, lai_time)), 1) #Last index of end year

  #Extract MODIS time steps matching site
  copernicus_ts_for_site   <- lai_clim_anomalies[start_ind:end_ind]
  copernicus_time_for_site <- lai_time[start_ind:end_ind]


  #Repeat Copernicus time series to create a time series matching site time step
  copernicus_tseries <- vector()


  #Loop through time steps
  for (t in 1:length(copernicus_time_for_site)) {

    #Last time step
    if (t == length(copernicus_time_for_site)) {

      #Use the number of time steps that ensures final time series matches the length of site data
      copernicus_tseries <- append(copernicus_tseries, rep(copernicus_ts_for_site[t],
                                   length(site_time) - length(copernicus_tseries)))

      #All other time steps
    } else {

      time_diff <- copernicus_time_for_site[t+1] - copernicus_time_for_site[t]

      #Repeat each days estimate by the number of days and time steps per day
      copernicus_tseries <- append(copernicus_tseries, rep(copernicus_ts_for_site[t],
                                   as.numeric(time_diff * site_tstep_size)))

    }

  }


  #Check that the number of time steps match
  if (length(copernicus_tseries) != length(site_time)) stop("MODIS and site time steps don't match")

  #CHANGED: linear interpolation. Each monthly value is placed at the middle of
  #its month (the half-hour halfway through the block it was repeated over above),
  #then straight lines are drawn between neighbouring months. Before the first
  #and after the last mid-month point the value is held flat (rule=2).
  if (interpolate_lai) {
      n_per_mon <- c(sapply(1:(length(copernicus_time_for_site)-1), function(t)
                   as.numeric(copernicus_time_for_site[t+1] - copernicus_time_for_site[t]) * site_tstep_size))
    n_per_mon <- c(n_per_mon, length(site_time) - sum(n_per_mon))
    mid_step  <- cumsum(n_per_mon) - n_per_mon + n_per_mon/2 + 0.5
    copernicus_tseries <- approx(x=mid_step, y=copernicus_ts_for_site,
                                 xout=1:length(site_time), rule=2)$y
  }

  #CHANGED: also check no missing values (as Step2 does)
  if (any(is.na(copernicus_tseries))) stop("Missing values in final BOM LAI time series")



  ######################################
  ### Add LAI time series to NC file ###
  ######################################

  #Save to file

  # Define variable:   CHANGED: name LAI_BOM and BOM long name
  laivar = ncvar_def('LAI_BOM', '-', list(site_nc[[s]]$dim[[1]], site_nc[[s]]$dim[[2]], site_nc[[s]]$dim[[3]]),
                     missval=-9999, longname='BOM 5 km monthly leaf area index (Jul2025 release)')
  # Add variable and then variable data:
  #CHANGED: only add the variable if it is not already there (overwrite_lai_bom)
  if (!lai_bom_exists) site_nc[[s]] = ncvar_add(site_nc[[s]], laivar)
  ncvar_put(site_nc[[s]], 'LAI_BOM', copernicus_tseries)

  #Close file handle
  nc_close(site_nc[[s]])

  summary_table <- rbind(summary_table, data.frame(site=site_codes[s], years=paste0(startyr, "-", endyr),
                                                   min=round(min(copernicus_tseries), 2),
                                                   max=round(max(copernicus_tseries), 2)))

}

print(summary_table, row.names=FALSE)
