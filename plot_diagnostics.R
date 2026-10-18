# plot_diagnostics.R: Anna's Plotting_diagnostic_plots_for_site_selection.R with
# paths set for Jamie's processed/ folder. One PNG per site with time series of
# the met, flux and LAI variables (grey = gap-filled), in
# processed/PLUMBER2_site_selection_figs/. LAI_BOM (BOM 5 km, replaces Copernicus)
# is added by add_bom_lai_to_met.R.

library(ncdf4)
library(FluxnetLSM)


#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "C:/Users/jamie/extending-plumber2-flux-sites/processed/"

#Source function
source("C:/Users/jamie/extending-plumber2-flux-sites/PLUMBER2/functions/plot_timeseries.R")

#Helpers plot_timeseries.R uses that live inside FluxnetLSM
getMonthDays  <- FluxnetLSM:::getMonthDays
stripFilename <- function(f) basename(f)


#Variables to plot
vars <- c("SWdown", "LWdown", "Precip", "Tair", "Qair", "Wind", "CO2air", 
          "Rnet", "Qle", "Qh", "Qg", "Ebal", "LAI_MODIS", "LAI_BOM")
  
#Variable types
var_type <- c("Met", "Met", "Met", "Met", "Met", "Met", "Met",
              "Flux", "Flux", "Flux", "Flux", "Flux", "Met", "Met") 



#Find sites
site_files <- list.files(paste0(path, "/FLUXNET2015_sites/Nc_files/Met_with_LAI/"),
                          full.names=TRUE)

#Open file handles
site_nc_met <- lapply(site_files, nc_open)

#Get site codes
site_codes <- sapply(site_nc_met, function(x) ncatt_get(x, varid=0, "site_code")$value)


#Set output directory
outdir <- paste0(path, "/PLUMBER2_site_selection_figs/")
dir.create(outdir, showWarnings = FALSE)



#Loop through sites
for (s in 1:length(site_codes)) {

  cat(sprintf("[%d/%d] %s ... ", s, length(site_codes), site_codes[s]))
  
  
  #Get matching flux file
  nc_flux <- nc_open(list.files(paste0(path, "/FLUXNET2015_sites/Nc_files/Flux"), 
                                pattern=site_codes[s], full.names=TRUE))
  
  
  #Get time (used to var lengths etc)
  time <- ncvar_get(nc_flux, "time")
  
  #Get timing information
  timing <- GetTimingNcfile(nc_flux)  
  
  #QC information
  QC_measured <- 0
  if(substr(site_codes[s], 1, 3) == "AU-") QC_measured <- append(QC_measured, 10)
  
  
  #Set up file
  outfile <- paste0(outdir, "/", site_codes[s], "_timeseries_plot.png")
  png(outfile, height=1300, width=3500, res=50, pointsize=20)
  
  par(mfcol=c(3,5))
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    #Get data
    
    #Energy balance
    if (vars[v] == "Ebal") {
      
      rnet <- tryCatch(ncvar_get(nc_flux, "Rnet"), error=function(e) rep(0, length(time)))
      qle  <- tryCatch(ncvar_get(nc_flux, "Qle"), error=function(e) rep(0, length(time)))
      qh   <- tryCatch(ncvar_get(nc_flux, "Qh"), error=function(e) rep(0, length(time)))
      
      #Look for ground heat flux, set to zero if not found
      qg <- tryCatch(ncvar_get(nc_flux, "Qg"), error=function(e) rep(0, length(time)))
      
      #Calculate ratio (Rnet-G) / (Qle + Qh)
      var_data <- (rnet - qg) / (qle + qh)
      
      qc_data <- NA
      
      data_unit <- "-"
      
      
      #Cap values to 0-5
      var_data[var_data > 1.5] <- 1.5
      var_data[var_data < 0]  <- 0
      var_data[which(is.infinite(var_data))] <- 1.5
      
      
      
    #All other variables
    } else {
      
      #Set file for met variable
      if (var_type[v] == "Met") {
        nc <- site_nc_met[[s]]
        
      #Set file for  flux variable  
      } else if (var_type[v] == "Flux") { 
        nc <- nc_flux
      }
      
        #Get variable data
        var_data <- tryCatch(ncvar_get(nc, vars[v]), error=function(e) rep(NA, length(time)))
        
        #Get QC data
        qc_data <- tryCatch(ncvar_get(nc, paste0(vars[v], "_qc")), error=function(e) NA)
        
        #Get data units
        data_unit <- tryCatch(ncatt_get(nc, varid=vars[v], "units")$value, error=function(e) NA)
    
    }
      
    #Convert precip and Tair units
    if (vars[v] == "Precip") {
      
      var_data <- var_data * (time[2] - time[1])
      data_unit <- paste("mm/", (time[2] - time[1])/60, "min", sep="")
    }
    
    if (vars[v] == "Tair") {
      
      var_data <- var_data - 273.15
      data_unit <- "deg C"
    }
    
    

    #Extract QC data and replace all gap-filled values with 0
    # and measured with 1 (opposite to Fluxnet but what PALS expects)
    if (any(!is.na(qc_data))) {
      
      
      qc_data[!(qc_data %in% QC_measured)]  <- 2 #replace gap-filled values with a temporary value
      qc_data[qc_data %in% QC_measured]     <- 1 #set measured to 1
      qc_data[qc_data == 2] <- 0 #set gap-filled to 0
      
      #If first value missing, set to measured (to avoid an error when PALS
      #checks if first value -1)
      if(is.na(qc_data[1])){ qc_data[1] <- 0}
      
      
      #Else set to PALS option corresponding to no QC data
    } else {
      qc_data <- matrix(-1, nrow = 1, ncol = 1)
    }
    
    
  
    #Plot file
    Timeseries(obslabel=site_codes[s], tsdata=as.matrix(var_data),
               varname=vars[v],
               ytext=paste(vars[v], " (", data_unit, ")", sep=""), 
               legendtext=vars[v],
               plotcex=2 , timing=timing, 
               smoothed = FALSE, winsize = 1, 
               plotcolours="black",
               vqcdata = as.matrix(qc_data),
               na.rm=TRUE)
  
  
  
  } #vars
  
  dev.off()
  nc_close(nc_flux)
  cat("saved", basename(outfile), "\n")
  
} #sites

cat("\nDone. Figures in", outdir, "\n")
