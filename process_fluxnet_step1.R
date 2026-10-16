# run_step1_fluxnet.R
#
# The FLUXNET2015 section of Anna's Step1_Process_all_available_flux_sites_for_PLUMBER2.R
# (lines 46-160), run on the new FLUXNET Shuttle zips in Jamie's FLUXNET folder.
#
# What it does:
#   Part 0 (new, needed because the Shuttle files have new names):
#     - opens each TERN_<site>_FLUXNET_..._v1.3_r1.zip and takes out only the two
#       files Step1 uses: FLUXMET_HH (flux + met data) and ERA5_HH (ERA5 weather,
#       used to gap-fill met). The zips are not changed.
#     - saves them under the FLUXNET2015 names Anna's functions look for, e.g.
#       FLX_AU-Ync_FLUXNET2015_FULLSET_HH_2012-2025_1-3.csv
#       FLX_AU-Ync_FLUXNET2015_ERAI_HH_1981-2025_1-3.csv
#   Part 1 (Anna's code, same functions, same order, same settings):
#     - get_fluxnet_files, get_path_site_code, get_fluxnet_version_no,
#       get_fluxnet_erai_files, then convert_fluxnet_to_netcdf with
#       met_gapfill = "ERAinterim", flux_gapfill = "statistical", 100% thresholds,
#       min_yrs = 1, check_range_action = "warn", include_all_eval = TRUE
#   Differences from Anna's run (processing is identical):
#     - one site at a time instead of a 12-core cluster, with a progress bar
#     - errors are printed and collected instead of silently dropped (Anna's
#       tryCatch returns NULL), and a summary table is saved at the end
#     - Hourly and Halfhourly, like Anna's tstep loop: hourly sites (AU-Tum,
#       AU-Otw) only have FLUXMET_HR/ERA5_HR files, which go to the Hourly folders
#
# Run with: source("run_step1_fluxnet.R", echo = TRUE)

library(FluxnetLSM)  # convert_fluxnet_to_netcdf

#clear R environment
rm(list=ls(all=TRUE))


################
### Settings ###
################

#Folder with the TERN_*_FLUXNET_*.zip files
zip_path <- "C:/Users/jamie/extending-plumber2-flux-sites/FLUXNET"

#Set output path for all data
out_path <- "C:/Users/jamie/extending-plumber2-flux-sites/processed"

#Unzipped and renamed inputs go here (as Anna: one folder per time step)
#CHANGED: Hourly added
tstep    <- c("Hourly", "Halfhourly")
res_code <- c(Hourly="HR", Halfhourly="HH")
in_paths  <- paste0(out_path, "/FLUXNET_Original_data/", tstep)
era_paths <- paste0(out_path, "/FLUXNET_Original_data/EraInterim/", tstep)
names(in_paths) <- names(era_paths) <- tstep

#Outputs will be saved to this directory
out_path_flx <- paste0(out_path, "/FLUXNET2015_sites/")

#Run only some sites, e.g. "AU-Ync" or c("AU-APL", "AU-Sil"); NA runs all sites
only_site <- c("AU-Tum", "AU-Otw")

#Anna's thresholds (Step1 lines 25-32)
missing_met  <- 100   #max. percent missing (must be set)
missing_flux <- 100
gapfill_met_tier1 <- 100  #max. gapfilled percentage
gapfill_met_tier2 <- 100
gapfill_flux      <- 100
min_yrs      <- 1   #min. number of consecutive years


##########################################
### Part 0: unzip and rename the files ###
##########################################

for (d in c(in_paths, era_paths)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

zips <- list.files(zip_path, pattern="^TERN_.*_FLUXNET_.*\\.zip$", full.names=TRUE)
if (length(zips) == 0) stop("No TERN_*_FLUXNET_*.zip files in ", zip_path)

#Site code from a file name, e.g. "AU-Ync"
code_of <- function(f) regmatches(basename(f), regexpr("[A-Z]{2}-[A-Za-z0-9]{3}", basename(f)))

#Old FLUXNET2015 name, e.g. FLX_AU-Ync_FLUXNET2015_FULLSET_HH_2012-2025_1-3.csv
old_name <- function(f, subset, res) {
  years <- regmatches(basename(f), regexpr("[0-9]{4}-[0-9]{4}", basename(f)))
  paste0("FLX_", code_of(f), "_FLUXNET2015_", subset, "_", res, "_", years, "_1-3.csv")
}

if (!all(is.na(only_site))) zips <- zips[sapply(zips, code_of) %in% only_site]
if (length(zips) == 0) stop("No zip for only_site ", paste(only_site, collapse=", "))

cat("Part 0: extracting FLUXMET and ERA5 (HH, or HR for hourly sites) from", length(zips), "zip(s)\n")
pb <- txtProgressBar(min=0, max=length(zips), style=3)
problems <- c()

for (z in seq_along(zips)) {

  inside <- unzip(zips[z], list=TRUE)$Name
  #CHANGED: half-hourly if the zip has HH files, otherwise hourly (HR)
  k <- if (any(grepl("_FLUXMET_HH_.*\\.csv$", inside))) "Halfhourly" else "Hourly"
  flux_f <- inside[grepl(paste0("_FLUXMET_", res_code[k], "_.*\\.csv$"), inside)]
  era_f  <- inside[grepl(paste0("_ERA5_", res_code[k], "_.*\\.csv$"), inside)]

  if (length(flux_f) != 1 || length(era_f) != 1) {
    problems <- c(problems, paste0(basename(zips[z]), ": found ", length(flux_f),
                                   " FLUXMET and ", length(era_f), " ERA5 files"))
    setTxtProgressBar(pb, z); next
  }

  #Extract (skipped if already done) and rename
  for (x in list(c(flux_f, in_paths[k], "FULLSET"), c(era_f, era_paths[k], "ERAI"))) {
    target <- file.path(x[2], old_name(x[1], x[3], res_code[k]))
    if (!file.exists(target)) {
      tmp <- unzip(zips[z], files=x[1], exdir=tempdir(), junkpaths=TRUE)
      file.copy(tmp, target); unlink(tmp)
    }
  }
  setTxtProgressBar(pb, z)
}
close(pb)
if (length(problems) > 0) { cat("\nSkipped zips:\n"); print(problems) }


################################################
### Part 1: Anna's FLUXNET2015 section       ###
################################################

#Remove old outputs (only this site's if testing one site)
if (all(is.na(only_site))) {
  unlink(out_path_flx, recursive = TRUE)
} else {
  for (o in only_site) unlink(list.files(out_path_flx, pattern=o, full.names=TRUE, recursive=TRUE))
}
dir.create(out_path_flx, recursive=TRUE, showWarnings=FALSE)

# Input Fluxnet data files (using FULLSET in this example, se R/Helpers.R for details)
#CHANGED: files from both the Hourly and Halfhourly folders (Anna loops over tstep)
infiles <- unlist(lapply(in_paths, get_fluxnet_files))
if (!all(is.na(only_site))) infiles <- infiles[sapply(infiles, get_path_site_code) %in% only_site]

#Retrieve site codes
site_codes <- sapply(infiles, get_path_site_code)

#Retrieve dataset versions
datasetversions <- sapply(infiles, get_fluxnet_version_no)

# Find ERA-files corresponding to site codes
#CHANGED: ERA file from the same time-step folder as the site's flux file
infile_tstep  <- ifelse(grepl("/Hourly$", dirname(infiles)), "Hourly", "Halfhourly")
ERA_files     <- lapply(seq_along(site_codes), function(i)
                   get_fluxnet_erai_files(era_paths[infile_tstep[i]], site_code=site_codes[i]))

#Stop if didn't find ERA files
if(any(sapply(ERA_files, length)==0)){
  stop("No ERA files found, amend input path")
}

#Sites FluxnetLSM has no metadata for (it then tries the web; may fail)
meta <- read.csv(system.file("extdata", "Site_metadata.csv", package="FluxnetLSM"))
if (any(!(site_codes %in% meta$SiteCode))) {
  cat("\nNot in installed Site_metadata.csv:",
      paste(site_codes[!(site_codes %in% meta$SiteCode)], collapse=", "), "\n")
}


### Process files ###

cat("\nPart 1: converting", length(infiles), "site(s) to PLUMBER2 NetCDF\n")
results <- data.frame(site=site_codes, status=NA, minutes=NA, message=NA, row.names=NULL)
t0 <- Sys.time()

for (s in seq_along(infiles)) {

  elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  eta <- if (s > 1) round(elapsed / (s - 1) * (length(infiles) - s + 1)) else "?"
  cat(sprintf("\n[%d/%d, %.0f%%, %.0f min so far, ~%s min left] %s ... ",
              s, length(infiles), 100 * (s - 1) / length(infiles), elapsed, eta, site_codes[s]))
  t_site <- Sys.time()

  results$status[s] <- tryCatch({
    convert_fluxnet_to_netcdf(infile=infiles[s], site_code=site_codes[s], out_path=out_path_flx,
                              datasetversion=datasetversions[s], met_gapfill="ERAinterim",
                              flux_gapfill="statistical", era_file=ERA_files[[s]],
                              missing_met=missing_met, missing_flux=missing_flux,
                              gapfill_met_tier1=gapfill_met_tier1,
                              gapfill_met_tier2=gapfill_met_tier2,
                              gapfill_flux=gapfill_flux, min_yrs=min_yrs,
                              check_range_action="warn",
                              include_all_eval=TRUE)
    "ok"
  }, error = function(e) { results$message[s] <<- conditionMessage(e); "FAILED" })

  results$minutes[s] <- round(as.numeric(difftime(Sys.time(), t_site, units="mins")), 1)

  #Check output files exist
  if (results$status[s] == "ok" &&
      length(list.files(paste0(out_path_flx, "/Nc_files/Met"), pattern=site_codes[s])) == 0) {
    results$status[s] <- "NO OUTPUT"; results$message[s] <- "see Logs folder"
  }
  cat(results$status[s], "(", results$minutes[s], "min )",
      ifelse(is.na(results$message[s]), "", results$message[s]))
}


###############
### Summary ###
###############

cat("\n\n======== SUMMARY ========\n")
print(results, right=FALSE)
write.csv(results, paste0(out_path, "/step1_fluxnet_summary.csv"), row.names=FALSE)
cat("\nOutputs:", paste0(out_path_flx, "Nc_files/{Met,Flux}"),
    "\nLogs:   ", paste0(out_path_flx, "Logs"),
    "\nSummary:", paste0(out_path, "/step1_fluxnet_summary.csv"), "\n")
