# get_modis_lai.R
#
# Downloads MODIS LAI (MCD15A2H: Lai_500m, LaiStdDev_500m, FparLai_QC) for each
# site: same product, bands and 1 km box as Anna's Get_MODIS_LAI_for_all_sites.R.
#
# Difference: instead of one huge batch request per band (which NASA's server
# often answers with "502" errors), it downloads one YEAR at a time, saves each
# year to its own small file, and retries failed years. If it stops, or you stop
# it, just run it again: years already downloaded are skipped.
#
# At the end it joins the years into one CSV per site and band, named the way
# add_lai_to_met.R (Anna's Step2) expects, e.g. AU-Cum_MCD15A2H_Lai_500m_2002-2026.csv

library(MODISTools)
library(ncdf4)

# Run several copies at once: in each window set a different `part`, e.g.
#   Rscript get_modis_lai.R 1 4     (window 1 of 4; then 2 4, 3 4, 4 4)
# or in an R console:  part <- 1; n_parts <- 4; source("get_modis_lai.R", echo = TRUE)
# Each copy downloads every n_parts-th site. When all are done, run the script
# once more with no part: everything is already downloaded, so it just joins
# the files and checks them.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) { part <- as.integer(args[1]); n_parts <- as.integer(args[2]) }
rm(list = setdiff(ls(all = TRUE), c("part", "n_parts")))

path   <- "C:/Users/jamie/extending-plumber2-flux-sites/processed/"
outdir <- paste0(path, "/MODIS_LAI_time_series/Raw_data/")
chunks <- paste0(path, "/MODIS_LAI_time_series/yearly_chunks/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(chunks, recursive = TRUE, showWarnings = FALSE)

product   <- "MCD15A2H"
bands     <- c("Lai_500m", "LaiStdDev_500m", "FparLai_QC")
km        <- 1
# MCD15A2H (Terra+Aqua combined) only starts in July 2002, so earlier years
# have no data. Step2 fills years before MODIS with the mean seasonal cycle.
years     <- 2002:as.numeric(format(Sys.Date(), "%Y"))
max_tries <- 5


### Site codes and coordinates from the Step1 Met files ###

site_files <- list.files(paste0(path, "/FLUXNET2015_sites/Nc_files/Met/"), full.names = TRUE)
if (length(site_files) == 0) stop("No Met files found; run run_step1.R first")

sites <- do.call(rbind, lapply(site_files, function(f) {
  nc <- nc_open(f)
  out <- data.frame(site = ncatt_get(nc, varid = 0, "site_code")$value,
                    lat  = as.numeric(ncvar_get(nc, "latitude")),
                    lon  = as.numeric(ncvar_get(nc, "longitude")))
  nc_close(nc)
  out
}))
sites <- sites[!duplicated(sites$site), ]
if (exists("part") && exists("n_parts")) {
  sites <- sites[seq_len(nrow(sites)) %% n_parts == part %% n_parts, ]
  cat("Part", part, "of", n_parts, ":", nrow(sites), "sites\n")
}
print(sites)


### Download, one site / band / year at a time ###

n_total <- nrow(sites) * length(bands) * length(years)
n_done  <- 0
n_new   <- 0
t_start <- Sys.time()

for (s in seq_len(nrow(sites))) {
  for (band in bands) {
    for (y in years) {

      n_done <- n_done + 1
      chunk_file <- paste0(chunks, sites$site[s], "_", band, "_", y, ".csv")
      if (file.exists(chunk_file)) next   # already downloaded

      # Progress: request number, % done, and a time estimate from the
      # average time per download so far in this run
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      eta <- if (n_new > 0) round(elapsed / n_new * (n_total - n_done + 1)) else NA
      cat(sprintf("[%d/%d, %.0f%%, %.0f min so far, ~%s min left] %s %s %d ... ",
                  n_done, n_total, 100 * n_done / n_total, elapsed,
                  ifelse(is.na(eta), "?", eta), sites$site[s], band, y))
      x <- NULL
      for (try in 1:max_tries) {
        err <- NULL
        x <- tryCatch(
          mt_subset(product = product, band = band,
                    lat = sites$lat[s], lon = sites$lon[s],
                    start = paste0(y, "-01-01"), end = paste0(y, "-12-31"),
                    km_lr = km, km_ab = km, site_name = sites$site[s],
                    internal = TRUE, progress = FALSE),
          error = function(e) { err <<- conditionMessage(e); NULL })
        # Check it is a full download: 25 pixels for every date, no blanks
        if (!is.null(x) && nrow(x) > 0) {
          n_dates <- length(unique(x$calendar_date))
          if (nrow(x) == 25 * n_dates && !any(is.na(x$value))) break
          err <- sprintf("incomplete: %d rows for %d dates", nrow(x), n_dates)
        }
        x <- NULL
        cat("[", err, "] retry", try, "... ")
        Sys.sleep(10 * try)   # wait a bit longer each time
      }

      if (!is.null(x)) {
        write.csv(x, chunk_file, row.names = FALSE)
        n_new <- n_new + 1
        cat("ok,", length(unique(x$calendar_date)), "dates\n")
      } else {
        cat("FAILED (run the script again later)\n")
      }
    }
  }
}


### Join the years into one file per site and band ###

missing <- c()
for (s in seq_len(nrow(sites))) {
  for (band in bands) {
    files <- paste0(chunks, sites$site[s], "_", band, "_", years, ".csv")
    have  <- file.exists(files)
    if (!all(have)) {
      missing <- c(missing, paste(sites$site[s], band, years[!have]))
      next
    }
    all <- do.call(rbind, lapply(files, read.csv))
    # Step2 takes each pixel's rows in order, so sort by pixel then date
    all <- all[order(all$pixel, as.Date(all$calendar_date)), ]
    all <- all[!duplicated(all[, c("pixel", "calendar_date")]), ]
    out <- paste0(outdir, sites$site[s], "_", product, "_", band, "_",
                  min(years), "-", max(years), ".csv")
    write.csv(all, out, row.names = FALSE)
    cat("Wrote", basename(out), "\n")
  }
}

### Final check: every site has all 3 bands, same dates, 25 pixels ###

cat("\n==== Check of Raw_data files ====\n")
ok_all <- length(missing) == 0
for (site in sites$site) {
  f <- list.files(outdir, pattern = paste0("^", site, "_", product, "_"), full.names = TRUE)
  if (length(f) != 3) { cat(site, ": expected 3 files, found", length(f), "\n"); ok_all <- FALSE; next }
  d <- lapply(f, read.csv)
  dates <- lapply(d, function(x) sort(unique(x$calendar_date)))
  same  <- all(sapply(dates, identical, dates[[1]]))
  px    <- all(sapply(d, function(x) nrow(x) == 25 * length(unique(x$calendar_date))))
  cat(sprintf("%s: %d dates, %s to %s, same dates in all bands: %s, 25 pixels each: %s\n",
              site, length(dates[[1]]), dates[[1]][1], tail(dates[[1]], 1), same, px))
  if (!same || !px) ok_all <- FALSE
}

if (length(missing) > 0) {
  cat("\nStill missing, run the script again:\n")
  print(missing)
} else if (!ok_all) {
  cat("\nSomething above is FALSE: send me the output before going on.\n")
} else {
  cat("\nAll good. Next: source('add_lai_to_met.R', echo = TRUE)\n")
}
