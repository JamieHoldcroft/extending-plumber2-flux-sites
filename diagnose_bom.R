# diagnose_bom.R: runs the LAI steps of add_bom_lai_to_met.R (same code) for EVERY
# site in Nc_files/Met and saves one figure per site showing the series after each
# step, plus the climatology. Prints a summary table at the end.
# The BOM version of diagnose.R (which does the same for MODIS).
# Figures: processed/BOM_LAI_diagnose/<site>_BOM_LAI_steps.png

library(ncdf4)
library(zoo)
rm(list = ls(all = TRUE))

path     <- "C:/Users/jamie/extending-plumber2-flux-sites/processed/"
met_dir  <- paste0(path, "/FLUXNET2015_sites/Nc_files/Met/")
lai_path <- "C:/Users/jamie/extending-plumber2-flux-sites/BOM_LAI"
outdir   <- paste0(path, "/BOM_LAI_diagnose/")
dir.create(outdir, showWarnings = FALSE)

# Set to "2003-01-01" to drop the years before it; NA keeps everything (2000-2024)
start_date <- NA

lai_files <- list.files(lai_path, pattern = "lai_bom_.*_Jul2025.nc", full.names = TRUE)
if (length(lai_files) != 25) stop("incorrect LAI files found")

# Site codes and coordinates from the met files
met_files <- list.files(met_dir, pattern = "_Met.nc", full.names = TRUE)
sites <- lat <- lon <- c()
for (f in met_files) {
  nc <- nc_open(f)
  sites <- c(sites, ncatt_get(nc, varid = 0, "site_code")$value)
  lat <- c(lat, as.numeric(ncvar_get(nc, "latitude")))
  lon <- c(lon, as.numeric(ncvar_get(nc, "longitude")))
  nc_close(nc)
}

# BOM grid and the pixel each site falls in (same as add_bom_lai_to_met.R)
nc_lai   <- nc_open(lai_files[1])
lai_lats <- ncvar_get(nc_lai, "latitude")
lai_lons <- ncvar_get(nc_lai, "longitude")
nc_close(nc_lai)
site_i <- sapply(lon, function(x) which.min(abs(lai_lons - x)))
site_j <- sapply(lat, function(x) which.min(abs(lai_lats - x)))

# Read every site's pixel, one year at a time (each file opened once)
cat("Reading BOM files\n")
raw <- lapply(sites, function(x) vector())
used_px <- rep("own", length(sites))
pb <- txtProgressBar(min = 0, max = length(lai_files), style = 3)
for (y in 1:length(lai_files)) {
  nc_lai <- nc_open(lai_files[y])
  for (s in 1:length(sites)) {
    vals <- ncvar_get(nc_lai, "lai", start = c(site_i[s], site_j[s], 1), count = c(1, 1, -1))
    # Empty pixel: nearest non-empty pixel within 3 degrees (same as add_bom_lai_to_met.R)
    if (all(is.na(vals))) {
      i_win <- which(abs(lai_lons - lon[s]) <= 3)
      j_win <- which(abs(lai_lats - lat[s]) <= 3)
      block <- ncvar_get(nc_lai, "lai", start = c(min(i_win), min(j_win), 1),
                         count = c(length(i_win), length(j_win), -1))
      has_data <- which(apply(block, 3, function(m) any(!is.na(m))))[1]
      ok   <- which(!is.na(block[, , has_data]), arr.ind = TRUE)
      dist <- (lai_lons[i_win][ok[, 1]] - lon[s])^2 + (lai_lats[j_win][ok[, 2]] - lat[s])^2
      best <- ok[which.min(dist), ]
      vals <- block[best[1], best[2], ]
      used_px[s] <- "nearest"
    }
    raw[[s]] <- append(raw[[s]], vals)
  }
  nc_close(nc_lai)
  setTxtProgressBar(pb, y)
}
close(pb)

all_time <- seq.Date(as.Date("2000-01-01"), by = "month", length.out = length(raw[[1]]))
summary <- data.frame()

for (s in 1:length(sites)) {

  site <- sites[s]
  cat(sprintf("[%d/%d] %s\n", s, length(sites), site))

  keep <- if (is.na(start_date)) rep(TRUE, length(all_time)) else all_time >= as.Date(start_date)
  time   <- all_time[keep]
  lai_ts <- raw[[s]][keep]
  n <- length(lai_ts)

  # Gaps
  r <- rle(is.na(lai_ts))
  lead_na <- if (r$values[1]) r$lengths[1] else 0

  # Gap-fill with the calendar-month mean (same as add_bom_lai_to_met.R)
  month_of_step <- as.numeric(format(time, "%m"))
  filled <- lai_ts
  for (m in which(is.na(lai_ts))) {
    filled[m] <- mean(lai_ts[month_of_step == month_of_step[m]], na.rm = TRUE)
  }

  # Smoothing (same call as add_bom_lai_to_met.R / Anna's Step3)
  spl_raw <- smooth.spline(time, filled)$y
  smo <- pmax(spl_raw, 0)
  lin <- na.approx(lai_ts, rule = 2)
  clim <- tapply(smo, month_of_step, mean)

  summary <- rbind(summary, data.frame(
    site = site, pixel = used_px[s], first_date = time[1], missing_at_start = lead_na,
    missing_total = sum(is.na(lai_ts)), longest_gap = max(c(0, r$lengths[r$values])),
    raw_min = round(min(lai_ts, na.rm = TRUE), 2), raw_max = round(max(lai_ts, na.rm = TRUE), 2),
    spline_min = round(min(spl_raw), 2), spline_max = round(max(spl_raw), 2),
    clim_min = round(min(clim), 2), clim_max = round(max(clim), 2)))

  out <- paste0(outdir, site, "_BOM_LAI_steps.png")
  png(out, width = 1600, height = 1300, res = 110)
  layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(1, 1, 1, 1))
  par(mar = c(3, 4, 2, 1))
  plot(time, lai_ts, pch = 20, cex = 0.6, ylab = "LAI",
       main = paste(site, "- 1. BOM pixel at the site (raw, monthly)"))
  plot(time, filled, type = "l", col = "red", ylab = "LAI",
       main = "2. After calendar-month-mean gap-fill (red) vs raw (black)")
  points(time, lai_ts, pch = 20, cex = 0.4); abline(h = 0, lty = 2)
  plot(time, smo, type = "l", col = "blue", ylab = "LAI",
       main = "3. After smoothing (blue) vs linear gap-fill alternative (green)")
  lines(time, lin, col = "darkgreen"); points(time, lai_ts, pch = 20, cex = 0.4)
  plot(as.numeric(names(clim)), clim, type = "l", xlab = "Month", ylab = "LAI",
       main = "4. Climatology (mean of smoothed LAI per month)")
  dev.off()
}

cat("\n==== SUMMARY ====\n")
print(summary, row.names = FALSE)
cat("\npixel = nearest means the site's own pixel was empty (sea/masked) and the nearest non-empty one was used.\n")
cat("A site has the AU-Sno problem if missing_at_start > 0 and spline_max is far above raw_max.\n")
cat("Figures in", outdir, "\n")
