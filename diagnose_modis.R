# diagnose.R: runs the LAI steps of Anna's Step2 (same code) for EVERY site in
# Raw_data and saves one figure per site showing the series after each step,
# plus the climatology. Prints a summary table at the end.
# Figures: processed/PLUMBER2_site_selection_figs/<site>_LAI_steps.png

library(zoo)
rm(list = ls(all = TRUE))

path   <- "C:/Users/jamie/extending-plumber2-flux-sites/processed/"
d      <- paste0(path, "/MODIS_LAI_time_series/Raw_data/")
outdir <- paste0(path, "/MODIS_LAI_diagnose/")
dir.create(outdir, showWarnings = FALSE)

# Set to "2003-01-01" to drop 2002 (the AU-Sno test); NA keeps everything
start_date <- NA

# Optional: drop the FIRST YEAR of data for the sites listed here (same as in
# add_lai_to_met.R). Uncomment this line and the matching line in read() below.
drop_first_year_sites <- c("AU-Sno", "AU-APL", "AU-ASM")

pixel_no <- c(7, 8, 9, 12, 13, 14, 17, 18, 19)
qc_flags <- c(0, 2, 24, 26, 32, 34, 56, 58)

sites <- unique(sub("_MCD15A2H.*", "", list.files(d, pattern = "_MCD15A2H_Lai_500m_")))
summary <- data.frame()

for (site in sites) {

  cat(sprintf("[%d/%d] %s\n", match(site, sites), length(sites), site))

  read <- function(band) {
    x <- read.csv(list.files(d, pattern = paste0(site, "_MCD15A2H_", band), full.names = TRUE))
    if (!is.na(start_date)) x <- x[x$calendar_date >= start_date, ]
    if (site %in% drop_first_year_sites)
      x <- x[x$calendar_date >= paste0(min(as.numeric(substr(x$calendar_date, 1, 4))) + 1, "-01-01"), ]
    x
  }
  lai <- read("Lai_500m_"); qc <- read("FparLai_QC"); sd <- read("LaiStdDev_500m_")

  n <- min(nrow(lai), nrow(sd), nrow(qc)) / max(lai$pixel)
  time <- as.Date(lai$calendar_date[lai$pixel == pixel_no[1]])[1:n]
  get <- function(df, scale) sapply(pixel_no, function(p) df$value[df$pixel == p][1:n] * scale)
  lai_px <- get(lai, lai$scale[1]); sd_px <- get(sd, sd$scale[1]); qc_px <- get(qc, 1)

  # Same masking as Step2
  lai_px[!(qc_px %in% qc_flags)] <- NA; sd_px[!(qc_px %in% qc_flags)] <- NA
  sd_px[sd_px < 0.1] <- NA; lai_px[is.na(sd_px)] <- NA
  lai_px[lai_px > 10] <- NA

  # SD-weighted mean of the 9 pixels per date
  lai_ts <- sapply(1:n, function(t) {
    if (all(is.na(lai_px[t, ]))) return(NA)
    w <- (1 / sd_px[t, ]^2) / sum(1 / sd_px[t, ]^2, na.rm = TRUE)
    weighted.mean(lai_px[t, ], w = w, na.rm = TRUE)
  })

  # Gaps
  r <- rle(is.na(lai_ts))
  lead_na <- if (r$values[1]) r$lengths[1] else 0

  # Spline gap-fill and smoothing (same calls as Step2)
  x <- 1:n
  spl_raw <- splinefun(x = x, y = lai_ts, method = "fmm", ties = mean)(x)
  spl <- pmax(spl_raw, 0)
  smo <- pmax(smooth.spline(x, spl)$y, 0)
  lin <- na.approx(lai_ts, rule = 2)
  clim <- tapply(smo, format(time, "%j"), mean)

  summary <- rbind(summary, data.frame(
    site = site, first_date = time[1], missing_at_start = lead_na,
    missing_total = sum(is.na(lai_ts)), longest_gap = max(c(0, r$lengths[r$values])),
    raw_min = round(min(lai_ts, na.rm = TRUE), 2), raw_max = round(max(lai_ts, na.rm = TRUE), 2),
    spline_min = round(min(spl_raw), 2), spline_max = round(max(spl_raw), 2),
    clim_min = round(min(clim), 2), clim_max = round(max(clim), 2)))

  out <- paste0(outdir, site, "_LAI_steps.png")
  png(out, width = 1600, height = 1300, res = 110)
  layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(1, 1, 1, 1))
  par(mar = c(3, 4, 2, 1))
  plot(time, lai_ts, pch = 20, cex = 0.6, ylab = "LAI",
       main = paste(site, "- 1. Weighted mean of good pixels (raw)"))
  plot(time, spl_raw, type = "l", col = "red", ylab = "LAI",
       main = "2. After Step2 spline gap-fill (red) vs raw (black)")
  points(time, lai_ts, pch = 20, cex = 0.4); abline(h = 0, lty = 2)
  plot(time, smo, type = "l", col = "blue", ylab = "LAI",
       main = "3. After smoothing (blue) vs linear gap-fill alternative (green)")
  lines(time, lin, col = "darkgreen"); points(time, lai_ts, pch = 20, cex = 0.4)
  plot(as.numeric(names(clim)), clim, type = "l", xlab = "Day of year", ylab = "LAI",
       main = "4. Climatology (mean of smoothed LAI per 8-day date)")
  dev.off()
}

cat("\n==== SUMMARY ====\n")
print(summary, row.names = FALSE)
cat("\nA site has the AU-Sno problem if missing_at_start > 0 and spline_max is far above raw_max.\n")
cat("Figures in", outdir, "\n")
