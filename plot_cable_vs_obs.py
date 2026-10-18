#!/usr/bin/env python3
"""Plot observed vs CABLE Rnet, Qle, Qh, GPP and NEE at every site.

For each site in cable_results/AU-xxx/AU-xxx_out.nc this finds the matching
processed/FLUXNET2015_sites/Nc_files/Flux/AU-xxx_*_Flux.nc, aligns the two on
time, averages to daily means and draws one PNG per site with five stacked
panels (r, RMSE and bias in each panel).

Outputs (in --out-dir, default processed/FLUXNET2015_sites/Figures/CABLE_vs_obs):
    AU-xxx_cable_vs_obs.png   one figure per site
    stats.csv                 n days, r, RMSE, bias per site and variable

Usage (from the repo root):
    python plot_cable_vs_obs.py
    python plot_cable_vs_obs.py --measured-only   # drop poorly gap-filled obs (qc > 1)
    python plot_cable_vs_obs.py --sites AU-Cum AU-How

Needs numpy, netCDF4 and matplotlib.
"""
import argparse
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
CABLE_DIR = os.path.join(ROOT, "cable_results")
FLUX_DIR = os.path.join(ROOT, "processed", "FLUXNET2015_sites", "Nc_files", "Flux")
OUT_DIR = os.path.join(ROOT, "processed", "FLUXNET2015_sites", "Figures", "CABLE_vs_obs")
MIN_DAY_FRACTION = 0.5  # a day needs at least this fraction of valid obs steps


def read_series(path, var, qc_max=None):
    """Return (datetimes, values) with fill values and (optionally) poor qc as NaN."""
    with netCDF4.Dataset(path) as ds:
        t = ds["time"]
        times = netCDF4.num2date(t[:], t.units, getattr(t, "calendar", "standard"),
                                 only_use_cftime_datetimes=False)
        x = np.ma.filled(ds[var][:].astype("f8"), np.nan).reshape(len(times), -1)[:, 0]
        x[(x <= -9999) | (np.abs(x) > 1e30)] = np.nan
        if qc_max is not None and f"{var}_qc" in ds.variables:
            qc = np.ma.filled(ds[f"{var}_qc"][:].astype("f8"), np.nan).ravel()
            x[~(qc <= qc_max)] = np.nan
    return np.array(times, dtype="datetime64[s]"), x


def daily_mean(times, x):
    """Daily means; days with too few valid steps are NaN."""
    days = times.astype("datetime64[D]")
    uniq, idx = np.unique(days, return_inverse=True)
    valid = ~np.isnan(x)
    n_all = np.bincount(idx)
    n_ok = np.bincount(idx, weights=valid)
    total = np.bincount(idx, weights=np.where(valid, x, 0.0))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / n_ok
    mean[n_ok < MIN_DAY_FRACTION * n_all] = np.nan
    return uniq, mean


def stats(obs, mod):
    ok = ~np.isnan(obs) & ~np.isnan(mod)
    if ok.sum() < 3:
        return dict(n_days=int(ok.sum()), r=np.nan, rmse=np.nan, bias=np.nan)
    o, m = obs[ok], mod[ok]
    return dict(n_days=int(ok.sum()), r=float(np.corrcoef(o, m)[0, 1]),
                rmse=float(np.sqrt(np.mean((m - o) ** 2))), bias=float(np.mean(m - o)))


VARS = [  # (name, label, units)
    ("Rnet", "Net radiation", "W m$^{-2}$"),
    ("Qle", "Latent heat", "W m$^{-2}$"),
    ("Qh", "Sensible heat", "W m$^{-2}$"),
    ("GPP", "GPP", "µmol m$^{-2}$ s$^{-1}$"),
    ("NEE", "NEE", "µmol m$^{-2}$ s$^{-1}$"),
]


def site_data(site, qc_max):
    """Daily obs and CABLE series for every variable in VARS: {var: (days, obs, mod)}."""
    cable_file = os.path.join(CABLE_DIR, site, f"{site}_out.nc")
    flux_files = glob.glob(os.path.join(FLUX_DIR, f"{site}_*_Flux.nc"))
    if not os.path.exists(cable_file) or not flux_files:
        return None
    out = {}
    for var, _, _ in VARS:
        try:
            t_mod, mod = read_series(cable_file, var)
            t_obs, obs = read_series(flux_files[0], var, qc_max)
        except (KeyError, IndexError):
            continue  # variable missing from one of the files
        # align on common timestamps (both files are in local standard time)
        common, i_obs, i_mod = np.intersect1d(t_obs, t_mod, return_indices=True)
        if common.size == 0:
            continue
        days, obs_d = daily_mean(common, obs[i_obs])
        _, mod_d = daily_mean(common, mod[i_mod])
        out[var] = (days, obs_d, mod_d)
    return out or None


def plot_site(site, data):
    fig, axes = plt.subplots(len(VARS), 1, figsize=(13, 11), sharex=True)
    for ax, (var, label, units) in zip(axes, VARS):
        ax.set_ylabel(f"{label}\n({units})", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.25)
        if var not in data:
            ax.text(0.5, 0.5, "not available", transform=ax.transAxes, ha="center", color="grey")
            continue
        days, obs, mod = data[var]
        d = days.astype(object)
        ax.plot(d, obs, color="k", lw=0.6, label="Observed")
        ax.plot(d, mod, color="tab:red", lw=0.6, alpha=0.75, label="CABLE")
        s = stats(obs, mod)
        ax.text(0.005, 0.97, f"r = {s['r']:.2f}   RMSE = {s['rmse']:.2f}   bias = {s['bias']:+.2f}",
                transform=ax.transAxes, va="top", fontsize=8.5,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))
    axes[0].legend(loc="upper right", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"{site}: CABLE vs observed (daily means)", fontsize=13, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--sites", nargs="*", help="site codes (default: every folder in cable_results)")
    p.add_argument("--measured-only", action="store_true",
                   help="use only measured and good-quality gap-filled obs (qc <= 1, where a qc flag exists)")
    p.add_argument("--out-dir", default=OUT_DIR)
    a = p.parse_args()

    sites = a.sites or sorted(os.path.basename(d) for d in glob.glob(os.path.join(CABLE_DIR, "AU-*"))
                              if os.path.isdir(d))
    qc_max = 1 if a.measured_only else None
    os.makedirs(a.out_dir, exist_ok=True)

    rows = []
    for site in sites:
        data = site_data(site, qc_max)
        if data is None:
            print(f"{site}: skipped (missing CABLE output or Flux file, or no overlapping times)")
            continue
        fig = plot_site(site, data)
        fig.savefig(os.path.join(a.out_dir, f"{site}_cable_vs_obs.png"), dpi=130)
        plt.close(fig)
        for var, (_, obs, mod) in data.items():
            s = stats(obs, mod)
            rows.append(dict(site=site, variable=var,
                             **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()}))
        print(f"{site}: done")

    with open(os.path.join(a.out_dir, "stats.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["site", "variable", "n_days", "r", "rmse", "bias"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nFigures and stats.csv written to {a.out_dir}")


if __name__ == "__main__":
    main()
