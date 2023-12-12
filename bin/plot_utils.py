#!/usr/bin/env python


import os
import sys
import fitsio
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator
from desimodel.footprint import radec2pix
from desiutil.plots import init_sky, plot_healpix_map
from astropy import units
from astropy.coordinates import SkyCoord
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
from astropy import units
from astropy.coordinates import SkyCoord


nest = True


def get_desfoot_fn():
    return os.path.join(
        os.getenv("DESI_ROOT"), "survey", "observations", "misc", "des_footprint.txt"
    )


def get_desi14kfoot_fn():
    return os.path.join(
        os.getenv("DESI_ROOT"),
        "users",
        "raichoor",
        "desi-14k-footprint",
        "desi-14k-footprint-dark.ecsv",
    )


def get_survprogs(prod):
    fn = os.path.join(
        os.getenv("DESI_ROOT"), "spectro", "redux", prod, f"tiles-{prod}.csv"
    )
    d = Table.read(fn)
    survprogs = {}
    for survey in np.unique(d["SURVEY"]):
        progs = np.unique(d["PROGRAM"][d["SURVEY"] == survey])
        survprogs[survey] = progs
    return survprogs


# This is an array wrapper object just so that we can set .vmin and .vmax
# which plot_healpix_map uses for color scaling.
# see https://numpy.org/doc/stable/user/basics.subclassing.html#simple-example-adding-an-extra-attribute-to-ndarray
# and https://stackoverflow.com/questions/67509913/add-an-attribute-to-a-numpy-array-in-runtime
class ArrayVMinMax(np.ndarray):
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        obj.vmin = 0.0
        obj.vmax = np.max(input_array)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.vmin = getattr(obj, "vmin", None)
        self.vmax = getattr(obj, "vmax", None)


def init_hpd_table(nside, nest=True):
    hpd = Table()
    hpd.meta["HPXNSIDE"], hpd.meta["HPXNEST"] = nside, nest
    npix = hp.nside2npix(nside)
    hpd["hpxpixel"] = np.arange(npix, dtype=int)

    thetas, phis = hp.pix2ang(nside, hpd["hpxpixel"], nest=nest)
    hpd["ra"], hpd["dec"] = np.degrees(phis), 90.0 - np.degrees(thetas)

    cs = SkyCoord(hpd["ra"] * units.degree, hpd["dec"] * units.degree, frame="icrs")
    hpd["galb"], hpd["gall"] = cs.galactic.b.value, cs.galactic.l.value

    return hpd


def get_pix_densities(prod, survey, progs, nside, debug=False):

    npix = hp.nside2npix(nside)
    pixarea = hp.nside2pixarea(nside, degrees=True)

    # read all requested programs
    zcatdir = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux", prod, "zcatalog")
    densities = []
    for prog in progs:
        fn = os.path.join(zcatdir, f"zpix-{survey}-{prog}.fits")
        hdr = fitsio.read_header(fn, "ZCATALOG")
        n = hdr["NAXIS2"]
        if debug:
            np.random.seed(1234)
            rows = np.random.choice(n, size=int(0.01 * n), replace=False)
        else:
            rows = np.arange(n, dtype=int)
        zcat = fitsio.read(
            fn,
            "ZCATALOG",
            rows=rows,
            columns=("TARGETID", "TARGET_RA", "TARGET_DEC", "ZWARN", "OBJTYPE"),
        )
        # - Basic quality cuts
        ii = (zcat["ZWARN"] == 0) & (zcat["OBJTYPE"] == "TGT")
        zcat = zcat[ii]
        assert np.unique(zcat["TARGETID"]).size == len(zcat)

        density = np.zeros(npix)
        tgtpix = radec2pix(nside, zcat["TARGET_RA"], zcat["TARGET_DEC"])
        for h, c in zip(*np.unique(tgtpix, return_counts=True)):
            density[h] = c / pixarea
        densities.append(density)

    return densities


def get_prod_pix_densities(prod, nside, debug=False, verbose=True):

    # grab all {survey,program}
    survprogs = get_survprogs(prod)

    npix = hp.nside2npix(nside)
    t = init_hpd_table(nside, nest=True)

    if verbose:
        print("# SURVEY MEDIAN 95TH_PERC MAX")

    for survey in survprogs:
        # grab the density per {survey, prog}
        progs = survprogs[survey]
        densities = get_pix_densities(prod, survey, progs, nside, debug=debug)
        # add all programs
        alldensity = np.zeros(npix)
        for density in densities:
            alldensity += density
        # fill the Table()
        for prog, density in zip(progs, densities):
            key = f"{survey}_{prog}"
            t[key] = density
            if verbose:
                ii = t[key] > 0
                p50, p95, pmax = np.percentile(t[key][ii], (50, 95, 100))
                print(f"{key}\t{p50:5.0f}\t{p95:5.0f}\t{pmax:5.0f}")
        t[survey] = alldensity
        if verbose:
            key = survey
            ii = t[key] > 0
            if np.sum(ii) > 0:
                p50, p95, pmax = np.percentile(t[key][ii], (50, 95, 100))
                print(f"{key}\t{p50:5.0f}\t{p95:5.0f}\t{pmax:5.0f}")

    return t


# - rescale sequential color maps so that they don't start at white
def get_custom_cmaps():
    istart = 65
    cmap = plt.get_cmap("Greens")
    cmap_greens = cmap.from_list("Greens2", cmap(np.arange(istart, 256)))
    cmap = plt.get_cmap("Blues")
    cmap_blues = cmap.from_list("Blues2", cmap(np.arange(istart, 256)))
    cmap = plt.get_cmap("Oranges")
    cmap_oranges = cmap.from_list("Orange2", cmap(np.arange(istart, 256)))
    return cmap_greens, cmap_blues, cmap_oranges


def ar_plot_sky(ax, ras, decs, **kwargs):
    ii = ax.projection_ra(ras).argsort()
    _ = ax.plot(ax.projection_ra(ras[ii]), ax.projection_dec(decs[ii]), **kwargs)


def ar_init_sky(gp=True, des=True, nsbound=True, desi=True, ra_center=120):
    ax = init_sky(
        galactic_plane_color="none", ecliptic_plane_color="none", ra_center=ra_center
    )
    npt = 1000
    # gp
    if gp:
        cs = SkyCoord(
            l=np.linspace(0, 360, npt) * units.deg,
            b=np.zeros(npt) * units.deg,
            frame="galactic",
        )
        ras, decs = cs.icrs.ra.degree, cs.icrs.dec.degree
        ar_plot_sky(ax, ras, decs, c="k", lw=0.5, zorder=10)
    # des
    if des:
        fn = get_desfoot_fn()
        ras, decs = np.loadtxt(fn, unpack=True)
        _ = ax.plot(
            ax.projection_ra(ras), ax.projection_dec(decs), color="k", lw=0.5, zorder=10
        )
    # north-south boundary
    if nsbound:
        ras = np.arange(100, 280)
        decs = 32.375 + 0 * ras
        _ = ax.plot(
            ax.projection_ra(ras),
            ax.projection_dec(decs),
            color="k",
            lw=0.5,
            ls="--",
            zorder=10,
        )
    # desi
    if desi:
        fn = get_desi14kfoot_fn()
        d = Table.read(fn)
        for cap in ["NGC", "SGC"]:
            sel = d["CAP"] == cap
            _ = ax.plot(
                ax.projection_ra(d["RA"][sel]),
                ax.projection_dec(d["DEC"][sel]),
                color="k",
                lw=2,
                zorder=10,
            )
    #
    return ax


def ar_plot_sky_circles(ax, ra_center, dec_center, field_of_view, **kwargs):
    if isinstance(ra_center, int) | isinstance(ra_center, float):
        ra_center, dec_center = [ra_center], [dec_center]
    proj_edge = ax._ra_center - 180
    while proj_edge < 0:
        proj_edge += 360
    dRA = field_of_view / np.cos(np.radians(dec_center))
    edge_dist = np.fabs(np.fmod(ra_center - proj_edge, 360))
    wrapped = np.minimum(edge_dist, 360 - edge_dist) < 1.05 * 0.5 * dRA
    es = []
    for ra, dec, dra in zip(
        ax.projection_ra(ra_center[~wrapped]),
        ax.projection_dec(dec_center[~wrapped]),
        dRA[~wrapped],
    ):
        es.append(
            Ellipse(
                (ra, dec),
                np.radians(dra),
                np.radians(field_of_view),
            )
        )
    ax.add_collection(PatchCollection(es, **kwargs))


# adapted from https://github.com/desihub/desiutil/blob/5735fdc34c4e77c7fda84c92c32b9ac41158b8e1/py/desiutil/plots.py#L735-L857
def ar_sky_cbar(ax, sc, label, extend=None, mloc=None):
    cbar = plt.colorbar(
        sc,
        ax=ax,
        location="top",
        orientation="horizontal",
        spacing="proportional",
        extend=extend,
        extendfrac=0.025,
        pad=0.11,
        fraction=0.035,
        aspect=40,
    )
    cbar.ax.xaxis.set_ticks_position("bottom")
    cbar.set_label(label, labelpad=10)
    if mloc is not None:
        cbar.ax.xaxis.set_major_locator(MultipleLocator(mloc))


# climmloc = (clim[0], clim[1], mloc)
def plot_skymap(outpng, density, cmap=None, climmloc=None):

    # clim, mloc, extend
    if climmloc is None:
        clim = [np.nanmin(density), np.nanmax(density)]
        mloc = None
    else:
        clim = (climmloc[0], climmloc[1])
        mloc = climmloc[2]
    islo = density < clim[0]
    ishi = density > clim[1]
    extend = "neither"
    if (islo.sum() > 0) & (ishi.sum() == 0):
        extend = "min"
    if (islo.sum() == 0) & (ishi.sum() > 0):
        extend = "max"
    if (islo.sum() > 0) & (ishi.sum() > 0):
        extend = "both"
    density[density < clim[0]] = clim[0]
    density[density > clim[1]] = clim[1]

    fig = plt.figure(figsize=(6.0, 5.0), dpi=200)
    ax = ar_init_sky()
    plot_healpix_map(density, cmap=cmap, nest=nest, ax=ax, colorbar=False)

    # colorbar (dummy values under the desi edge contour)
    dummy_ras, dummy_decs, dummy_denss = np.array([240]), np.array([0]), np.array([0])
    sc = ax.scatter(
        ax.projection_ra(dummy_ras),
        ax.projection_dec(dummy_decs),
        c=dummy_denss,
        cmap=cmap,
        clim=clim,
        zorder=-1,
    )
    ar_sky_cbar(
        ax, sc, "Density of confident redshifts [deg$^{-2}$]", mloc=mloc, extend=extend
    )
    plt.savefig(outpng, bbox_inches="tight")
    plt.close()
