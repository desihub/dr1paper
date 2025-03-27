import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap
from astropy.cosmology import LambdaCDM
from astropy.coordinates import SkyCoord
import astropy.units as u
cosmo = LambdaCDM(H0=69.7, Om0=0.284, Ode0=0.702)

def trim_to_complete_region(catalog,ra_min=90, ra_max = 280, dec_min=-7, dec_max=7, z_min=0.0, z_max=.2):
    # Trim to complete region
    catalog = catalog[(catalog['RA']>ra_min) & (catalog['RA']<ra_max) & (catalog['DEC']>dec_min) & (catalog['DEC']<dec_max) &
                      (catalog['Z']>z_min) & (catalog['Z']<z_max)]
    return catalog

def get_comoving_coords(catalog):
    try:
        catalog['dist_comoving']
    except KeyError:
        # center so that the LOS is along the x axis
        catalog['RA_rot'] = catalog['RA'] - (np.max(catalog['RA']) + np.min(catalog['RA']))/2
        catalog['DEC_rot'] = catalog['DEC'] - (np.max(catalog['DEC']) + np.min(catalog['DEC']))/2
        # Convert to comoving coordinates but only if the catalog doesn't already have them
        catalog['dist_comoving'] = cosmo.comoving_distance(catalog['Z'])
    catalog['RA_rot'] -= 31 # just for plotting, remove later
    c = SkyCoord(ra=catalog['RA_rot'], dec=catalog['DEC_rot'], distance=catalog['dist_comoving'])
    catalog['x_comoving'] = c.cartesian.x
    catalog['y_comoving'] = c.cartesian.y
    catalog['z_comoving'] = c.cartesian.z
    return catalog

def prep_catalog(catalog, ra_min=150, ra_max=180, dec_min=-5, dec_max=5, z_min=0.0, z_max=1.6, sort='True'):
    try:
        catalog.keep_columns(['Z_not4clus', 'RA', 'DEC'])
        catalog.rename_column('Z_not4clus', 'Z')
    except KeyError:
        pass
    catalog = trim_to_complete_region(catalog, ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, z_min=z_min, z_max=z_max)
    catalog = get_comoving_coords(catalog)
    #catalog = commoving_trim(catalog)
    catalog.sort('z_comoving') # sort by z-comoving for depth when plotting
    return catalog

# setting up custom colormap "lacerta" from cmastro

cmap_file = 'lacerta.csv'

cmaps = dict()
kw = dict()
kw['delimiter'] = ','

rgb_data = np.loadtxt(cmap_file, **kw)
from pathlib import Path

cmap_name = f"cma:{Path(cmap_file).stem}"

cmaps[cmap_name] = ListedColormap(rgb_data, name=cmap_name)
cmaps[f"{cmap_name}_r"] = ListedColormap(rgb_data[::-1],
                                            name=f"{cmap_name}_r")

for cmap in cmaps.values():
    cm.register_cmap(cmap=cmap)