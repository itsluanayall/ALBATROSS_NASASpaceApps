"""
Visualizes Ncdf4 file with the xarray library and generates 
a geographic map with the cartopy library.
"""

import xarray as xr
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy
import os

#Load data in xarray DataSet.
homepath = os.getenv("HOME")
path = os.path.join(homepath, "ALBATROSS_NasaSpaceAppsChallenge/oscar_vel2014.nc.gz.nc4")
ds = xr.open_dataset(path)

#Variables.
meridional_velocity = ds.v
zonal_velocity = ds.u
sum_vector = numpy.sqrt(ds.u**2 + ds.v**2) #Norm of the vector sum of v and u.
sum_vector.attrs = {'units': 'm/s', 'long_name': 'Ocean Surface Current Speed'}

#Plotting the maps.
for t in ds.time:
    speed = sum_vector.sel(time=t, method='nearest')
    fig = plt.figure(figsize=(40,25))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND, zorder=100, color='k')
    speed.plot(ax=ax, transform=ccrs.PlateCarree(), vmin=0, vmax=1, cbar_kwargs={'shrink': 0.7})

    name_fig = str(t) + '.jpeg'
    plt.tight_layout()
    plt.show()
