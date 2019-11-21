import os
import cmocean as cm

import cartopy as cpy
import cartopy.crs as ccrs

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

import xarray as xr
import regionmask

from matplotlib.backends.backend_pdf import PdfPages
from cams_visu import utils as u

generate_daily = False
ofig = './figures'
idir = '/local/AIX/tristan.harmel/project/ardyna/cams/data'

land_feat = cpy.feature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face', facecolor=cpy.feature.COLORS['land'])

crs = ccrs.NearsidePerspective(100, 71)
extent = [-180, 180, 65, 90]  # [90,180, 71, 78]

cmap = cm.tools.crop_by_percent(cm.cm.delta, 30, which='both')
cice = cm.tools.crop_by_percent(cm.cm.ice, 25, which='min')
caero = cm.tools.crop_by_percent(cmap, 20, which='max')
year = 2014

aspect_ratio = 1
years = np.arange(2004, 2018)
rows = len(years)
fig_clima, axs_clima = plt.subplots(nrows=rows, ncols=3, figsize=(15, 4 * rows * aspect_ratio))

# ---------------------------------------------------
# generate averaged values for anomaly computations
# ---------------------------------------------------

files = os.path.join(idir, 'cams_artic_jul_aug_*.nc')
ds = xr.open_mfdataset(files, concat_dim='time')
param='aod550'

climatology_mean = ds[param].groupby('time.month').mean('time')
climatology_mean.to_netcdf('aod_climatology.nc')
climatology_mean.plot(x='longitude', y='latitude', col='month', col_wrap=3)


plt.figure(figsize=(15, 7))
ax = plt.subplot(1, 2, 1, projection=crs)
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_feature(land_feat)
ax.grid()
ax.gridlines()
ax.coastlines('50m', linewidth=0.5)
ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
ax.coastlines()
p = climatology_mean.isel(month=1).plot(ax=ax, transform=ccrs.PlateCarree(), cmap=caero,
                                cbar_kwargs=dict(pad=.1, aspect=20, shrink=0.6))

climatology_std = ds.groupby('time.month').std('time')
stand_anomalies = xr.apply_ufunc(lambda x, m, s: (x - m) / s,
                                ds.groupby('time.month'),
                                climatology_mean, climatology_std)

stand_anomalies.mean('location').to_dataframe()[['tmin', 'tmax']].plot()

file_aod_ave = 'aod_ave_'+str(years[0])+'_'+str(years[-1])+'.nc'
file_ice_ave = 'ice_ave_'+str(years[0])+'_'+str(years[-1])+'.nc'

first = True
for idx, year in enumerate(years):
    print(idx,year)
    file = os.path.join(idir, 'cams_artic_jul_aug_' + str(year) + '.nc')
    file_ice = os.path.join(idir, 'era5_ice_artic_jul_aug_' + str(year) + '.nc')

    ds = xr.open_dataset(file)
    ds_ice = xr.open_dataset(file_ice)

    # ----- aod
    aod_mean = ds.aod550.mean(dim=('time'))
    # aod_sum = ds.aod550.sum(dim=('time'))

    # ----- ice
    ice_mean = ds_ice.siconc.mean(dim=('time'))


    if first:
        aod_ave = aod_mean
        ice_ave =ice_mean
        first=False
    else:
        aod_ave = aod_ave+aod_mean
        ice_ave = ice_ave+ice_mean

aod_ave_ = aod_ave / (idx+1)
ice_ave_ = ice_ave / (idx+1)

aod_ave_.to_netcdf(file_aod_ave)
ice_ave_.to_netcdf(file_ice_ave)

# ---------------------------------------------------
# plot timeseries for absolute and anomaly values
# --------------------------------------------------
for idx, year in enumerate(years):
    print(idx,year)
    file = os.path.join(idir, 'cams_artic_jul_aug_' + str(year) + '.nc')
    file_ice = os.path.join(idir, 'era5_ice_artic_jul_aug_' + str(year) + '.nc')

    ds = xr.open_dataset(file)

    ds_ice = xr.open_dataset(file_ice)

    if (generate_daily):
        for i in range(0, ds.time.shape[0], 4):
            print(i)
            plt.figure(figsize=(15, 7))
            ax = plt.subplot(1, 2, 1, projection=crs)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.add_feature(land_feat)
            ax.grid()
            ax.gridlines()
            ax.coastlines('50m', linewidth=0.5)
            ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
            ax.coastlines()
            p = ds.aod550.isel(time=i).plot(ax=ax, transform=ccrs.PlateCarree(), cmap=caero,
                                            cbar_kwargs=dict(pad=.1, aspect=20, shrink=0.6))
            p.set_clim(0, 1)

            ax = plt.subplot(1, 2, 2, projection=crs)
            ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
            ax.coastlines()
            ds.siconc.isel(time=i).plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cice,
                                        cbar_kwargs=dict(pad=.1, aspect=20, shrink=0.6))

            plt.savefig(os.path.join(ofig, 'aot_ice_from_cams_artic_' + str(ds.time[i].values)[:10] + '.png'))
            plt.close()

    # ----------------------------------
    # plot climatology data
    # ----------------------------------

    # ----- aod
    aod_mean = ds.aod550.mean(dim=('time'))
    aod_sum = ds.aod550.sum(dim=('time'))

    # ----- ice
    ice_mean = ds_ice.siconc.mean(dim=('time'))

    # ----------------------------------
    # format data into timeseries
    # ----------------------------------
    lat_min1, lat_max1, lon_min1, lon_max1 = 70, 80, 90, 160
    lat_min2, lat_max2, lon_min2, lon_max2 = 80, 87, 90, 170

    roi1 = [[lon_min1, lat_min1], [lon_max1, lat_min1], [lon_max1, lat_max1], [lon_min1, lat_max1]]
    roi2 = [[lon_min2, lat_min2], [lon_max2, lat_min2], [lon_max2, lat_max2], [lon_min2, lat_max2]]
    roi = [roi1, roi2]
    id = [0, 1]
    name = ['roi_laptev', 'roi_laptev_north']
    abbrev = ['roi1', 'roi2']
    mask = regionmask.Regions_cls('roi', id, name, abbrev, roi)
    mask_ = mask.mask(ds, lon_name='longitude', lat_name='latitude')
    mask_ice = mask.mask(ds_ice, lon_name='longitude', lat_name='latitude')

    # ----- aod
    # ts_aod = ds.aod550.mean(dim=('latitude','longitude'))
    aod = ds.aod550.where(mask_ == 0)
    ts_aod_roi1 = aod.mean(dim=('latitude', 'longitude'))
    ts_aod25_roi1 = aod.quantile(0.25, dim=('latitude', 'longitude'))
    ts_aod75_roi1 = aod.quantile(0.75, dim=('latitude', 'longitude'))
    aod = ds.aod550.where(mask_ == 1)
    ts_aod25_roi2 = aod.quantile(0.25, dim=('latitude', 'longitude'))
    ts_aod75_roi2 = aod.quantile(0.75, dim=('latitude', 'longitude'))
    ts_aod_roi2 = aod.mean(dim=('latitude', 'longitude'))

    # ----- ice
    ice = ds_ice.siconc.where(mask_ice == 0)
    ts_ice_roi1 = ice.mean(dim=('latitude', 'longitude'))
    ice = ds_ice.siconc.where(mask_ice == 1)
    ts_ice_roi2 = ice.mean(dim=('latitude', 'longitude'))
    # ts_ice25_roi2 = ice.quantile(0.25, dim=('latitude', 'longitude'))
    # ts_ice75_roi2 = ice.quantile(0.75, dim=('latitude', 'longitude'))

    # ----------------------------------
    #  plot mean and time series*
    # ----------------------------------

    plt.figure(figsize=(15, 7))
    # create subplot grid
    G = gridspec.GridSpec(2, 8, left=0.01, wspace=0.25, hspace=0.25)

    # ----- aod
    ax = plt.subplot(G[0, :3], projection=crs)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    # ax.add_feature(land_feat)
    ax.grid()
    ax.gridlines()
    ax.coastlines('50m', linewidth=0.5)

    mask.plot(ax=ax, regions=[0, 1], add_ocean=False, coastlines=False, label='abbrev', )
    p = aod_mean.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=caero, cbar_kwargs=dict(pad=.01, aspect=20, shrink=0.8))
    p.set_clim(0, 0.6)

    ax = plt.subplot(G[0, 3:])

    ts_aod_roi1.plot(label='roi1')
    plt.fill_between(ts_aod_roi1.time.values, ts_aod25_roi1.values, ts_aod75_roi1.values, alpha=.4)
    ts_aod_roi2.plot(label='roi2')
    plt.fill_between(ts_aod_roi2.time.values, ts_aod25_roi2.values, ts_aod75_roi2.values, alpha=.4)
    plt.legend(ncol=2)
    plt.ylim([0,1])
    # ----- ice
    ax = plt.subplot(G[1, :3], projection=crs)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(land_feat)
    ax.grid()
    ax.gridlines()
    ax.coastlines('50m', linewidth=0.5)


    # mask.plot(ax=ax, regions=[0, 1], add_ocean=False, coastlines=False, label='abbrev', )
    p = ice_mean.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cice, cbar_kwargs=dict(pad=.01, aspect=20, shrink=0.8))
    # p.set_clim(0, 0.6)

    ax = plt.subplot(G[1, 3:])
    ts_ice_roi1.plot(label='roi1')
    ts_ice_roi2.plot(label='roi2')
    # plt.fill_between(ts_ice_roi2.time.values, ts_ice25_roi2.values, ts_ice75_roi2.values, alpha=.4)
    plt.legend(ncol=2)
    plt.ylim([0, 1])
    plt.suptitle('Aerosol and ice; Jul-Aug ' + str(year))
    plt.savefig(os.path.join(ofig, 'aot_ice_from_cams_era5_artic_' + str(year) + '.png'), dpi=300)
    plt.close()

aod_mean_ = aod_mean_/(idx+1)
ice_mean_ = ice_mean_/(idx+1)