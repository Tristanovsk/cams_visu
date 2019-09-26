import os
import cmocean as cm

import cartopy as cpy
import cartopy.crs as ccrs

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import xarray as xr
import regionmask

from matplotlib.backends.backend_pdf import PdfPages

generate_daily=False
ofig = '/local/AIX/tristan.harmel/project/ardyna/cams/fig/'

land_feat = cpy.feature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face', facecolor=cpy.feature.COLORS['land'])

crs = ccrs.NearsidePerspective(100, 71)
extent = [-180, 180, 65, 90]  # [90,180, 71, 78]

cmap = cm.tools.crop_by_percent(cm.cm.delta,30,which='both')
cice = cm.tools.crop_by_percent(cm.cm.ice, 25, which='min')
caero = cm.tools.crop_by_percent(cmap, 20, which='max')
year = 2014

aspect_ratio =1
years = np.arange(2004, 2018)
rows = len(years)
fig_clima, axs_clima = plt.subplots(nrows=rows, ncols=3, figsize=(15, 4 * rows * aspect_ratio))

for idx, year in enumerate(years):
    file = '/local/AIX/tristan.harmel/project/ardyna/cams/data/cams_artic_jul_aug_' + str(year) + '.nc'
    file_ice = '/local/AIX/tristan.harmel/project/ardyna/cams/data/era5_ice_artic_jul_aug_' + str(year) + '.nc'

    ds = xr.open_dataset(file)

    ds_ice = xr.open_dataset(file_ice)

    if(generate_daily):
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
    G = gridspec.GridSpec(2, 8, left=0.01, wspace=0.25,hspace=0.25)

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

# ----- ice
    ax = plt.subplot(G[1, :3], projection=crs)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(land_feat)
    ax.grid()
    ax.gridlines()
    ax.coastlines('50m', linewidth=0.5)

    #mask.plot(ax=ax, regions=[0, 1], add_ocean=False, coastlines=False, label='abbrev', )
    p = ice_mean.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cice, cbar_kwargs=dict(pad=.01, aspect=20, shrink=0.8))
    #p.set_clim(0, 0.6)

    ax = plt.subplot(G[1, 3:])
    ts_ice_roi1.plot(label='roi1')
    ts_ice_roi2.plot(label='roi2')
    #plt.fill_between(ts_ice_roi2.time.values, ts_ice25_roi2.values, ts_ice75_roi2.values, alpha=.4)
    plt.legend(ncol=2)
    plt.suptitle('Aerosol and ice; Jul-Aug '+str(year))
    plt.savefig(os.path.join(ofig, 'aot_ice_from_cams_era5_artic_' + str(year) + '.png'), dpi=300)
    plt.close()



