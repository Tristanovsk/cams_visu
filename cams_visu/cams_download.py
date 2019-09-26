import os, sys
import numpy as np

def download_erainterim(target, date, time='00:00:00', grid='0.125/0.125',
                        param='137.128/151.128/206.210/207.210/213.210/214.210/215.210/216.210',
                        area=None, data_type=''):
    ''' This function open a connexion through an existing CAMS/ERAIterim account and download the requested data.


    Arguments:

        * ``target`` -- path name for output file
        * ``time`` -- time of forecast or reanalysis
        * ``date`` -- acquisition date (ex: '2015-01-31')
        * ``grid`` -- resolution in degrees (ex: '0.125/0.125')
        * ``param`` -- ID numbers for the requsted data
        * ``area`` -- option to restrict area (for global dataset None or '90/-180/-90/180')


    Notes:
        * ``server`` -- connexion au server via API ECMWF
        * ``step`` -- forecast time hours
    '''

    from cams_visu.api import ECMWFDataServer

    server = ECMWFDataServer()
    if area is None:
        area = '90/-180/40/180'

    step = '0'
    if data_type == 'cams_forecast':
        class_ = 'mc'
        dataset = 'cams_nrealtime'
        time = '00:00:00'
        step = '0/6/12/18'
        type = 'fc'


    elif data_type == 'cams_reanalysis':
        class_ = 'mc'
        dataset = 'cams_reanalysis'
        time = '00:00:00/06:00:00/12:00:00/18:00:00'
        type = 'an'

    elif data_type == 'interim':
        class_ = 'ei'
        dataset = 'interim'
        type = 'an'
    else:
        #print('Error: not appropriate dataset for ecmwf/cams download')
        class_ = 'mc'
        dataset = data_type
        time = '00:00:00/06:00:00/12:00:00/18:00:00'
        type = 'an'
        #sys.exit()
    server.retrieve({ \
            'class': class_, \
            'dataset': dataset, \
            'date': date, \
            'grid': grid, \
            'levtype': 'sfc', \
            'param': param, \
            'step': step, \
            'stream': 'oper', \
            'time': time, \
            'type': type, \
            'format': 'netcdf',
            'area': area,
            'target': target})
    try:
        server.retrieve({ \
            'class': class_, \
            'dataset': dataset, \
            'date': date, \
            'grid': grid, \
            'levtype': 'sfc', \
            'param': param, \
            'step': step, \
            'stream': 'oper', \
            'time': time, \
            'type': type, \
            'format': 'netcdf',
            'area': area,
            'target': target})
    except:
        print('Error: not appropriate cams settings for download:')
        print({ \
            'class': class_, \
            'dataset': dataset, \
            'date': date, \
            'grid': grid, \
            'levtype': 'sfc', \
            'param': param, \
            'step': step, \
            'stream': 'oper', \
            'time': time, \
            'type': type, \
            'format': 'netcdf',
            'area': area,
            'target': target})
        sys.exit()

    server = None
    return

for year in np.arange(2004,2018):

    # target='/local/AIX/tristan.harmel/project/ardyna/cams/data/cams_artic_jul_aug_'+str(year)+'.nc'
    # if os.path.isfile(target): continue
    # print('downloading CAMS files...' + str(target))
    # startDate = '%04d%02d%02d' % (year,7, 1)
    # lastDate = '%04d%02d%02d' % (year,8,31)
    # param='125.210/137.128/151.128/165.128/166.128/167.128/206.128/207.210/213.210/214.210/215.210/216.210/31.128/34.128'
    # data_type='cams_reanalysis'
    # requestDates = startDate + '/TO/' + lastDate
    # download_erainterim(str(target), requestDates, param=param,data_type=data_type)


    target='/local/AIX/tristan.harmel/project/ardyna/cams/data/era5_ice_artic_jul_aug_'+str(year)+'.nc'
    if os.path.isfile(target): continue
    print('downloading ERA5 files...' + str(target))

    param='sea_ice_cover'
    data_type='reanalysis-era5-single-levels'


    import cdsapi
    time = '00:00:00/06:00:00/12:00:00/18:00:00'
    c = cdsapi.Client()
    c.retrieve(data_type,
               {
                   "variable": param,
                   "product_type": "reanalysis",
                   'month':[
            '07','08'
            ],
            'day':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12',
                '13','14','15',
                '16','17','18',
                '19','20','21',
                '22','23','24',
                '25','26','27',
                '28','29','30',
                '31'
            ],
            'year': str(year),

                   "time": time,
                   "format": "netcdf"
               },
               target)
