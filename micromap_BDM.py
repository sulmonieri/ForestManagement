# -*- coding: utf-8 -*-
"""
Desc: 
Created on 03.05.22 17:31
@author: malle
"""

import pandas as pd
import numpy as np
import datetime
from osgeo import gdal
import ray
import psutil
import time
import statsmodels.formula.api as smf
from pathlib import Path
import matplotlib.pyplot as plt
import os


def raster2array(geotif_file):
    """
    Convert raster to array

    Parameters
    ==========
    geotif_file : path to .tif
        geotif of CHM

    Returns
    =======
    asp_array : array
        new array of all points in .tif
    chm_array_metadata: string
        all metadata used to convert raster2array
    """

    metadata = {}
    dataset = gdal.Open(geotif_file)
    metadata['array_rows'] = dataset.RasterYSize
    metadata['array_cols'] = dataset.RasterXSize
    metadata['bands'] = dataset.RasterCount
    metadata['driver'] = dataset.GetDriver().LongName
    metadata['projection'] = dataset.GetProjection()
    metadata['gt_ch2018'] = dataset.GetGeoTransform()

    mapinfo = dataset.GetGeoTransform()
    metadata['pixelWidth'] = mapinfo[1]
    metadata['pixelHeight'] = mapinfo[5]

    metadata['ext_dict'] = {}
    metadata['ext_dict']['xMin'] = mapinfo[0]
    metadata['ext_dict']['xMax'] = mapinfo[0] + dataset.RasterXSize*mapinfo[1]
    metadata['ext_dict']['yMin'] = mapinfo[3] + dataset.RasterYSize*mapinfo[5]
    metadata['ext_dict']['yMax'] = mapinfo[3]

    metadata['extent'] = (metadata['ext_dict']['xMin'], metadata['ext_dict']['xMax'],
                          metadata['ext_dict']['yMin'], metadata['ext_dict']['yMax'])

    if metadata['bands'] == 1:
        raster = dataset.GetRasterBand(1)
        metadata['noDataValue'] = raster.GetNoDataValue()
        metadata['scaleFactor'] = raster.GetScale()

        # band statistics
        metadata['bandstats'] = {}  # make a nested dictionary to store band stats in same
        stats = raster.GetStatistics(True, True)
        metadata['bandstats']['min'] = round(stats[0], 2)
        metadata['bandstats']['max'] = round(stats[1], 2)
        metadata['bandstats']['mean'] = round(stats[2], 2)
        metadata['bandstats']['stdev'] = round(stats[3], 2)

        array = dataset.GetRasterBand(1).ReadAsArray(0, 0, metadata['array_cols'],
                                                     metadata['array_rows']).astype(float)
        array[array == (metadata['noDataValue'])] = np.nan
        array = array/metadata['scaleFactor']
        return array, metadata

    elif metadata['bands'] > 1:
        print('More than one band ... need to modify function for case of multiple bands')


if __name__ == '__main__':
    TSTART = datetime.datetime.now()

    transm_use = 1
    date_int = '-06-06'
    year_int = '2046'
    scenario = 2  # set to 1 if scenario is RCP8.5, 0 if scenario is RCP2.6, 2 if scenario is today
    forest_scenario = 'summer'

    # read in dataset
    dfs_low = pd.read_csv('/home/malle/eric_micromap/sum_autumn_mod_dfs_low.csv')

    if transm_use == 0:
        mod = smf.mixedlm('Tmax ~ Tmax_meteo + topo_index + topo_wetness + vegh + aspect_n + slope', dfs_low,
                          groups=dfs_low['region'])
    else:
        mod = smf.mixedlm('Tmax ~ Tmax_meteo * trans_Tmax + topo_index + topo_wetness + vegh + aspect_n + slope',
                          dfs_low, groups=dfs_low['region'])

        mod_tmin = smf.mixedlm('Tmin ~ Tmin_meteo + skyview + topo_index + topo_wetness + vegh + aspect_n + slope',
                          dfs_low, groups=dfs_low['region'])
        result_tmin=mod_tmin.fit()

    result = mod.fit()
    print(result.summary())

    print("Forest Scenario: ",forest_scenario)

    if scenario == 0:
        print("Climate Scenario: RCP2.6")
    elif scenario ==1:
        print("Climate Scenario: RCP8.5")
    elif scenario == 2:
        print("Climate Scenario: Today")

    wrk_dir = '/home/malle/slfhome/Postdoc2/experiment_sites_select'
    sites = os.listdir(wrk_dir)
    site = 'BDM_1'  # os.path.basename(sites[0])

    rasters_bf = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/'+site+'/PredRasters')
    topo_index_file = rasters_bf / 'tpi.tif'
    topo_wetness_file = rasters_bf / 'twi.tif'

    if forest_scenario == 'summer_075_fm1_0mbuffer':
        vegh_file = rasters_bf / 'chm_random_0.75_fm1_buffer0m_10m.tif'
    elif forest_scenario == 'summer_025_fm1_0mbuffer':
        vegh_file = rasters_bf / 'chm_random_0.25_fm1_buffer0m_10m.tif'
    elif forest_scenario == 'summer_050_fm1_0mbuffer':
        vegh_file = rasters_bf / 'chm_random_0.5_fm1_buffer0m_10m.tif'
    else:
   #     vegh_file = rasters_bf / 'vhm.tif'
        vegh_file = rasters_bf / 'CHM_noPycrown_10m.tif'  # also use "my" canopy height model for reference


    aspect_n_file = rasters_bf / 'aspect_n.tif'
    slope_file = rasters_bf / 'slope.tif'
    if scenario == 0:
        temp_file = rasters_bf / 'CH2018_corrected_tasmax_DMI-HIRHAM_ECEARTH_EUR11_RCP26_QMgrid_summer10yrs5_2046_06_06.tif'
    elif scenario == 1:
        temp_file = rasters_bf / 'CH2018_corrected_tasmax_DMI-HIRHAM_ECEARTH_EUR11_RCP85_QMgrid_summer10yrs5_2046_06_06.tif'
    elif scenario == 2:
        temp_file = rasters_bf / 'temp_20210605.tif'

    fi = gdal.Open(str(topo_index_file))
    band = fi.GetRasterBand(1)
    topo_index1 = band.ReadAsArray()
    topo_index = topo_index1.reshape(1, topo_index1.shape[0], topo_index1.shape[1])

    fi = gdal.Open(str(topo_wetness_file))
    band = fi.GetRasterBand(1)
    topo_wetness1 = band.ReadAsArray()
    topo_wetness = topo_wetness1.reshape(1, topo_wetness1.shape[0], topo_wetness1.shape[1])

    fi = gdal.Open(str(vegh_file))
    band = fi.GetRasterBand(1)
    vegh1 = band.ReadAsArray()

    if forest_scenario == 'summer_075_fm1_0mbuffer' or forest_scenario == 'summer_025_fm1_0mbuffer' or forest_scenario == 'summer_050_fm1_0mbuffer' or forest_scenario == 'summer':
        vegh2 = vegh1[20:-20,20:-20]
        vegh = vegh2.reshape(1, vegh2.shape[0], vegh2.shape[1])*100
    else:
        vegh = vegh1.reshape(1, vegh1.shape[0], vegh1.shape[1])

    fi = gdal.Open(str(aspect_n_file))
    band = fi.GetRasterBand(1)
    aspect_n1 = band.ReadAsArray()
    aspect_n = aspect_n1.reshape(1, aspect_n1.shape[0], aspect_n1.shape[1])

    fi = gdal.Open(str(slope_file))
    band = fi.GetRasterBand(1)
    slope1 = band.ReadAsArray()
    slope = slope1.reshape(1, slope1.shape[0], slope1.shape[1])

    fi = gdal.Open(str(temp_file))
    gt = fi.GetGeoTransform()
    proj = fi.GetProjection()
    band = fi.GetRasterBand(1)
    temp1 = band.ReadAsArray()
    temp = temp1.reshape(1, temp1.shape[0], temp1.shape[1])

    if transm_use != 1:
        stack = np.vstack((temp, topo_index, topo_wetness, vegh, aspect_n, slope))
    else:
        transm_bf = Path('/home/malle/transm_calcs/'+site+'/Output_CR_10m_'+forest_scenario+'/OutTifs')
        transm_all = list(transm_bf.glob('**/*.tif'))
        str_match = list(filter(lambda x: date_int in str(x), transm_all))
        transm_file = str_match[0]   # transm_bf / 'transm2020-06-20.tif'
        fi = gdal.Open(str(transm_file))
        band = fi.GetRasterBand(1)
        transm1 = band.ReadAsArray()
        transm = transm1.reshape(1, transm1.shape[0], transm1.shape[1])/100  # output is 0-100, but in model it is 0-1
        stack = np.vstack((temp, transm, topo_index, topo_wetness, vegh, aspect_n, slope))

    id_array = np.transpose(np.array(np.meshgrid(range(stack.shape[1]),
                                                 range(stack.shape[2]))).T.reshape(-1, 2), (1, 0))

    combo_all = np.zeros((stack.shape[0] + id_array.shape[0], id_array.shape[1]), dtype=int)
    combo_all[:2, :] = id_array
    combo_all[2:, :] = stack[:, id_array[0, :], id_array[1, :]]

    # option 1: just 1 function that does it all
    @ray.remote
    def f(lme_result, predictors):
        micro_map_test = (
                100 * mod.predict(lme_result.fe_params, exog=np.insert(predictors[2:] / 100, 0, 1))).astype(np.int16)
        return micro_map_test, predictors[0], predictors[1]

    @ray.remote
    def f_test(lme_result, predictors):
        tt = np.insert(predictors[2:] / 100, 0, 1)
        micro_map_test = (
                100 * lme_result.predict(exog=dict(Tmax_meteo=tt[0], topo_index=tt[1], topo_wetness=tt[2], vegh=tt[3],
                                                   aspect_n=tt[4], slope=tt[5]))).astype(np.int16)
        return micro_map_test, predictors[0], predictors[1]

    @ray.remote
    def f2(lme_result, predictors):
        tt = np.insert(predictors[2:] / 100, 0, 1)
        micro_map_test = (
                100 * lme_result.predict(exog=dict(Tmax_meteo=tt[0], trans_Tmax=tt[1], topo_index=tt[2],
                                                   topo_wetness=tt[3], vegh=tt[4], aspect_n=tt[5],
                                                   slope=tt[6]))).astype(np.int16)
        return micro_map_test, predictors[0], predictors[1]

    @ray.remote
    def f2_test(lme_result, predictors):
        tt1 = np.insert(predictors[2:] / 100, 0, 1)
        tt = np.insert(tt1, 3, 1)
        micro_map_test = (100 * mod.predict(lme_result.fe_params, exog=tt)).astype(np.int16)
        return micro_map_test, predictors[0], predictors[1]

    num_cpus = psutil.cpu_count(logical=False)
    start = time.time()
    ray.init(num_cpus=num_cpus)
    if transm_use == 1:
        micro_map_out = ray.get([f2_test.remote(result, combo_all[:, i]) for i in range(combo_all.shape[1])])
    else:
        micro_map_out = ray.get([f.remote(result, combo_all[:, i]) for i in range(combo_all.shape[1])])

    ray.shutdown()
    print("duration all 1 function test =", time.time() - start)

    micro_map = np.zeros((temp1.shape[0], temp1.shape[1]), dtype=int)

    for num in range(len(micro_map_out)):
        aaa1 = micro_map_out[num]
        micro_map[aaa1[1], aaa1[2]] = aaa1[0]

    if transm_use == 0:
        out_folder = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/'+site+'/MicroMaps/excl_transm')
    else:
        out_folder = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site + '/MicroMaps/incl_transm')

    out_folder.mkdir(parents=True, exist_ok=True)

    if scenario == 0:
        out_micro_map = out_folder / Path('temp_'+year_int+date_int+'_transm' + str(transm_use) + '_RCP26_'+forest_scenario+'.tif')
        out_micro_map_plot = out_folder / Path('temp_'+year_int+date_int+'_transm' + str(transm_use) + '_RCP26_'+forest_scenario+'.png')
    elif scenario == 1:
        out_micro_map = out_folder / Path('temp_'+year_int+date_int+'_transm' + str(transm_use) + '_RCP85_'+forest_scenario+'.tif')
        out_micro_map_plot = out_folder / Path('temp_'+year_int+date_int+'_transm' + str(transm_use) + '_RCP85_'+forest_scenario+'.png')
    elif scenario == 2:
        out_micro_map = out_folder / Path('temp_'+year_int+date_int+'_transm' + str(transm_use) + '_today_'+forest_scenario+'.tif')
        out_micro_map_plot = out_folder / Path('temp_'+year_int+date_int+'_transm' + str(transm_use) + '_today_'+forest_scenario+'.png')

    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    # create new clipped temp file:
    out_ds = driver.Create(str(out_micro_map), micro_map.shape[1], micro_map.shape[0], 1, gdal.GDT_Int16,
                           options=['COMPRESS=ZSTD', 'PREDICTOR=2', 'ZSTD_LEVEL=1'])
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(proj)
    outband = out_ds.GetRasterBand(1)
    outband.WriteArray(micro_map)
    outband.SetNoDataValue(np.nan)
    outband.FlushCache()
    outband = None  # close properly! (important!!)
    out_ds = None  # close properly! (important!!)

    temp_array, temp_array_metadata = raster2array(str(out_micro_map))

    # 1b - plot aspect
    fig1 = plt.figure(figsize=(13, 4))
    ax_old = fig1.add_subplot(111)
    ax_old.set_title('Temp @ '+str(site), fontsize=9)
    plt.imshow(micro_map/100, extent=temp_array_metadata['extent'], alpha=1)
    cbar1 = plt.colorbar(fraction=0.046, pad=0.04)
    cbar1.set_label('Temp [deg C]', rotation=270, labelpad=20)
    fig1.savefig(str(out_micro_map_plot), dpi=350, bbox_inches='tight')
    plt.close()
