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
import math


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


def search_in_file(path_in, searchstring_in, filelist1):
    with open(path_in, 'r') as file:
        if searchstring_in in file.name:
            filelist1.append(path_in.name)
        return filelist1


@ray.remote
def f(lme_result, predictors):
    micro_map_test = (
            100 * mod_tmax.predict(lme_result.fe_params,
                                   exog=np.insert(predictors[2:] / 100, 0, 1))).astype(np.int16)
    return micro_map_test, predictors[0], predictors[1]


@ray.remote
def ray_tmax(lme_result, predictors_all):
    micro_map_pred = np.zeros(predictors_all.shape[1])
    for idd in range(predictors_all.shape[1]):
        predictors = predictors_all[:, idd]
        tt1 = np.insert(predictors[2:], 0, 1)
        tt = np.insert(tt1, 3, 1)
        micro_map_pred[idd] = (mod_tmax.predict(lme_result.fe_params, exog=tt)).astype(np.int16)
    return micro_map_pred


@ray.remote
def ray_tmean(lme_result, predictors_all):
    micro_map_pred = np.zeros(predictors_all.shape[1])
    for idd in range(predictors_all.shape[1]):
        predictors = predictors_all[:, idd]
        tt1 = np.insert(predictors[2:], 0, 1)
        tt = np.insert(tt1, 3, 1)
        micro_map_pred[idd] = (mod_tmean.predict(lme_result.fe_params, exog=tt)).astype(np.int16)
    return micro_map_pred


@ray.remote
def ray_tmin(lme_result, predictors_all):
    micro_map_pred = np.zeros(predictors_all.shape[1])
    for idd in range(predictors_all.shape[1]):
        predictors = predictors_all[:, idd]
        tt = np.insert(predictors[2:], 0, 1)
        micro_map_pred[idd] = (mod_tmin.predict(lme_result.fe_params, exog=tt)).astype(np.int16)
    return micro_map_pred



if __name__ == '__main__':
    TSTART = datetime.datetime.now()
    transm_use = 1
    forest_scenario_all = ['summer','summer_manual']
    model_int = ["high", "low", "soil"]
    site_names = ["BDM_2", "BDM_3", "BDM_1"]
    numbers_sc = [1, 0.5, 0]
    # set to 1 if scenario is RCP8.5, 0.5 if RCP4.5 and 0 if scenario is RCP2.6, 2 if scenario is today


    for forest_scenario in forest_scenario_all:

        for model in model_int:
            start_model = time.time()
            dfs_model = pd.read_csv('/home/malle/eric_micromap/master_mod_dfs_'+model+'.csv')

            if transm_use == 0:
                mod_tmax = smf.mixedlm('Tmax ~ Tmax_meteo + topo_index + topo_wetness + vegh + aspect_n + slope', dfs_model,
                                       groups=dfs_model['region'])

                mod_tmean = smf.mixedlm('Tmean ~ Tmean_meteo + topo_index + topo_wetness + vegh + aspect_n + slope',
                                        dfs_model, groups=dfs_model['region'])

                mod_tmin = smf.mixedlm('Tmin ~ Tmin_meteo + topo_index + topo_wetness + vegh + aspect_n + slope',
                                       dfs_model, groups=dfs_model['region'])

            else:
                mod_tmax = smf.mixedlm('Tmax ~ Tmax_meteo * trans_Tmax + topo_index + topo_wetness + vegh + aspect_n + '
                                       'slope', dfs_model, groups=dfs_model['region'])

                mod_tmean = smf.mixedlm('Tmean ~ Tmean_meteo * trans_Tmax + topo_index + topo_wetness + vegh + aspect_n + '
                                        'slope', dfs_model, groups=dfs_model['region'])

                mod_tmin = smf.mixedlm('Tmin ~ Tmin_meteo + skyview + topo_index + topo_wetness + vegh + aspect_n + slope',
                                       dfs_model, groups=dfs_model['region'])

            lme_tmax = mod_tmax.fit()
            lme_tmean = mod_tmean.fit()
            lme_tmin = mod_tmin.fit()

            for site in site_names:
                start_site = time.time()
                for scenario in numbers_sc:
                    print("Site = " + site + ", Model = " + model + ", Forest Scenario = " + forest_scenario)
                    start_szenario = time.time()
                    if scenario == 0:
                        print("Climate Scenario: RCP2.6")
                        sce = 'RCP26'
                    elif scenario == 0.5:
                        print("Climate Scenario: RCP4.5")
                        sce = 'RCP45'
                    elif scenario == 1:
                        print("Climate Scenario: RCP8.5")
                        sce = 'RCP85'
                    elif scenario == 2:
                        print("Climate Scenario: Today")
                        sce = 'None'
                    else:
                        sce = 'None'

                    wrk_dir = '/home/malle/slfhome/Postdoc2/experiment_sites_select'

                    rasters_bf = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/'+site+'/PredRasters')
                    topo_index_file = rasters_bf / 'tpi.tif'
                    topo_wetness_file = rasters_bf / 'twi.tif'
                    aspect_n_file = rasters_bf / 'aspect_n.tif'
                    slope_file = rasters_bf / 'slope.tif'
                    svf_file = Path('/home/malle/transm_calcs/' + site + '/Output_CR_10m_' + forest_scenario +
                                    '/Plots/SVF.tif')

                    if forest_scenario == 'summer_075_fm1_0mbuffer':
                        vegh_file = rasters_bf / 'chm_random_0.75_fm1_buffer0m_10m.tif'
                    elif forest_scenario == 'summer_025_fm1_0mbuffer':
                        vegh_file = rasters_bf / 'chm_random_0.25_fm1_buffer0m_10m.tif'
                    elif forest_scenario == 'summer_050_fm1_0mbuffer':
                        vegh_file = rasters_bf / 'chm_random_0.5_fm1_buffer0m_10m.tif'
                    elif forest_scenario == 'summer_manual':
                        vegh_file = rasters_bf / 'chm_manual_fm1_10m.tif'
                    elif forest_scenario == 'summer':
                        vegh_file = rasters_bf / 'CHM_noPycrown_10m.tif'  # also use "my" canopy height model for reference

                    fi = gdal.Open(str(topo_index_file))
                    band = fi.GetRasterBand(1)
                    topo_index1 = band.ReadAsArray()
                    topo_index = topo_index1.reshape(1, topo_index1.shape[0], topo_index1.shape[1]) / 100

                    fi = gdal.Open(str(topo_wetness_file))
                    band = fi.GetRasterBand(1)
                    topo_wetness1 = band.ReadAsArray()
                    topo_wetness = topo_wetness1.reshape(1, topo_wetness1.shape[0], topo_wetness1.shape[1]) / 100

                    fi = gdal.Open(str(svf_file))
                    band = fi.GetRasterBand(1)
                    svf1 = band.ReadAsArray()
                    skyview = svf1.reshape(1, svf1.shape[0], svf1.shape[1]) / 100

                    fi = gdal.Open(str(vegh_file))
                    band = fi.GetRasterBand(1)
                    vegh1 = band.ReadAsArray()

                    if forest_scenario == 'summer_075_fm1_0mbuffer' or forest_scenario == 'summer_025_fm1_0mbuffer' or \
                            forest_scenario == 'summer_050_fm1_0mbuffer' or forest_scenario == 'summer' or forest_scenario == 'summer_manual':
                        vegh2 = vegh1[20:-20, 20:-20]
                        vegh = vegh2.reshape(1, vegh2.shape[0], vegh2.shape[1])
                    else:
                        vegh = vegh1.reshape(1, vegh1.shape[0], vegh1.shape[1]) / 100

                    fi = gdal.Open(str(aspect_n_file))
                    band = fi.GetRasterBand(1)
                    aspect_n1 = band.ReadAsArray()
                    aspect_n = aspect_n1.reshape(1, aspect_n1.shape[0], aspect_n1.shape[1]) / 100

                    fi = gdal.Open(str(slope_file))
                    band = fi.GetRasterBand(1)
                    slope1 = band.ReadAsArray()
                    slope = slope1.reshape(1, slope1.shape[0], slope1.shape[1]) / 100

                    transm_bf = Path('/home/malle/transm_calcs/' + site + '/Output_CR_10m_' + forest_scenario +
                                     '/OutTifs_monthly')
                    transm_all = list(transm_bf.glob('**/*.tif'))

                    temp_loc_overall = Path('/home/malle/johanna_micromap/input_BDM_CH2018/')

                    num_cpus = psutil.cpu_count(logical=False)
                    if ray.is_initialized() is False:
                        ray.init(num_cpus=num_cpus)

                    for date_num in range(0, 1529, 1):

                        yr_id_frac, yr_id_whole = math.modf(date_num/153)
                        if date_num > 764:
                            real_date = datetime.date(2041+math.floor(yr_id_whole)+44, 6, 1) +\
                                        datetime.timedelta(days=153*yr_id_frac)
                        else:
                            real_date = datetime.date(2041+math.floor(yr_id_whole), 6, 1) +\
                                        datetime.timedelta(days=153*yr_id_frac)

                        if transm_use == 0:
                            out_folder_tmax = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                   '/MicroMaps/excl_transm/' + model + '/tmax')
                            out_folder_tmean = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                    '/MicroMaps/excl_transm/' + model + '/tmean')
                            out_folder_tmin = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                   '/MicroMaps/excl_transm/' + model + '/tmin')

                        else:
                            out_folder_tmax = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                   '/MicroMaps/incl_transm/' + model + '/' + forest_scenario + '/tmax')
                            out_folder_tmean = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                    '/MicroMaps/incl_transm/' + model + '/' + forest_scenario + '/tmean')
                            out_folder_tmin = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                   '/MicroMaps/incl_transm/' + model + '/' + forest_scenario + '/tmin')

                            out_folder_plot_tmax = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                        '/MicroMaps/incl_transm/' + model + '/' +forest_scenario + '/tmax_plot')
                            out_folder_plot_tmean = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                         '/MicroMaps/incl_transm/' + model + '/' +forest_scenario + '/tmean_plot')
                            out_folder_plot_tmin = Path('/home/malle/slfhome/Postdoc2/experiment_sites_select/' + site +
                                                        '/MicroMaps/incl_transm/' + model + '/' + forest_scenario + '/tmin_plot')

                        out_folder_tmax.mkdir(parents=True, exist_ok=True)
                        out_folder_tmean.mkdir(parents=True, exist_ok=True)
                        out_folder_tmin.mkdir(parents=True, exist_ok=True)
                        out_folder_plot_tmax.mkdir(parents=True, exist_ok=True)
                        out_folder_plot_tmean.mkdir(parents=True, exist_ok=True)
                        out_folder_plot_tmin.mkdir(parents=True, exist_ok=True)

                        out_micromap_map_tmax = out_folder_tmax / Path('microtemp_' + sce + '_' +
                                                                       real_date.strftime('%-Y-%m-%d')+'.tif')
                        out_micromap_map_tmean = out_folder_tmean / Path('microtemp_' + sce + '_' +
                                                                         real_date.strftime('%-Y-%m-%d')+'.tif')
                        out_micromap_map_tmin = out_folder_tmin / Path('microtemp_' + sce + '_' +
                                                                       real_date.strftime('%-Y-%m-%d')+'.tif')

                        out_micromap_plot_tmax = out_folder_plot_tmax / Path('microtemp_' + sce + '_' +
                                                                             real_date.strftime('%-Y-%m-%d')+'.png')
                        out_micromap_plot_tmean = out_folder_plot_tmean / Path('microtemp_' + sce+'_' +
                                                                               real_date.strftime('%-Y-%m-%d')+'.png')
                        out_micromap_plot_tmin = out_folder_plot_tmin / Path('microtemp_' + sce + '_' +
                                                                             real_date.strftime('%-Y-%m-%d')+'.png')
                        if out_micromap_map_tmax.is_file() and out_micromap_map_tmean.is_file() and \
                                out_micromap_map_tmin.is_file():
                            pass
                        else:
                            if date_num % 50 == 0:  # only write date every xth iteration
                                print(real_date)
                                ray.shutdown()
                                ray.init(num_cpus=num_cpus)

                            date_int = '-'+f"{real_date.month:02d}"+'-'+f"{real_date.day:02d}"

                            str_match = list(filter(lambda x: date_int[2:3] in str(x), transm_all))

                            transm_file = str_match[0]  # transm_bf / 'transm2020-06-20.tif'
                            fi = gdal.Open(str(transm_file))
                            band = fi.GetRasterBand(1)
                            transm1 = band.ReadAsArray()
                            transm = transm1.reshape(1, transm1.shape[0],
                                                     transm1.shape[1]) / 100  # output is 0-10000, but in model it is 0-100

                            searchstring = '_EUR11_'+sce+'_QMgrid_summer10yrs_ext'+str(date_num+1)+'bdm'+site[-1:]+'.tif'

                            temp_int = ["bdm_tas", "bdm_tasmax", "bdm_tasmin"]
                            for sz in temp_int:
                                temp_loc = temp_loc_overall / Path(sz)
                                dir_content = sorted(temp_loc.iterdir())
                                filelist = []
                                for path in dir_content:
                                    if not path.is_dir():
                                        search_in_file(path, searchstring, filelist)

                                temp_all_stack = []
                                for temp_rast_num in filelist:
                                    if temp_rast_num[-3:] == 'tif':
                                        temp_rast = temp_loc / temp_rast_num

                                        fi = gdal.Open(str(temp_rast))
                                        gt = fi.GetGeoTransform()
                                        proj = fi.GetProjection()
                                        band = fi.GetRasterBand(1)
                                        temp1 = band.ReadAsArray()
                                        temp_all = temp1.reshape(1, temp1.shape[0], temp1.shape[1])
                                        temp_all_stack.append(temp_all)

                                AS = np.stack(temp_all_stack, axis=0)
                                del temp_all_stack
                                temp_used = np.mean(AS, axis=0)
                                del AS

                                if sz == 'bdm_tas':
                                    temp_mean = temp_used
                                elif sz == 'bdm_tasmax':
                                    temp_max = temp_used
                                else:
                                    temp_min = temp_used

                            if transm_use != 1:
                                stack_tmax = np.vstack((temp_max, topo_index, topo_wetness, vegh, aspect_n, slope))
                                stack_tmean = np.vstack((temp_mean, topo_index, topo_wetness, vegh, aspect_n, slope))
                                stack_tmin = np.vstack((temp_min, topo_index, topo_wetness, vegh, aspect_n, slope))
                            else:
                                stack_tmax = np.vstack((temp_max, transm, topo_index, topo_wetness, vegh, aspect_n, slope))
                                stack_tmean = np.vstack((temp_mean, transm, topo_index, topo_wetness, vegh, aspect_n, slope))
                                stack_tmin = np.vstack((temp_min, skyview, topo_index, topo_wetness, vegh, aspect_n, slope))

                            id_array_tmax = np.transpose(np.array
                                                         (np.meshgrid(range(stack_tmax.shape[1]), range(stack_tmax.shape[2])))
                                                         .T.reshape(-1, 2), (1, 0))
                            combo_all_tmax = np.zeros((stack_tmax.shape[0] + id_array_tmax.shape[0], id_array_tmax.shape[1]),
                                                      dtype=float)
                            combo_all_tmax[:2, :] = id_array_tmax
                            combo_all_tmax[2:, :] = stack_tmax[:, id_array_tmax[0, :], id_array_tmax[1, :]]

                            id_array_tmean = np.transpose(np.array
                                                          (np.meshgrid(range(stack_tmean.shape[1]),
                                                                       range(stack_tmean.shape[2]))).T.reshape(-1, 2), (1, 0))
                            combo_all_tmean = np.zeros((stack_tmean.shape[0] + id_array_tmean.shape[0],
                                                        id_array_tmean.shape[1]), dtype=float)
                            combo_all_tmean[:2, :] = id_array_tmean
                            combo_all_tmean[2:, :] = stack_tmean[:, id_array_tmean[0, :], id_array_tmean[1, :]]

                            id_array_tmin = np.transpose(
                                np.array(np.meshgrid(range(stack_tmin.shape[1]),
                                                     range(stack_tmin.shape[2]))).T.reshape(-1, 2), (1, 0))
                            combo_all_tmin = np.zeros((stack_tmin.shape[0] + id_array_tmin.shape[0], id_array_tmin.shape[1]),
                                                      dtype=float)
                            combo_all_tmin[:2, :] = id_array_tmin
                            combo_all_tmin[2:, :] = stack_tmin[:, id_array_tmin[0, :], id_array_tmin[1, :]]

                            start = time.time()
                            if transm_use == 1:
                                num_split = 1000
                                micro_map_out_tmax_raw = [ray_tmax.remote(lme_tmax, combo_all_tmax[:, i:i+num_split])
                                                          for i in range(0, combo_all_tmax.shape[1], num_split)]
                                micro_map_out_tmean_raw = [ray_tmean.remote(lme_tmean, combo_all_tmean[:, i:i+num_split])
                                                           for i in range(0, combo_all_tmean.shape[1], num_split)]
                                micro_map_out_tmin_raw = [ray_tmin.remote(lme_tmin, combo_all_tmin[:, i:i+num_split])
                                                          for i in range(0, combo_all_tmin.shape[1], num_split)]

                                micro_map_out_tmax_comb = np.zeros(combo_all_tmax.shape[1])
                                micro_map_out_tmean_comb = np.zeros(combo_all_tmax.shape[1])
                                micro_map_out_tmin_comb = np.zeros(combo_all_tmax.shape[1])
                                idx_tmax = np.zeros(combo_all_tmax.shape[1])
                                idy_tmax = np.zeros(combo_all_tmax.shape[1])

                                for mii in range(int(combo_all_tmax.shape[1]/num_split)):
                                    part_micromap = ray.get(micro_map_out_tmax_raw[mii])
                                    micro_map_out_tmax_comb[mii*num_split:mii*num_split+num_split] = part_micromap
                                    idx_tmax[mii*num_split:mii*num_split+num_split] = \
                                        combo_all_tmax[0, mii*num_split:mii*num_split+num_split]
                                    idy_tmax[mii*num_split:mii*num_split+num_split] = \
                                        combo_all_tmax[1, mii*num_split:mii*num_split+num_split]

                                    part_micromap = ray.get(micro_map_out_tmean_raw[mii])
                                    micro_map_out_tmean_comb[mii*num_split:mii*num_split+num_split] = part_micromap
                                    part_micromap = ray.get(micro_map_out_tmin_raw[mii])
                                    micro_map_out_tmin_comb[mii*num_split:mii*num_split+num_split] = part_micromap

                                micro_map_out_tmax = list(zip(micro_map_out_tmax_comb.astype(np.int16),
                                                              idx_tmax.astype(np.int16),
                                                              idy_tmax.astype(np.int16)))
                                micro_map_out_tmean = list(zip(micro_map_out_tmean_comb.astype(np.int16),
                                                               idx_tmax.astype(np.int16),
                                                               idy_tmax.astype(np.int16)))
                                micro_map_out_tmin = list(zip(micro_map_out_tmin_comb.astype(np.int16),
                                                              idx_tmax.astype(np.int16),
                                                              idy_tmax.astype(np.int16)))
                            else:
                                micro_map_out_tmax = ray.get([f.remote(lme_tmax, combo_all_tmax[:, i])
                                                              for i in range(combo_all_tmax.shape[1])])
                                micro_map_out_tmean = ray.get([f.remote(lme_tmean, combo_all_tmean[:, i])
                                                               for i in range(combo_all_tmean.shape[1])])
                                micro_map_out_tmin = ray.get([f.remote(lme_tmin, combo_all_tmin[:, i])
                                                              for i in range(combo_all_tmin.shape[1])])

                            # print("Processing time micromap creation day =", time.time() - start)

                            micro_map_tmax = np.zeros((temp1.shape[0], temp1.shape[1]), dtype=int)
                            for num in range(len(micro_map_out_tmax)):
                                aaa1 = micro_map_out_tmax[num]
                                micro_map_tmax[aaa1[1], aaa1[2]] = aaa1[0]

                            micro_map_tmean = np.zeros((temp1.shape[0], temp1.shape[1]), dtype=int)
                            for num in range(len(micro_map_out_tmean)):
                                aaa1 = micro_map_out_tmean[num]
                                micro_map_tmean[aaa1[1], aaa1[2]] = aaa1[0]

                            micro_map_tmin = np.zeros((temp1.shape[0], temp1.shape[1]), dtype=int)
                            for num in range(len(micro_map_out_tmin)):
                                aaa1 = micro_map_out_tmin[num]
                                micro_map_tmin[aaa1[1], aaa1[2]] = aaa1[0]

                            del micro_map_out_tmax, micro_map_out_tmean, micro_map_out_tmin, micro_map_out_tmax_comb, \
                                micro_map_out_tmean_comb, micro_map_out_tmin_comb, part_micromap, micro_map_out_tmax_raw, \
                                micro_map_out_tmean_raw, micro_map_out_tmin_raw

                            driver = gdal.GetDriverByName('GTiff')
                            driver.Register()
                            out_ds = driver.Create(str(out_micromap_map_tmax), micro_map_tmax.shape[1], micro_map_tmax.shape[0],
                                                   1, gdal.GDT_Int16, options=['COMPRESS=ZSTD', 'PREDICTOR=2', 'ZSTD_LEVEL=1'])
                            out_ds.SetGeoTransform(gt)
                            out_ds.SetProjection(proj)
                            outband = out_ds.GetRasterBand(1)
                            outband.WriteArray(micro_map_tmax)
                            outband.SetNoDataValue(np.nan)
                            outband.FlushCache()
                            outband = None  # close properly! (important!!)
                            out_ds = None  # close properly! (important!!)

                            if date_int == '-06-20':
                                temp_array, temp_array_metadata = raster2array(str(out_micromap_map_tmax))
                                fig1 = plt.figure(figsize=(13, 4))
                                ax_old = fig1.add_subplot(111)
                                ax_old.set_title('Temp @ '+str(site), fontsize=9)
                                plt.imshow(micro_map_tmax / 100, extent=temp_array_metadata['extent'], alpha=1)
                                cbar1 = plt.colorbar(fraction=0.046, pad=0.04)
                                cbar1.set_label('Temp [deg C]', rotation=270, labelpad=20)
                                fig1.savefig(str(out_micromap_plot_tmax), dpi=350, bbox_inches='tight')
                                plt.close()

                            # TMEAN
                            out_ds = driver.Create(str(out_micromap_map_tmean), micro_map_tmean.shape[1],
                                                   micro_map_tmean.shape[0],
                                                   1, gdal.GDT_Int16, options=['COMPRESS=ZSTD', 'PREDICTOR=2', 'ZSTD_LEVEL=1'])
                            out_ds.SetGeoTransform(gt)
                            out_ds.SetProjection(proj)
                            outband = out_ds.GetRasterBand(1)
                            outband.WriteArray(micro_map_tmean)
                            outband.SetNoDataValue(np.nan)
                            outband.FlushCache()
                            outband = None  # close properly! (important!!)
                            out_ds = None  # close properly! (important!!)
                            if date_int == '-06-20':
                                temp_array, temp_array_metadata = raster2array(str(out_micromap_map_tmean))
                                fig1 = plt.figure(figsize=(13, 4))
                                ax_old = fig1.add_subplot(111)
                                ax_old.set_title('Temp @ '+str(site), fontsize=9)
                                plt.imshow(micro_map_tmean / 100, extent=temp_array_metadata['extent'], alpha=1)
                                cbar1 = plt.colorbar(fraction=0.046, pad=0.04)
                                cbar1.set_label('Temp [deg C]', rotation=270, labelpad=20)
                                fig1.savefig(str(out_micromap_plot_tmean), dpi=350, bbox_inches='tight')
                                plt.close()

                            # TMIN
                            out_ds = driver.Create(str(out_micromap_map_tmin), micro_map_tmin.shape[1],
                                                   micro_map_tmin.shape[0], 1,
                                                   gdal.GDT_Int16, options=['COMPRESS=ZSTD', 'PREDICTOR=2', 'ZSTD_LEVEL=1'])
                            out_ds.SetGeoTransform(gt)
                            out_ds.SetProjection(proj)
                            outband = out_ds.GetRasterBand(1)
                            outband.WriteArray(micro_map_tmin)
                            outband.SetNoDataValue(np.nan)
                            outband.FlushCache()
                            outband = None  # close properly! (important!!)
                            out_ds = None  # close properly! (important!!)

                            if date_int == '-06-20':
                                temp_array, temp_array_metadata = raster2array(str(out_micromap_map_tmin))
                                fig1 = plt.figure(figsize=(13, 4))
                                ax_old = fig1.add_subplot(111)
                                ax_old.set_title('Temp @ '+str(site), fontsize=9)
                                plt.imshow(micro_map_tmin / 100, extent=temp_array_metadata['extent'], alpha=1)
                                cbar1 = plt.colorbar(fraction=0.046, pad=0.04)
                                cbar1.set_label('Temp [deg C]', rotation=270, labelpad=20)
                                fig1.savefig(str(out_micromap_plot_tmin), dpi=350, bbox_inches='tight')
                                plt.close()


                    print("Processing time for szenario " + sce + " : ", time.time() - start_szenario)
                print("Processing time for site " + site + " : ", time.time() - start_site)
            print("Processing time for model " + model + " : ", time.time() - start_model)

        TEND4 = datetime.datetime.now()
        print(f'Total processing time: {TEND4-TSTART} [HH:MM:SS]')
