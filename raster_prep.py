# -*- coding: utf-8 -*-
"""
Calculate and save DTM, DSM, CHM from CH-wide datasets for a given spatial extent
Created on 17.01.22 17:37
@author: malle
"""

import os
import csv
from osgeo import gdal
from pathlib import Path
import shutil
import subprocess
import matplotlib.pyplot as plt
import numpy as np


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


# somehow without this it can't find the proj.db file...
os.environ['PROJ_LIB'] = '/home/malle/miniconda3/envs/SetUp_sites/share/proj'

# define paths etc.
wrk_dir = '/home/malle/slfhome/Postdoc2/experiment_sites_select'
plot_dir = '/home/malle/slfhome/Postdoc2/BD_sites/CHMs_select'
sites = os.listdir(wrk_dir)
pycrown_dir = '/home/malle/pycrown/experiment_sites_select'

# swiss-wide chm & dtm can be loaded outside of loop:
# swiss-wide CHM:
chm_ch = '/media/malle/LaCie/Canopy_height_model/VHM/VHM_2021-08-12.tif'
in_chm = gdal.Open(chm_ch)
in_band_chm = in_chm.GetRasterBand(1)  # chm only has 1 band
in_gt_chm = in_chm.GetGeoTransform()  # get coordinates etc.
inv_gt_chm = gdal.InvGeoTransform(in_gt_chm)  # needed to convert real-world coord. to image coord.
if inv_gt_chm is None:
    raise RuntimeError('Inverse gt_ch2018 failed')

# swiss-wide DTM:
dtm_ch = '/media/malle/LaCie/Canopy_height_model/swissALTI3D_5M_CHLV95_LN02_2020.tif'
in_dtm = gdal.Open(dtm_ch)
in_band_dtm = in_dtm.GetRasterBand(1)  # chm only has 1 band
in_gt_dtm = in_dtm.GetGeoTransform()  # get coordinates etc.
inv_gt_dtm = gdal.InvGeoTransform(in_gt_dtm)
if inv_gt_dtm is None:
    raise RuntimeError('Inverse gt_ch2018 failed')

# swiss-wide forest outline:
forest_mask_file = '/media/malle/LaCie/Canopy_height_model/OSHDForestMask_10m_EPSG2056.tif'
in_fm = gdal.Open(forest_mask_file)
in_band_fm = in_fm.GetRasterBand(1)  # chm only has 1 band
in_gt_fm = in_fm.GetGeoTransform()  # get coordinates etc.
inv_gt_fm = gdal.InvGeoTransform(in_gt_fm)
if inv_gt_fm is None:
    raise RuntimeError('Inverse gt_ch2018 failed')

# needed to create TIFFs later:
gtiff_driver = gdal.GetDriverByName('GTiff')
buffer_m = 200

# all BD sites:
for site in sites:
    site = os.path.basename(site)
    chm_plot = os.path.join(plot_dir, "CHM_"+str(site)+".png")
    # def. all output files/locations:
    chm_cut = os.path.join(wrk_dir, site, "CHM.tif")
    chm_cut_neg = os.path.join(wrk_dir, site, "CHM_incl_neg.tif")
    dtm_cut_orig = os.path.join(wrk_dir, site, site + "_dtm_cut_orig.tif")
    dtm_cut = os.path.join(wrk_dir, site, "DTM.tif")
    dsm_cut = os.path.join(wrk_dir, site, "DSM.tif")
    fm_cut_1m = os.path.join(wrk_dir, site, "Forest_mask_1m.tif")
    fm_cut = os.path.join(wrk_dir, site, "Forest_mask_10m.tif")
    fm_cut_5m = os.path.join(wrk_dir, site, "Forest_mask_5m.tif")
    fm_cut_20m = os.path.join(wrk_dir, site, "Forest_mask_20m.tif")
    fm_plot = os.path.join(wrk_dir, site, "Forest_mask_cut.png")

    out_pycrown = os.path.join(pycrown_dir, site, 'data')
    Path(out_pycrown).mkdir(parents=True, exist_ok=True)

    # read coordinate file
    coord_file = os.path.join(wrk_dir, site, "coord.csv")
    xys = []
    with open(coord_file) as fp:
        reader = csv.reader(fp)
        next(reader)
        for row in reader:
            xys.append([float(n) for n in row[:2]])
    upper_left, lower_right = xys[0], xys[1]
    ul_x, ul_y = upper_left[0] - buffer_m, upper_left[1] + buffer_m
    lr_x, lr_y = lower_right[0] + buffer_m, lower_right[1] - buffer_m

    ul_x_nobuf, ul_y_nobuf = upper_left[0], upper_left[1]
    lr_x_nobuf, lr_y_nobuf = lower_right[0], lower_right[1]

    # 1) CHM:
    offsets_ul = gdal.ApplyGeoTransform(inv_gt_chm, ul_x, ul_y)
    offsets_lr = gdal.ApplyGeoTransform(inv_gt_chm, lr_x, lr_y)
    off_ulx, off_uly = map(int, offsets_ul)
    off_lrx, off_lry = map(int, offsets_lr)
    # number rows/columns of new extent
    rows_chm = off_lry - off_uly
    columns_chm = off_lrx - off_ulx
    # create new clipped CHM file:
    out_ds = gtiff_driver.Create(chm_cut_neg, columns_chm, rows_chm, 1, gdal.GDT_Float32)
    out_ds.SetProjection(in_chm.GetProjection())
    subset_ulx, subset_uly = gdal.ApplyGeoTransform(in_gt_chm, off_ulx, off_uly)
    out_gt = list(in_gt_chm)
    out_gt[0] = subset_ulx
    out_gt[3] = subset_uly
    out_ds.SetGeoTransform(out_gt)
    out_band = out_ds.GetRasterBand(1)
    data_chm = in_band_chm.ReadAsArray(off_ulx, off_uly, columns_chm, rows_chm)
    out_band.WriteArray(data_chm)
    out_band.FlushCache()   # write to disk
    out_band = None  # close properly! (important!!)
    out_ds = None  # close properly! (important!!)

    # this includes negative CHM values -> actually set to 0!
    if os.path.exists(chm_cut):
        os.remove(chm_cut)
    # use gdal lib (file can't exist yet!)
    args1 = ['gdal_calc.py', '-A', chm_cut_neg, '--outfile', chm_cut, '--calc', 'A*(A>=0)', '--NoDataValue', '0']
    result1 = subprocess.call(args1)

    chm_array, chm_array_metadata = raster2array(chm_cut)
    chm_array[chm_array < 1] = np.nan

    # 1b - plot CHM
    fig1 = plt.figure(figsize=(13, 4))
    ax_old = fig1.add_subplot(111)
    ax_old.set_title('Original CHM @ '+str(site), fontsize=9)
    plt.imshow(chm_array, extent=[ul_x_nobuf, lr_x_nobuf, lr_y_nobuf, ul_y_nobuf], alpha=1)
    cbar1 = plt.colorbar(fraction=0.046, pad=0.04)
    cbar1.set_label('Canopy height [m]', rotation=270, labelpad=20)
    fig1.savefig(chm_plot, dpi=350, bbox_inches='tight')
    plt.close()

    # 2) DTM:
    offsets_ul = gdal.ApplyGeoTransform(inv_gt_dtm, ul_x, ul_y)
    offsets_lr = gdal.ApplyGeoTransform(inv_gt_dtm, lr_x, lr_y)
    off_ulx, off_uly = map(int, offsets_ul)
    off_lrx, off_lry = map(int, offsets_lr)
    rows = off_lry - off_uly
    columns = off_lrx - off_ulx
    # create new clipped DTM file:
    out_ds_dtm = gtiff_driver.Create(dtm_cut_orig, columns, rows, 1, gdal.GDT_Float32)
    out_ds_dtm.SetProjection(in_dtm.GetProjection())
    subset_ulx, subset_uly = gdal.ApplyGeoTransform(in_gt_dtm, off_ulx, off_uly)
    out_gt = list(in_gt_dtm)
    out_gt[0] = subset_ulx
    out_gt[3] = subset_uly
    out_ds_dtm.SetGeoTransform(out_gt)
    out_band_dtm = out_ds_dtm.GetRasterBand(1)
    data_dtm = in_band_dtm.ReadAsArray(off_ulx, off_uly, columns, rows)
    out_band_dtm.WriteArray(data_dtm)
    out_band_dtm.FlushCache()   # write to disk
    out_band_dtm = None  # close properly! (important!!)
    out_ds_dtm = None  # close properly! (important!!)

    # 2a) Forest mask (without buffer!):
    offsets_ul = gdal.ApplyGeoTransform(inv_gt_fm, ul_x_nobuf, ul_y_nobuf)
    offsets_lr = gdal.ApplyGeoTransform(inv_gt_fm, lr_x_nobuf, lr_y_nobuf)
    off_ulx, off_uly = map(int, offsets_ul)
    off_lrx, off_lry = map(int, offsets_lr)
    rows = off_lry - off_uly
    columns = off_lrx - off_ulx
    # create new clipped forest mask file:
    out_ds_fm = gtiff_driver.Create(fm_cut, columns, rows, 1, gdal.GDT_Float32)
    out_ds_fm.SetProjection(in_fm.GetProjection())
    subset_ulx, subset_uly = gdal.ApplyGeoTransform(in_gt_fm, off_ulx, off_uly)
    out_gt = list(in_gt_fm)
    out_gt[0] = subset_ulx
    out_gt[3] = subset_uly
    out_ds_fm.SetGeoTransform(out_gt)
    out_band_fm = out_ds_fm.GetRasterBand(1)
    data_fm = in_band_fm.ReadAsArray(off_ulx, off_uly, columns, rows)
    out_band_fm.WriteArray(data_fm)
    out_band_fm.FlushCache()   # write to disk
    out_band_fm = None  # close properly! (important!!)
    out_ds_fm = None  # close properly! (important!!)

    chm_array, chm_array_metadata = raster2array(fm_cut)

    # 2b - plot forest_mask
    fig1 = plt.figure(figsize=(13, 4))
    ax_old = fig1.add_subplot(111)
    ax_old.set_title('Forest mask @ '+str(site), fontsize=9)
    plt.imshow(chm_array, extent=chm_array_metadata['extent'], alpha=1)
    fig1.savefig(fm_plot, dpi=350, bbox_inches='tight')
    plt.close()

    # 3) FM resampled to higher 5m res :
    if os.path.exists(fm_cut_5m):
        os.remove(fm_cut_5m)
    # use gdal lib (file can't exist yet!)
    args = ['gdalwarp', '-tr', '5.0', '5.0', '-r', 'near', fm_cut, fm_cut_5m]
    result = subprocess.call(args)

    if os.path.exists(fm_cut_1m):
        os.remove(fm_cut_1m)
    # use gdal lib (file can't exist yet!)
    args = ['gdalwarp', '-tr', '1.0', '1.0', '-r', 'near', fm_cut, fm_cut_1m]
    result = subprocess.call(args)

    if os.path.exists(fm_cut_20m):
        os.remove(fm_cut_20m)
    # use gdal lib (file can't exist yet!)
    args = ['gdalwarp', '-tr', '20.0', '20.0', '-r', 'near', fm_cut, fm_cut_20m]
    result = subprocess.call(args)

    # 3) DTM resampled to higher CHM res. :
    if os.path.exists(dtm_cut):
        os.remove(dtm_cut)
    # use gdal lib (file can't exist yet!)
    args = ['gdalwarp', '-tr', '1.0', '1.0', '-r', 'bilinear', dtm_cut_orig, dtm_cut]
    result = subprocess.call(args)

    # 4) DSM creation by adding DTM+CHM :
    if os.path.exists(dsm_cut):
        os.remove(dsm_cut)
    # use gdal lib (file can't exist yet!)
    args1 = ['gdal_calc.py', '-A', dtm_cut, '-B', chm_cut, '--outfile', dsm_cut, '--calc', 'A+B']
    result1 = subprocess.call(args1)

    # 5) now copy chm,dtm,dsm pycrown data directory
    shutil.copy(dsm_cut, out_pycrown)
    shutil.copy(dtm_cut, out_pycrown)
    shutil.copy(chm_cut, out_pycrown)
    shutil.copy(fm_cut_1m, out_pycrown)
