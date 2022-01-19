import subprocess
import os
import csv
from osgeo import gdal
from pathlib import Path
import shutil

os.environ['PROJ_LIB'] = '/home/malle/miniconda3/envs/SetUp_sites/share/proj' #somehow without this it can't find the proj.db file...

#define paths etc.
wrk_dir = '/home/malle/slfhome/Postdoc2/experiment_sites'
sites  = os.listdir(wrk_dir)
pycrown_dir = '/home/malle/pycrown/experiment_sites'

#swiss wide chm & dtm can be loaded outside of loop:
chm_ch = '/media/malle/LaCie/Canopy_height_model/VHM/VHM_2021-08-12.tif'
dtm_ch = '/media/malle/LaCie/Canopy_height_model/swissALTI3D_5M_CHLV95_LN02_2020.tif'

in_chm = gdal.Open(chm_ch)
in_band_chm = in_chm.GetRasterBand(1)  # chm only has 1 band
in_gt_chm = in_chm.GetGeoTransform()  # get coordinates etc.
inv_gt_chm = gdal.InvGeoTransform(in_gt_chm)
if inv_gt_chm is None:
    raise RuntimeError('Inverse geotransform failed')

in_dtm = gdal.Open(dtm_ch)
in_band_dtm = in_dtm.GetRasterBand(1)  #chm only has 1 band
in_gt_dtm = in_dtm.GetGeoTransform() #get coordinates etc.
inv_gt_dtm = gdal.InvGeoTransform(in_gt_dtm)
if inv_gt_dtm is None:
    raise RuntimeError('Inverse geotransform failed')

for site in sites:
    site = os.path.basename(site)
    chm_cut = os.path.join(wrk_dir, site,"CHM.tif")
    dtm_cut_orig = os.path.join(wrk_dir, site, site + "_dtm_cut_orig.tif")
    dtm_cut = os.path.join(wrk_dir, site,"DTM.tif")
    dsm_cut = os.path.join(wrk_dir, site,"DSM.tif")

    out_pycrown = os.path.join(pycrown_dir, site,'data')
    Path(out_pycrown).mkdir(parents=True, exist_ok=True)

    #read coordinate file
    coord_file = os.path.join(wrk_dir, site, "coord.csv")
    xys = []
    with open(coord_file) as fp:
        reader = csv.reader(fp)
        next(reader)
        for row in reader:
            xys.append([float(n) for n in row[:2]])
    upper_left, lower_right = xys[0], xys[1]
    ul_x, ul_y = upper_left[0],upper_left[1]
    lr_x, lr_y = lower_right[0],lower_right[1]

    offsets_ul = gdal.ApplyGeoTransform(inv_gt_chm, ul_x, ul_y)
    offsets_lr = gdal.ApplyGeoTransform(inv_gt_chm, lr_x, lr_y)
    off_ulx, off_uly = map(int, offsets_ul)
    off_lrx, off_lry = map(int, offsets_lr)

    # number rows/columns of new extent
    rows_chm = off_lry - off_uly
    columns_chm = off_lrx - off_ulx

    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(chm_cut, columns_chm,rows_chm,1,gdal.GDT_Float32)

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
    out_band = None  # close properly!
    out_ds = None # close properly!

    # now same for dtm
    offsets_ul = gdal.ApplyGeoTransform(inv_gt_dtm, ul_x, ul_y)
    offsets_lr = gdal.ApplyGeoTransform(inv_gt_dtm, lr_x, lr_y)
    off_ulx, off_uly = map(int, offsets_ul)
    off_lrx, off_lry = map(int, offsets_lr)

    rows = off_lry - off_uly
    columns = off_lrx - off_ulx
    out_ds_dtm = gtiff_driver.Create(dtm_cut_orig, columns, rows, 1,gdal.GDT_Float32)
    out_ds_dtm.SetProjection(in_dtm.GetProjection())
    subset_ulx, subset_uly = gdal.ApplyGeoTransform(in_gt_dtm, off_ulx, off_uly)
    out_gt = list(in_gt_dtm)
    out_gt[0] = subset_ulx
    out_gt[3] = subset_uly
    out_ds_dtm.SetGeoTransform(out_gt)

    out_band = out_ds_dtm.GetRasterBand(1)
    data_dtm = in_band_dtm.ReadAsArray(off_ulx, off_uly, columns, rows)
    out_band.WriteArray(data_dtm)
    out_band.FlushCache()   # write to disk
    out_band = None  # close properly!
    out_ds_dtm = None  # close properly!

    #now also resample dtm to higher res! - use gdal lib - need to make sure file does not exist yet!!
    if os.path.exists(dtm_cut):
        os.remove(dtm_cut)

    args = ['gdalwarp','-tr','1.0','1.0','-r','bilinear',dtm_cut_orig,dtm_cut]
    result = subprocess.call(args)

    if os.path.exists(dsm_cut):
        os.remove(dsm_cut)

    args1 = ['gdal_calc.py','-A',dtm_cut,'-B',chm_cut,'--outfile',dsm_cut,'--calc','A+B']
    result1 = subprocess.call(args1)

    # now copy chm,dtm,dsm pycrown data directory
    shutil.copy(dsm_cut,out_pycrown)
    shutil.copy(dtm_cut, out_pycrown)
    shutil.copy(chm_cut, out_pycrown)