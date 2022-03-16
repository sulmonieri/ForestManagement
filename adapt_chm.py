# -*- coding: utf-8 -*-
"""
Desc: Script to manipulate shapefiles (canopy height models) for land use change experiments
3 Options implemented to remove trees: 'manual' , 'auto', 'random'
'manual' allows user to select area of interest to remove all trees within selected perimeter (e.g. wind-throw events)
'auto' allows user to specify a number, based on which every xth tree will be removed
'random ' allows user to specify fraction of forest which will be randomly removed

Created on 20.01.22 09:28
@author: malle
"""

import gdal
import numpy as np
import geopandas as gpd
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path as Path1
import matplotlib.pyplot as plt
from shapely.geometry import Point
import rasterio
import fiona
import pandas as pd
import rasterio.mask
from osgeo import ogr
import random
import time
from pathlib import Path
import click

USE_CLI_ARGS = True  # set to True if running from the command line, set to False if running from PyCharm


class ManualSelect:

    def __init__(self, ax1, collection, alpha_other=0.1):
        self.canvas = ax1.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        self.fc = collection.get_facecolors()
        self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax1, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path1(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def update_chm(sel_pts, pycrown_out):
    """
    Find all crowns + tree tops that need to be eliminated and update/create TIFs accordingly

    Parameters
    ==========
    sel_pts : np.array
        selected points to cut from canopy height model

    pycrown_out : Path Object
        path to input/output folder

    Returns
    =======
    top_cor_cut : np.array
        new array of all points of top of the canopy
    crown_rast_cut: np.array
        new array of all crown polygons
    """

    driver = ogr.GetDriverByName("ESRI Shapefile")

    top_cor_new = pycrown_out / "tree_location_top_cor_new.shp"
    chm_pc_adapt = pycrown_out / "chm_new.tif"
    crown_rast_new = pycrown_out / "tree_crown_poly_raster_new.shp"

    crown_rast = pycrown_out / "tree_crown_poly_raster.shp"
    crown_rast_all = gpd.GeoDataFrame.from_file(str(crown_rast))

    top_cor = pycrown_out / "tree_location_top_cor.shp"
    top_cor_all = gpd.GeoDataFrame.from_file(str(top_cor))

    ids_crown = []
    ids_top = []
    for num in range(len(sel_pts)):
        idc = crown_rast_all.contains(Point(sel_pts[num, :]))
        idc1 = np.where(idc)[0]
        if idc1.size != 0:
            ids_crown.append(idc1[0])
        id_t = top_cor_all.contains(Point(sel_pts[num, :]))
        id_t1 = np.where(id_t)[0]
        ids_top.append(id_t1[0])

    # cut geopandas files to exclude eliminated trees
    top_cor_cut = top_cor_all.loc[~top_cor_all['DN'].isin(np.array(ids_top))]
    crown_rast_cut = crown_rast_all.loc[~crown_rast_all['DN'].isin(np.array(ids_crown))]

    # Remove output shapefiles if they already exist
    if crown_rast_new.exists:
        driver.DeleteDataSource(str(crown_rast_new))
    if top_cor_new.exists:
        driver.DeleteDataSource(str(top_cor_new))
    # Write new shapefiles:
    top_cor_cut.to_file(top_cor_new)
    crown_rast_cut.to_file(crown_rast_new)

    # cut CHM:
    crown_rast = pycrown_out / "tree_crown_poly_raster.shp"
    with fiona.open(crown_rast, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    shapes_series = pd.Series(shapes)

    chm_pc = pycrown_out / "chm.tif"
    with rasterio.open(chm_pc) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes_series[ids_crown], crop=False, invert=True)
        out_meta = src.meta

    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    with rasterio.open(chm_pc_adapt, "w", **out_meta) as dest:
        dest.write(out_image)

    return top_cor_cut, crown_rast_cut


def raster2array(geotif_file):
    """
    Convert raster to array

    Parameters
    ==========
    geotif_file : path to .tif
        geotif of CHM

    Returns
    =======
    chm_array : array
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
    metadata['geotransform'] = dataset.GetGeoTransform()

    mapinfo = dataset.GetGeoTransform()
    metadata['pixelWidth'] = mapinfo[1]
    metadata['pixelHeight'] = mapinfo[5]

    metadata['ext_dict'] = {}
    metadata['ext_dict']['xMin'] = mapinfo[0]
    metadata['ext_dict']['xMax'] = mapinfo[0] + dataset.RasterXSize/mapinfo[1]
    metadata['ext_dict']['yMin'] = mapinfo[3] + dataset.RasterYSize/mapinfo[5]
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
        array[array == int(metadata['noDataValue'])] = np.nan
        array = array/metadata['scaleFactor']
        return array, metadata

    elif metadata['bands'] > 1:
        print('More than one band ... need to modify function for case of multiple bands')


def plot_figs(top_cor_cut, crown_rast_all, crown_rast_cut, x, y, pycrown_out, fig_comp):

    chm_pc_adapt = pycrown_out / "chm_new.tif"
    chm_array_new, chm_array_metadata_new = raster2array(str(chm_pc_adapt))

    chm_pc = pycrown_out / "chm.tif"
    chm_array, chm_array_metadata = raster2array(str(chm_pc))

    fig1 = plt.figure(figsize=(13, 4))
    ax_old = fig1.add_subplot(121)
    ax_old.set_title('Original CHM \n #Trees = '+str(len(x)), fontsize=9)
    plt.imshow(chm_array, extent=chm_array_metadata['extent'], alpha=1)
    cbar1 = plt.colorbar(fraction=0.046, pad=0.04)
    cbar1.set_label('Canopy height [m]', rotation=270, labelpad=20)
    crown_rast_all.plot(ax=ax_old, color='w', alpha=0.6, edgecolor='r', linewidth=0.5)
    ax_old.scatter(x, y, s=1, marker='x', color='r', linewidth=0.2)

    ax_new = fig1.add_subplot(122)
    ax_new.set_title('Modified CHM\n #Trees = '+str(len(top_cor_cut['geometry'].x)), fontsize=9)
    plt.imshow(chm_array_new, extent=chm_array_metadata_new['extent'], alpha=1)
    cbar2 = plt.colorbar(fraction=0.046, pad=0.04)
    cbar2.set_label('Canopy height [m]', rotation=270, labelpad=20)
    crown_rast_cut.plot(ax=ax_new, color='w', alpha=0.6, edgecolor='r', linewidth=0.4)
    ax_new.scatter(top_cor_cut['geometry'].x, top_cor_cut['geometry'].y, s=1, marker='x', color='r', linewidth=0.2)
    fig1.savefig(fig_comp, dpi=350, bbox_inches='tight')
    plt.close()


def manual_cutting(pycrown_out, crown_rast_all, x, y, *_):

    chm_pc = pycrown_out / "chm.tif"
    chm_array, chm_array_metadata = raster2array(str(chm_pc))

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    plt.imshow(chm_array, extent=chm_array_metadata['extent'], alpha=0.6)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    crown_rast_all.plot(ax=ax, color='w', alpha=0.6, edgecolor='r', linewidth=0.4)
    pts = ax.scatter(x, y, s=5, marker='x', color='r')
    ax.axis('equal')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    selector = ManualSelect(ax, pts)

    def accept(event):
        if event.key == "enter":
            global selected_pts1
            selected_pts1 = np.array(selector.xys[selector.ind].data)
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()
            plt.close()
            plt.ioff()

    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Press enter to accept selected points.")

    plt.show()
    selected_pts = selected_pts1
    return selected_pts


def random_cutting(_0, _1, _2, _3, top_cor_all, random_fraction_cut, _4):
    all_trees = top_cor_all['geometry']
    x_coords = all_trees.x
    y_coords = all_trees.y

    ids_all = np.array(range(0, len(x_coords) - 1))
    num_to_select = int(len(x_coords) * random_fraction_cut)
    list_of_random_items = random.sample(list(ids_all), num_to_select)

    sel_pts_x = x_coords[np.array(list_of_random_items)]
    sel_pts_y = y_coords[np.array(list_of_random_items)]
    selected_pts = np.transpose(np.array([sel_pts_x, sel_pts_y]))
    return selected_pts


def auto_cutting(_0, _1, _2, _3, top_cor_all, _4, amount_trees_cut):
    all_trees = top_cor_all['geometry']
    x_coords = all_trees.x
    y_coords = all_trees.y
    sel_pts_x = x_coords[::amount_trees_cut]
    sel_pts_y = y_coords[::amount_trees_cut]
    selected_pts = np.transpose(np.array([sel_pts_x, sel_pts_y]))
    return selected_pts


def main(cut_trees_method, amount_trees_cut, random_fraction_cut, path_data):
    tt = time.time()
    timeit1 = 'Selected points to cut successfully [{:.3f}s]'
    timeit2 = 'Cut & output CHM/crowns/tops successfully [{:.3f}s]'
    timeit3 = 'Output figure generated and saved [total time = {:.3f}s]'

    driver = ogr.GetDriverByName("ESRI Shapefile")

    pycrown_out = Path(path_data)
    crown_rast = pycrown_out / "tree_crown_poly_raster.shp"
    crown_rast_all = gpd.GeoDataFrame.from_file(str(crown_rast))

    if cut_trees_method == 'manual':
        fig_comp = pycrown_out / ("comp_chm_"+cut_trees_method+".png")
        cutting_method = manual_cutting
    elif cut_trees_method == 'random':
        fig_comp = pycrown_out / ("comp_chm_"+cut_trees_method+"_"+str(random_fraction_cut)+".png")
        cutting_method = random_cutting
    elif cut_trees_method == 'auto':
        fig_comp = pycrown_out / ("comp_chm_"+cut_trees_method+"_"+str(amount_trees_cut)+".png")
        cutting_method = auto_cutting

    top_cor = pycrown_out / "tree_location_top_cor.shp"
    top_cor_all = gpd.GeoDataFrame.from_file(str(top_cor))
    datasource_tops = driver.Open(str(top_cor), 0)
    lyr_tops = datasource_tops.GetLayer()

    x = []
    y = []
    for row in lyr_tops:
        geom = row.geometry()
        x.append(geom.GetX())
        y.append(geom.GetY())

    selected_pts = cutting_method(pycrown_out, crown_rast_all, x, y, top_cor_all, random_fraction_cut, amount_trees_cut)

    print(timeit1.format(time.time() - tt))
    top_cor_cut, crown_rast_cut = update_chm(selected_pts, pycrown_out)

    print(timeit2.format(time.time() - tt))
    plot_figs(top_cor_cut, crown_rast_all, crown_rast_cut, x, y, pycrown_out, fig_comp)

    print(timeit3.format(time.time() - tt))
    return None


@click.command()
@click.option('--cut_trees_method', help='auto or random or manual [str]')
@click.option('--amount_trees_cut', default=None, type=float, help='only needs to be set if auto - '
                                                                   'every xth tree to be cut [float]')
@click.option('--random_fraction_cut', default=None, type=float, help='only needs to be set if random - '
                                                                      'fraction of dataset to be cut [float]')
@click.option('--path_in', help='input path [str]')
def cli(cut_trees_method, amount_trees_cut, random_fraction_cut, path_in):
    main(cut_trees_method, amount_trees_cut, random_fraction_cut, path_in)


if __name__ == '__main__':

    if USE_CLI_ARGS:
        cli()
    else:
        cut_trees_method = 'random'  # options: 'manual' , 'auto', 'random'-
        amount_trees_cut = 2  # if using auto setting - every xth tree to cut
        random_fraction_cut = 0.2  # if using random setting - which fraction of all trees should be cut?

        path_in = '/home/malle/pycrown/experiment_sites/BD1/result/dalponteCIRC_numba_10Mrad_ws4_chm4'

        main(cut_trees_method, amount_trees_cut, random_fraction_cut, path_in)
