# -*- coding: utf-8 -*-
"""
Desc: take pycrown output and manipulate shapefiles/chm
Created on 20.01.22 09:28
@author: malle
"""

import gdal
import numpy as np
import geopandas as gpd
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from shapely.geometry import Point
import rasterio
import fiona
import pandas as pd
import rasterio.mask
import os
import matplotlib.pyplot as plt
from osgeo import ogr
import random


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
        path = Path(verts)
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


def update_chm(sel_pts):
    # find all crowns + tree tops that need to be eliminated (don't have the same id...)
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
    if os.path.exists(crown_rast_new):
        driver.DeleteDataSource(crown_rast_new)
    if os.path.exists(top_cor_new):
        driver.DeleteDataSource(top_cor_new)
    # Write new shapefiles:
    top_cor_cut.to_file(top_cor_new)
    crown_rast_cut.to_file(crown_rast_new)

    # now deal with cutting chm too:
    with fiona.open(crown_rast, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    shapes_series = pd.Series(shapes)

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


def plot_figs(top_cor_cut, crown_rast_cut):

    chm_array_new, chm_array_metadata_new = raster2array(chm_pc_adapt)

    fig1 = plt.figure(figsize=(15, 9))
    print(fig1)
    ax_old = fig1.add_subplot(121)
    ax_old.set_title('Original CHM')
    plt.imshow(chm_array, extent=chm_array_metadata['extent'], alpha=1)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    crown_rast_all.plot(ax=ax_old, color='w', alpha=0.6, edgecolor='r', linewidth=0.5)
    plt.scatter(x, y, s=1, marker='x', color='r', linewidth=0.2)

    ax_new = fig1.add_subplot(122)
    ax_new.set_title('Modified CHM')
    plt.imshow(chm_array_new, extent=chm_array_metadata_new['extent'], alpha=1)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    crown_rast_cut.plot(ax=ax_new, color='w', alpha=0.6, edgecolor='r', linewidth=0.4)
    ax_new.scatter(top_cor_cut['geometry'].x, top_cor_cut['geometry'].y, s=1, marker='x', color='r', linewidth=0.2)

    fig_comp = os.path.join(pycrown_out, "comp_chm"+cut_trees_method+".png")
    print(fig_comp)
    fig1.savefig(fig_comp, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':

    cut_trees_method = 'manual'  # need to update later but for now select here if
    # trees are cut automatically or manual selection .. options: 'manual' , 'auto' -
    # if set to auto also need to define how many trees to cut
    amount_trees_cut = 2  # if using auto setting - how many trees to cut?
    random_fraction_cut = 0.5  # if using random setting - which fraction should be cut?

    driver = ogr.GetDriverByName("ESRI Shapefile")

    pycrown_out = '/home/malle/pycrown/experiment_sites/BD1/result/dalponteCIRC_numba_10Mrad_ws4_chm4'
    chm_pc = os.path.join(pycrown_out, "chm.tif")
    chm_pc_adapt = os.path.join(pycrown_out, "chm_new.tif")
    crown_rast = os.path.join(pycrown_out, "tree_crown_poly_raster.shp")
    crown_rast_new = os.path.join(pycrown_out, "tree_crown_poly_raster_new.shp")
    crown_rast_all = gpd.GeoDataFrame.from_file(crown_rast)

    top_cor = os.path.join(pycrown_out, "tree_location_top_cor.shp")
    top_cor_new = os.path.join(pycrown_out, "tree_location_top_cor_new.shp")
    top_cor_all = gpd.GeoDataFrame.from_file(top_cor)
    fig_comp = os.path.join(pycrown_out, "comp_chm.png")

    dataSource_tops = driver.Open(top_cor, 0)
    dataSource_crown = driver.Open(crown_rast, 0)
    lyr_tops = dataSource_tops.GetLayer()
    lyr_crown = dataSource_crown.GetLayer()

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    chm_array, chm_array_metadata = raster2array(chm_pc)

    plt.imshow(chm_array, extent=chm_array_metadata['extent'], alpha=0.6)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)

    crown_rast_all.plot(ax=ax, color='w', alpha=0.6, edgecolor='r', linewidth=0.4)

    x = []
    y = []
    for row in lyr_tops:
        geom = row.geometry()
        x.append(geom.GetX())
        y.append(geom.GetY())
        # print(row.GetField("DN"), row.GetField("TH"))

    pts = plt.scatter(x, y, s=5, marker='x', color='r')

    plt.axis('equal')
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])

    if cut_trees_method == 'manual':

        selector = ManualSelect(ax, pts)

        def accept(event):
            if event.key == "enter":
                print("Selected points:")
                print(selector.xys[selector.ind])
                selected_pts1 = np.array(selector.xys[selector.ind].data)
                top_cor_cut, crown_rast_cut = update_chm(selected_pts1)
                selector.disconnect()
                ax.set_title("")
                fig.canvas.draw()
                plt.close()
                plt.ioff()
                plot_figs(top_cor_cut, crown_rast_cut)

        fig.canvas.mpl_connect("key_press_event", accept)
        ax.set_title("Press enter to accept selected points.")

        plt.show()
        plt.close()
        plt.ioff()

    elif cut_trees_method == 'auto':
        all_trees = top_cor_all['geometry']
        x_coords = all_trees.x
        y_coords = all_trees.y

        sel_pts_x = x_coords[::amount_trees_cut]
        sel_pts_y = y_coords[::amount_trees_cut]
        selected_pts = np.transpose(np.array([sel_pts_x, sel_pts_y]))
        top_cor_cut, crown_rast_cut = update_chm(selected_pts)
        plot_figs(top_cor_cut, crown_rast_cut)

    elif cut_trees_method == 'random':
        all_trees = top_cor_all['geometry']
        x_coords = all_trees.x
        y_coords = all_trees.y

        ids_all = np.array(range(0, len(x_coords)-1))
        num_to_select = int(len(x_coords) * random_fraction_cut)
        list_of_random_items = random.sample(list(ids_all), num_to_select)

        sel_pts_x = x_coords[np.array(list_of_random_items)]
        sel_pts_y = y_coords[np.array(list_of_random_items)]
        selected_pts = np.transpose(np.array([sel_pts_x, sel_pts_y]))
        top_cor_cut, crown_rast_cut = update_chm(selected_pts)
        plot_figs(top_cor_cut, crown_rast_cut)
