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

from osgeo import gdal
import numpy as np
import geopandas as gpd
from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path as Path1
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import rasterio
import pandas as pd
import rasterio.mask
from osgeo import ogr
import random
import time
from pathlib import Path
import click
import shutil
import rioxarray
import os
import subprocess


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

        self.poly = PolygonSelector(ax1, self.onselect,
                                    props=dict(color='g', alpha=1),
                                    handle_props=dict(mec='g', mfc='g', alpha=1))
        self.path = None
        self.ind = []

    def onselect(self, verts):
        path = Path1(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()
        self.path = path

    def disconnect(self):
        self.poly.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()



def update_chm(sel_pts, pycrown_out, name_chm, cut_trees_method1, buffer, forest_mask, poly_2cut, crown_rast_all,
               top_cor_all):
    """
    Find all crowns + tree tops that need to be eliminated and update/create TIFs accordingly
    Parameters
    ==========
    sel_pts : np.array
        selected points to cut from canopy height model
    pycrown_out : Path Object
        path to input/output folder
    name_chm : str
        add on for CHM name
    cut_trees_method1 : str
        manual/auto/random?
    Returns
    =======
    top_cor_cut : np.array
        new array of all points of top of the canopy
    crown_rast_cut: np.array
        new array of all crown polygons
    """
    
    driver = ogr.GetDriverByName("ESRI Shapefile")

    top_cor_new = pycrown_out / "tree_location_top_cor_new.shp"
    chm_pc_adapt = pycrown_out / Path("chm_" + name_chm + ".tif")
    crown_rast_new = pycrown_out / "tree_crown_poly_raster_new.shp"

    # cut CHM based on selected tree crown rasters:
    crown_rast = pycrown_out / "tree_crown_poly_raster.shp"
    shapefile = gpd.read_file(crown_rast)

    top_cor_all.reset_index(drop=True, inplace=True)
    crown_rast_all.reset_index(drop=True, inplace=True)

    df_coords = pd.DataFrame({'x': sel_pts[:, 0], 'y': sel_pts[:, 1]})
    df_coords['coords'] = list(zip(df_coords['x'], df_coords['y']))
    df_coords['coords'] = df_coords['coords'].apply(Point)
    no_for = gpd.GeoSeries(df_coords['coords'])
    to_mask_layer23 = crown_rast_all["geometry"].apply(lambda x: no_for.within(x).any())
    to_mask_layer231 = top_cor_all["geometry"].apply(lambda x: no_for.intersects(x).any())

    crown_rast_cut = crown_rast_all.drop(crown_rast_all.index[to_mask_layer23])
    top_cor_cut = top_cor_all.drop(top_cor_all.index[to_mask_layer231])

    # Remove output shapefiles if they already exist
    if crown_rast_new.is_file():
        driver.DeleteDataSource(str(crown_rast_new))
    if top_cor_new.is_file():
        driver.DeleteDataSource(str(top_cor_new))
    # Write new shapefiles:
    top_cor_cut.to_file(top_cor_new)
    crown_rast_cut.to_file(crown_rast_new)

    # if manual, lay a convex hull around outermost dimension of selected tree crown polygons and cut out everything
    # within it; if automatic/random, only cut out trees that are above 10m (selected in previous step), and apply a 1m
    # buffer around each tree top polygon
    if cut_trees_method1 == 'manual':
        a2 = ([[p[0], p[1]] for p in np.array(poly_2cut)[0]])
        hull_to_cut = pd.Series(Polygon(a2))
        cut_out_chm = hull_to_cut
    else:
        to_mask_layer23a = shapefile["geometry"].apply(lambda x: no_for.within(x).any())
        shapefile_cut_10m = shapefile.drop(shapefile.index[~to_mask_layer23a])
        shapefile_cut_10m_buffer = shapefile_cut_10m.buffer(buffer)
        cut_out_chm = pd.Series(shapefile_cut_10m_buffer)

    chm_pc = pycrown_out / "CHM_noPycrown.tif"
    with rasterio.open(chm_pc) as src:
        out_image, out_transform = rasterio.mask.mask(src, cut_out_chm, crop=False, invert=True)
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
#        metadata['scaleFactor'] = raster.GetScale()

        # band statistics
        metadata['bandstats'] = {}  # make a nested dictionary to store band stats in same
        stats = raster.GetStatistics(True, True)
        metadata['bandstats']['min'] = round(stats[0], 2)
        metadata['bandstats']['max'] = round(stats[1], 2)
        metadata['bandstats']['mean'] = round(stats[2], 2)
        metadata['bandstats']['stdev'] = round(stats[3], 2)

        array = dataset.GetRasterBand(1).ReadAsArray(0, 0, metadata['array_cols'],
                                                     metadata['array_rows']).astype(float)
#        array[array == int(metadata['noDataValue'])] = np.nan
#        array = array/metadata['scaleFactor']
        return array, metadata

    elif metadata['bands'] > 1:
        print('More than one band ... need to modify function for case of multiple bands')


def plot_figs(top_cor_cut, crown_rast_all, crown_rast_cut, x, y, pycrown_out, fig_comp, name_chm):

    chm_pc_adapt = pycrown_out / Path("chm_" + name_chm + ".tif")
    chm_array_new, chm_array_metadata_new = raster2array(str(chm_pc_adapt))
    chm_array_new[chm_array_new < 1] = np.nan
    ex1a = chm_array_metadata_new['extent']
    ex2a = ex1a[0] + 200, ex1a[1] - 200, ex1a[2] + 200, ex1a[3] - 200

    chm_pc = pycrown_out / "CHM_noPycrown.tif"
    chm_array, chm_array_metadata = raster2array(str(chm_pc))
    chm_array[chm_array < 1] = np.nan
    ex1 = chm_array_metadata['extent']
    ex2 = ex1[0] + 200, ex1[1] - 200, ex1[2] + 200, ex1[3] - 200

    dtm_pc = pycrown_out.parents[0] / "DTM.tif"
    dtm_array, dtm_array_metadata = raster2array(str(dtm_pc))
    ex1a = dtm_array_metadata['extent']
    ex2b = ex1a[0] + 200, ex1a[1] - 200, ex1a[2] + 200, ex1a[3] - 200

    fig1 = plt.figure(figsize=(14, 6))
    ax_old = fig1.add_subplot(121)
    ax_old.set_title('Original CHM \n #Trees = '+str(len(x)), fontsize=9)
    plt.imshow(chm_array[200:-200, 200:-200], extent=ex2, alpha=1)
    crown_rast_all.plot(ax=ax_old, facecolor='none', alpha=0.7, edgecolor='cyan', linewidth=0.4)
    plt.contour(np.flipud(dtm_array[200:-200, 200:-200]), extent=ex2b, linewidths=1.3, colors="red",
                levels=list(range(0, 3000, 100)))
    plt.contour(np.flipud(dtm_array[200:-200, 200:-200]), extent=ex2b, linewidths=0.7, colors="red",
                levels=list(range(0, 3000, 10)))
    # ax_old.scatter(x, y, s=1, marker='x', color='r', linewidth=0.2)

    ax_new = fig1.add_subplot(122)
    ax_new.set_title('Modified CHM\n #Trees = '+str(len(top_cor_cut['geometry'].x)), fontsize=9)
    plt.imshow(chm_array_new[200:-200, 200:-200], extent=ex2a, alpha=1)
    cbar2 = plt.colorbar(fraction=0.046, pad=0.04)
    cbar2.set_label('Canopy height [m]', rotation=270, labelpad=20)
    fig1.savefig(fig_comp, dpi=350, bbox_inches='tight')
    plt.close()
    # 3) for analysis I also need tif @10m resolution
    chm_pc_adapt_10m = pycrown_out / Path("chm_" + name_chm + "_10m.tif")
    if os.path.exists(chm_pc_adapt_10m):
        os.remove(chm_pc_adapt_10m)
    # use gdal lib (file can't exist yet!)
    print("10m res tif ready")
    args = ['gdalwarp', '-tr', '10.0', '10.0', '-r', 'near', chm_pc_adapt, chm_pc_adapt_10m]
    subprocess.call(args)


def manual_cutting(pycrown_out, crown_rast_all_in, x, y, *_):

    chm_pc = pycrown_out / "CHM_noPycrown.tif"
    chm_array, chm_array_metadata = raster2array(str(chm_pc))
    ex1 = chm_array_metadata['extent']
    ex2 = ex1[0] + 200, ex1[1] - 200, ex1[2] + 200, ex1[3] - 200

    dtm_pc = pycrown_out.parents[0] / "DTM.tif"
    dtm_array, dtm_array_metadata = raster2array(str(dtm_pc))
    ex1a = dtm_array_metadata['extent']
    ex2a = ex1a[0] + 200, ex1a[1] - 200, ex1a[2] + 200, ex1a[3] - 200

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    plt.imshow(chm_array[200:-200, 200:-200], extent=ex2, alpha=0.7)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    crown_rast_all_in.plot(ax=ax, facecolor='none', alpha=0.8, edgecolor='cyan', linewidth=0.4)
    pts = ax.scatter(x, y, s=5, marker='x', color='cyan')
    ax.axis('equal')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    plt.contour(np.flipud(dtm_array[200:-200, 200:-200]), extent=ex2a, linewidths=1.9, colors="red",
                levels=list(range(0, 3000, 100)))

    plt.contour(np.flipud(dtm_array[200:-200, 200:-200]), extent=ex2a, linewidths=1.1, colors="red",
                levels=list(range(0, 3000, 10)))

    selector = ManualSelect(ax, pts)
    print("Select points in the figure by enclosing them within a polygon.")
    print("Press the 'esc' key to start a new polygon.")
    print("Hold the 'shift' key to move all of the vertices.")
    print("Hold the 'ctrl' key to move a single vertex.")

    def accept(event):
        if event.key == "enter":
            global selected_pts1, poly_2cut1
            selected_pts1 = np.array(selector.xys[selector.ind].data)
            poly_2cut2 = selector.path
            poly_2cut1 = poly_2cut2.to_polygons(closed_only=True)
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()
            plt.close()
            plt.ioff()

    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Press enter to accept selected points.")
    plt.xlim([np.min(x), np.max(x)])
    plt.ylim([np.min(y), np.max(y)])
    plt.show()
    selected_pts = selected_pts1
    poly_2cut = poly_2cut1
    return selected_pts, poly_2cut


def random_cutting(_0, _1, _2, _3, top_cor_all, random_fraction_cut, _4):
    all_trees = top_cor_all[top_cor_all['TH'] > 8]['geometry']  # only select trees >10m
    all_trees.index = np.arange(len(all_trees))
    all_trees.reset_index()
    x_coords = all_trees.x
    y_coords = all_trees.y

    print(top_cor_all['TH'])
    print(x_coords)
    print(np.size(y_coords))
    print(np.size(top_cor_all))

    ids_all = np.array(range(len(x_coords)))
    num_to_select = int(len(x_coords) * random_fraction_cut)
    list_of_random_items = random.sample(list(ids_all), num_to_select)
    print(list_of_random_items)

    sel_pts_x = x_coords[np.array(list_of_random_items)]
    sel_pts_y = y_coords[np.array(list_of_random_items)]
    selected_pts = np.transpose(np.array([sel_pts_x, sel_pts_y]))
    return selected_pts, None


def auto_cutting(_0, _1, _2, _3, top_cor_all, _4, amount_trees_cut):
    all_trees = top_cor_all[top_cor_all['TH'] > 10]['geometry']  # only select trees >10m
    all_trees.index = np.arange(len(all_trees))
    x_coords = all_trees.x
    y_coords = all_trees.y
    sel_pts_x = x_coords[::amount_trees_cut]
    sel_pts_y = y_coords[::amount_trees_cut]
    selected_pts = np.transpose(np.array([sel_pts_x, sel_pts_y]))
    return selected_pts, None


def main(cut_trees_method, amount_trees_cut, random_fraction_cut, path_data, buffer, forest_mask, buffer_peri):
    #buffer = np.int32(buffer)
    #buffer_peri = np.int32(buffer_peri)
    tt = time.time()
    timeit1 = 'Selected points to cut successfully [{:.3f}s]'
    timeit2 = 'Cut & output CHM/crowns/tops successfully [{:.3f}s]'
    timeit3 = 'Output figure generated and saved [total time = {:.3f}s]'

    driver = ogr.GetDriverByName("ESRI Shapefile")

    pycrown_out = Path(path_data)

    orig_chm = pycrown_out.parents[0] / 'CHM.tif'
    shutil.copy(orig_chm, pycrown_out / "CHM_noPycrown.tif")

    chm_array1, chm_array_metadata1 = raster2array(str(pycrown_out / "CHM_noPycrown.tif"))
    chm_array1[chm_array1 < 1] = np.nan

    fig = plt.figure(figsize=(10, 10))
    plt.imshow(chm_array1, extent=chm_array_metadata1['extent'], alpha=1.0)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    fig_comp1 = pycrown_out / "CHM_noPycrown.png"
    fig.savefig(fig_comp1, dpi=350, bbox_inches='tight')
    plt.close()

    chm_array2, chm_array_metadata2 = raster2array(str(pycrown_out / "chm.tif"))
    chm_array2[chm_array2 < 1] = np.nan

    fig = plt.figure(figsize=(10, 10))
    plt.imshow(chm_array2, extent=chm_array_metadata2['extent'], alpha=1.0)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    fig_comp1 = pycrown_out / "CHM_Pycrown.png"
    fig.savefig(fig_comp1, dpi=350, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 10))
    plt.imshow(chm_array2-chm_array1, extent=chm_array_metadata2['extent'], alpha=1.0)
    plt.clim(np.min(chm_array2-chm_array1), np.max(chm_array2-chm_array1))
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Canopy height [m]', rotation=270, labelpad=20)
    fig_comp1 = pycrown_out / "CHM_Pycrown_diff.png"
    fig.savefig(fig_comp1, dpi=350, bbox_inches='tight')
    plt.close()

    crown_rast = pycrown_out / "tree_crown_poly_raster.shp"
    crown_rast_all = gpd.GeoDataFrame.from_file(str(crown_rast))

    top_cor = pycrown_out / "tree_location_top_cor.shp"
    top_cor_all = gpd.GeoDataFrame.from_file(str(top_cor))
    ext11 = chm_array_metadata2['extent']
    ext_nocrop = ext11[0] + buffer_peri, ext11[0] + buffer_peri, ext11[1] - buffer_peri, ext11[1] - buffer_peri, \
        ext11[2] + buffer_peri, ext11[3] - buffer_peri, ext11[3] - buffer_peri, ext11[2] + buffer_peri
    poly_cut = Polygon((list(zip(ext_nocrop[0:4], ext_nocrop[4:]))))

    if forest_mask == 1:
        forest_mask_in = pycrown_out.parents[0] / 'Forest_mask_10m.tif'
        xds_fm = rioxarray.open_rasterio(forest_mask_in)
        a45 = np.where(xds_fm.values == 128, np.nan, xds_fm.values)
        xds_fm.values = a45
        fm_df = xds_fm.to_dataframe(name="mask_value").reset_index()
        to_mask = fm_df[np.isnan(fm_df['mask_value'])]

        df_coords = pd.DataFrame({'x': to_mask['x'], 'y': to_mask['y']})
        df_coords['coords'] = list(zip(to_mask['x'], to_mask['y']))
        df_coords['coords'] = df_coords['coords'].apply(Point)

        no_for = gpd.GeoSeries(df_coords['coords'])
        no_for_buf = no_for.buffer(10)

        to_mask_layer23 = top_cor_all["geometry"].apply(lambda x: no_for_buf.contains(x).any())
        to_mask_layer231 = crown_rast_all["geometry"].apply(lambda x: no_for_buf.contains(x.centroid).any())

        crown_rast_all.drop(crown_rast_all.index[to_mask_layer231], inplace=True)
        top_cor_all.drop(top_cor_all.index[to_mask_layer23], inplace=True)
        print(np.size(top_cor_all))
        print(np.size(crown_rast_all))

    # cut extra buffer from input data: (needed so we actually only cut out trees in the area we are interested in...
    crown_rast_all1 = crown_rast_all[crown_rast_all.geometry.centroid.within(poly_cut)]
    top_cor_all1 = top_cor_all[top_cor_all.geometry.within(poly_cut)]
    x = list(top_cor_all1.geometry.apply(lambda p: p.x))
    y = list(top_cor_all1.geometry.apply(lambda p: p.y))

    if cut_trees_method == 'manual':
        name_chm = cut_trees_method+"_fm"+str(forest_mask)
        cutting_method = manual_cutting
    elif cut_trees_method == 'random':
        name_chm = cut_trees_method+"_"+str(random_fraction_cut)+"_fm"+str(forest_mask)+"_buffer"+str(buffer)+"m"
        cutting_method = random_cutting
    elif cut_trees_method == 'auto':
        name_chm = cut_trees_method+"_"+str(amount_trees_cut)+"_fm"+str(forest_mask)+"_buffer"+str(buffer)+"m"
        cutting_method = auto_cutting
    else:
        name_chm = None
        cutting_method = None

    fig_comp = pycrown_out / ("comp_chm_" + name_chm + ".png")
    selected_pts, selected_path = cutting_method(pycrown_out, crown_rast_all1, x, y, top_cor_all1, random_fraction_cut,
                                                 amount_trees_cut)

    print(timeit1.format(time.time() - tt))
    top_cor_cut, crown_rast_cut = update_chm(selected_pts, pycrown_out, name_chm, cut_trees_method, buffer, forest_mask,
                                             selected_path, crown_rast_all1, top_cor_all1)

    print(timeit2.format(time.time() - tt))
    plot_figs(top_cor_cut, crown_rast_all1, crown_rast_cut, x, y, pycrown_out, fig_comp, name_chm)

    print(timeit3.format(time.time() - tt))
    return None


@click.command()
@click.option('--cut_trees_method', help='auto or random or manual [str]')

@click.option('--amount_trees_cut', default=None, type=float, help='only needs to be set if auto - '
                                                                   'every xth tree to be cut [float]')

@click.option('--random_fraction_cut', default=None, type=float, help='only needs to be set if random - '
                                                                      'fraction of dataset to be cut [float]')

@click.option('--path_in', help='input path [str]')

@click.option('--buffer', type=float, help='buffer in meters around individual trees to cut (not necessary if pycrown delineation went well) [float]')

@click.option('--buffer_peri', type=float, help='buffer in meters to perimeter of aoi [float]')

@click.option('--forest_mask', help='use forest mask [1] or not [0] [int]')

def cli(cut_trees_method, amount_trees_cut, random_fraction_cut, path_in, buffer, forest_mask, buffer_peri):
    main(cut_trees_method, amount_trees_cut, random_fraction_cut, path_in, buffer, forest_mask, buffer_peri)


if __name__ == '__main__':

    if USE_CLI_ARGS:
        cli()
    else:
        cut_trees_method = 'random'  # options: 'manual' , 'auto', 'random'-
        amount_trees_cut = 3  # if using auto setting - every xth tree to cut (if random/manual - ignore)
        random_fraction_cut = 0.25  # if using random setting - which fraction of all trees should be cut? (if auto/manual - ignore)
        buffer = 0  # if wanting to add buffer around each individual tree crown [if pycrown delination is done well this should not be necessary]
        buffer_peri = 200  # meters added to perimeter of BDM site (not to be incorporated into this analysis, but for transmissivity calculations)
        forest_mask = 1  # set to 0 or 1 => select x percent of forest within forest mask only [default: 1]

        path_in = '/home/vagrant/ForestManagement/aoi/adapt_chm/results/' # path to pycrown output
        main(cut_trees_method, amount_trees_cut, random_fraction_cut, path_in, buffer, forest_mask, buffer_peri)
