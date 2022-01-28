# -*- coding: utf-8 -*-
"""
Desc: 
Created on 27.01.22 08:55
@author: malle
"""


from datetime import datetime

from pycrown import PyCrown

if __name__ == '__main__':

    TSTART = datetime.now()

    F_CHM = '/home/malle/pycrown/experiment_sites/BD1/data/CHM.tif'
    F_DTM = '/home/malle/pycrown/experiment_sites/BD1/data/DTM.tif'
    F_DSM = '/home/malle/pycrown/experiment_sites/BD1/data/DSM.tif'
    F_LAS = ''

    PC = PyCrown(F_CHM, F_DTM, F_DSM, outpath='/home/malle/pycrown/experiment_sites/BD1/result/dalponteCIRC_numba_10Mrad_ws4_chm4_thseed02')

    # Cut off edges
    # PC.clip_data_to_bbox((1802200, 1802400, 5467250, 5467450))

    # Smooth CHM with 5m median filter
    PC.filter_chm(4, ws_in_pixels=True)

    # Tree Detection with local maximum filter
    PC.tree_detection(PC.chm, ws=4, ws_in_pixels=True, hmin=3.)

    # Clip trees to bounding box (no trees on image edge)
    # original extent: 1802140, 1802418, 5467295, 5467490
    # PC.clip_trees_to_bbox(bbox=(1802150, 1802408, 5467305, 5467480))
    # PC.clip_trees_to_bbox(bbox=(1802160, 1802400, 5467315, 5467470))
    PC.clip_trees_to_bbox(inbuf=0.5)  # inward buffer of 11 metre

    # Crown Delineation
    PC.crown_delineation(algorithm='dalponteCIRC_numba', th_tree=3.,
                         th_seed=0.2, th_crown=0.55, max_crown=10.)

    # Correct tree tops on steep terrain
    PC.correct_tree_tops()

    # Calculate tree height and elevation
    PC.get_tree_height_elevation(loc='top')
    PC.get_tree_height_elevation(loc='top_cor')

    # Screen small trees
    PC.screen_small_trees(hmin=2., loc='top')

    # Convert raster crowns to polygons
    PC.crowns_to_polys_raster()
    # PC.crowns_to_polys_smooth(store_las=False)

    # Check that all geometries are valid
    PC.quality_control()

    # Export results
    PC.export_raster(PC.chm, PC.outpath / 'chm.tif', 'CHM')
    PC.export_tree_locations(loc='top')
    PC.export_tree_locations(loc='top_cor')
    PC.export_tree_crowns(crowntype='crown_poly_raster')
    # PC.export_tree_crowns(crowntype='crown_poly_smooth')

    TEND = datetime.now()

    print(f"Number of trees detected: {len(PC.trees)}")
    print(f'Processing time: {TEND-TSTART} [HH:MM:SS]')
