# ForestManagement
Set of scripts needed to prep for radiative transfer modelling at biodiversity sites in Switzerland. 

Use raster_prep.py to cut national CHM/DTM and create DSM for spatial extent of interest
User run_pycrown.py to delineate trees and create shapefiles for all tree crowns
Use adapt_chm.py to create a modified CHM based on user input (manual tree cutting, automated tree cutting or random tree cutting)

For adapt_chm.py:
The following python packages are needed: gdal numpy geopandas matplotlib shapely rasterio osgeo random time pathlib click shutil rioxarray os subprocess

If wanting to run from command line, set USE_CLI_ARGS=TRUE in top of script
Otherwise edit settings on bottom of the script
