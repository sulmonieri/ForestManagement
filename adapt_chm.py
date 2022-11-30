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

#### random removal of a specified percentage of trees ####
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
    tt = time.time()
    timeit1 = 'Selected points to cut successfully [{:.3f}s]'
    timeit2 = 'Cut & output CHM/crowns/tops successfully [{:.3f}s]'
    timeit3 = 'Output figure generated and saved [total time = {:.3f}s]'

    driver = ogr.GetDriverByName("ESRI Shapefile")

    pycrown_out = Path(path_data)

    orig_chm = pycrown_out.parents[1] / 'data' / 'CHM.tif'
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
        forest_mask_in = pycrown_out.parents[1] / 'data' / 'Forest_mask_10m.tif'
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

@click.option('--buffer', help='auto or random or manual [str]')
@click.option('--buffer_peri', help='auto or random or manual [str]')
@click.option('--forest_mask', help='auto or random or manual [str]')

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
