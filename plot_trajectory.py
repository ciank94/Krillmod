#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import sys
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import os
from analyse_trajectory import ssmu_target


#from Krillmod.get_trajectory import geo2grid, store_traj
#from Krillmod.import_settings import tr_file
#from io import BytesIO

def plot_transit(list_dir, sub_idx):
    # Arguments:
    # Pointer to subset of individuals (sub_key)
    # Use that pointer to identify folder
    # Use target function to identify target area (ssmu_target in this case)- loop over target areas (subs dictionary)
    key_list = list(sub_idx.keys())
    print('Key list for analysis: ')
    print(key_list)
    for sub_key in key_list:
        subs = ssmu_target(sub_key)  # Target regions stored in dictionary (there may be multiple target areas);
        for target_area in subs:
            filepath = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit.npy'
            save_path = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit.png'
            save_path2 = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit_map.png'
            if not os.path.exists(filepath):
                print('Error: Directory ' + filepath + ' with intermediate matrix does not exist, exiting')
                sys.exit()
            else:
                df = np.load(filepath)
                # if sub_key == 'SO' and target_area == 'SG':
                #     breakpoint()
                x = df[:, 0]
                y = df[:, 1]
                t = df[:, 2]
                # Calculate the fraction that reach SOI, the mean time, std, etc., plot histogram;
                x_arrive = x[x > 0.1]
                y_arrive = y[y > 0.1]
                t_arrive = t[t > 0.1]
                t_arrive = t_arrive/24
                if np.shape(t_arrive)[0]==0:
                    print('No particles arrived, skipping')
                else:
                    # mean_t = np.floor(np.mean(t_arrive))
                    # std_t = np.floor(np.std(t_arrive))
                    # print('mean_time:  ' + str(mean_t))
                    # print('standard_deviation_time:  ' + str(std_t))

                    # Transit time distribution:
                    axis1_title = 'number of particles N'
                    axis2_title = 'cumulative N in %'
                    fig, ax1 = plt.subplots()
                    ax2 = ax1.twinx()
                    ax1.set_ylabel(axis1_title, color='b', fontsize=13)
                    ax2.set_ylabel(axis1_title, color='b', fontsize=13)
                    hist1 = ax1.hist(t_arrive, bins=20, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5)
                    cum_vector = np.zeros(len(hist1[1]))
                    cum_val = np.cumsum(hist1[0] / np.sum(hist1[0])) * 100
                    cum_vector[1:len(cum_vector)] = cum_val
                    ax2.plot(hist1[1], cum_vector, 'r', linewidth=2)
                    ax1.set_ylabel(axis1_title, color='b', fontsize=13)
                    ax2.set_ylabel(axis2_title, color='r', fontsize=13)
                    ax1.tick_params(axis='x', labelsize=10)
                    ax1.tick_params(axis='y', labelsize=10)
                    ax2.tick_params(axis='y', labelsize=10)
                    ax1.set_xlabel('Transit days', fontsize=13)
                    plt.grid(alpha=0.45)  # nice and clean grid
                    plt.savefig(save_path, dpi=400)
                    plt.close(fig)

                    # worm plots
                    filepath2 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormx.npy'
                    filepath3 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormy.npy'
                    wormx = np.load(filepath2)
                    wormy = np.load(filepath3)
                    sub_id = np.sum(wormx, axis=1) > 0
                    x_end = df[:, 3]
                    y_end = df[:, 4]
                    xe = x_end[sub_id]
                    ye = y_end[sub_id]
                    wx = wormx[sub_id, :]
                    wy = wormy[sub_id, :]
                    if len(wx) < 1000:
                        sub = 1
                    else:
                        sub = 10
                    skip = (slice(None, None, sub), slice(None, None, sub))
                    skip2 = slice(None, None, sub)
                    wxp = wx[skip]
                    wyp = wy[skip]
                    xep = xe[skip2]
                    yep = ye[skip2]
                    #breakpoint()
                    wxp[wxp==0] = np.nan
                    wyp[wyp == 0] = np.nan

                    #wx[wx==0] = np.nan
                    #wy[wy == 0] = np.nan
                    #breakpoint()

                    grid_lims = dict()
                    depth = np.load(list_dir['depth_file'])
                    grid_lims['idlimx'] = [np.nanmax([np.nanmin(wxp) - 50,0]), np.nanmin([np.nanmax(wxp) + 20, np.shape(depth)[0]])]
                    grid_lims['idlimy'] = [np.nanmax([np.nanmin(wyp) - 100,0]), np.nanmin([np.nanmax(wyp) + 50, np.shape(depth)[1]])]
                    fig, ax = plot_depth(list_dir, grid_lims)


                    for i in range(0, np.shape(wxp)[0]):
                        ax.plot(wxp[i, :], wyp[i, :], '.r-', markersize = 1.5, linewidth=0.8)

                    ax.plot(x_arrive, y_arrive, '.k', markersize = 3, linestyle='None')#linewidth=0.0001)
                    ax.plot(xep, yep, '.w', markersize = 3,linestyle='None')# linewidth=0.0001)
                    plt.savefig(save_path2, dpi=400)
    return


def plot_dom_paths(list_dir, sub_idx):
    depth = np.load(list_dir['depth_file'])
    key_list = list(sub_idx.keys())
    for sub_key in key_list:
        dom_file = list_dir[sub_key + '_folder'] + 'dominant_paths.npy'
        df = np.load(dom_file)
        df[np.isnan(depth)] = np.nan
        list_ids = np.where([df > 0])

        grid_lims = dict()
        grid_lims['idlimx'] = [np.nanmax([np.nanmin(list_ids[2]) - 50, 0]),
                               np.nanmin([np.nanmax(list_ids[2]) + 20, np.shape(depth)[1]])]
        grid_lims['idlimy'] = [np.nanmax([np.nanmin(list_ids[1]) - 100, 0]),
                               np.nanmin([np.nanmax(list_ids[1]) + 50, np.shape(depth)[0]])]
        fig, ax = plot_background(list_dir, grid_lims)
        #breakpoint()
        # NOTE: Do a check for the existence of file
        # Modify for the limits of the plot- hard coding elsewhere- create a dictionary with options for plotting
        # tot_file = sv_dir + tr_folder + '/' + area_name + '_total_visits.npy'
        if sub_key == 'ALL':
            cmax = 0.3
            crange = 0.02
        else:
            cmax = 2.5
            crange = 0.1

        #Dominant pathways
        cmap1 = plt.get_cmap('OrRd')  # Oranges, Reds- sequential coolwarm= divergent, jet, seismic
        data_plot = ax.contourf(df, levels=np.arange(0, cmax, crange), extend='both',cmap=cmap1)
        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cbar = plt.colorbar(data_plot, cax=ax_cb)
        cbar.ax.set_ylabel('Probability (%)', loc='center', size=12.5, weight='bold')
        cbar.ax.tick_params(labelsize=9,rotation=0)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        file = list_dir[sub_key + '_folder'] + 'dominant_fig.png'
        plt.savefig(file, dpi=400)
        plt.close(fig)
    return




def plot_background(list_dir, grid_lims):
    depth_file = list_dir['save_folder'] + list_dir['sim_folder'] + '/depth.npy'
    depth = np.load(depth_file)
    depth[np.isnan(depth)] = -35000
    fig, ax = plt.subplots()
    cmap2 = plt.get_cmap('gray')
    land_plot = ax.contourf(depth, levels=[-40000, -20000], extend='both', cmap=cmap2)
    depth_plot = ax.contour(depth, np.linspace(0, 1500, 4), extend='both', colors='k', alpha=0.15)

    # add latlong:
    idxlimx = grid_lims['idlimx']
    idxlimy = grid_lims['idlimy']
    add_latlon(idxlimx, idxlimy)
    plt.ylim([idxlimy[0], idxlimy[1]])
    plt.xlim([idxlimx[0], idxlimx[1]])
    plt.grid(alpha=0.45)
    return fig, ax


def plot_depth(list_dir, grid_lims):
    depth_file = list_dir['save_folder'] + list_dir['sim_folder'] + '/depth.npy'
    depth = np.load(depth_file)
    depth[np.isnan(depth)] = -35000
    depth2 = np.load(depth_file)
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('Blues')
    cmap2 = plt.get_cmap('gray')
    cmax = 4400
    crange = 200

    land_plot = ax.contourf(depth, levels=[-40000, -20000], extend='both', cmap=cmap2)
    depth_plot = ax.contourf(depth2, levels=np.arange(0, cmax, crange), extend='both', cmap=cmap, alpha=1)
    depth_plot2 = ax.contour(depth, np.linspace(0, 2000, 6), extend='both', colors='k', alpha=0.15)

    # add latlong:
    idxlimx = grid_lims['idlimx']
    idxlimy = grid_lims['idlimy']
    add_latlon(idxlimx, idxlimy)
    plt.ylim([idxlimy[0], idxlimy[1]])
    plt.xlim([idxlimx[0], idxlimx[1]])
    plt.grid(alpha=0.45)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(depth_plot, cax=ax_cb)
    cbar.ax.set_ylabel('Depth (m)', loc='center', size=12.5, weight='bold')
    cbar.ax.tick_params(labelsize=9, rotation=0)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    return fig, ax


def add_latlon(idxlimx, idxlimy):
    from get_trajectory import geo2grid
    GridData2_IMAX = 780
    GridData2_JMAX = 825
    margv = 25
    for axfig in range(1, 3):
        if axfig == 1:
            a2 = np.floor(np.linspace(idxlimy[0] + margv, idxlimy[1] - margv, 5))
            a1 = np.floor(GridData2_IMAX / 2) * np.ones_like(a2)
        else:
            a1 = np.floor(np.linspace(idxlimx[0] + margv, idxlimx[1] - margv, 5))
            a2 = np.floor(GridData2_JMAX / 2) * np.ones_like(a1)
        # Convert grid coordinates to Lambert Conic coordinates
        f1, f2 = geo2grid(a1, a2, 'get_bl')

        v1 = np.floor(np.mean(np.abs(np.diff(f1))))
        v12 = (v1 * np.ones_like(f1))
        v12[0] = 0
        f1 = np.ceil(f1[0]) + np.cumsum(v12)

        v2 = np.floor(np.mean(np.abs(np.diff(f2))))
        v22 = (v2 * np.ones_like(f2))
        v22[0] = 0
        f2 = np.ceil(f2[0]) + np.cumsum(v22)

        # Generate tick labels
        tick_labels = []
        for i in range(len(f1)):
            if axfig == 1:
                vs = str(-1 * np.floor(f1[i]).astype(int))
                tick_labels.append(vs + '$^\circ$S')
            else:
                vs = str(-1 * np.floor(f2[i]).astype(int))
                tick_labels.append(vs + '$^\circ$W')
            # Clear ticks and labels

        # Set ticks and labels based on the axis
        if axfig == 1:
            plt.yticks([])
            plt.yticks(a2, tick_labels)
        else:
            plt.xticks([])
            plt.xticks(a1, tick_labels)
    return


