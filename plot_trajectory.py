# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
import sys
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import os
# from analyse_trajectory import ssmu_target
from matplotlib.animation import PillowWriter
import matplotlib.cm as cm
import matplotlib.animation as animation


class Plots:

    def __init__(self):
        # Land plot parameters
        self.land = None
        self.land_par = [-40000, -20000]
        self.depth_contours = np.linspace(0, 1500, 4)
        self.depth_colors = np.arange(0, 4500, 200)
        # todo: make different standard areas for plotting;
        self.xlim_standard = [50, 750]
        self.ylim_standard = [50, 750]

        # Define colormaps
        self.land_cmap = plt.get_cmap('gray')
        self.error_cmap = plt.get_cmap('bwr')
        self.depth_cmap = plt.get_cmap('Blues')  # Alternative = cmocean.cm.deep

    def plot_back(self, files):
        save_name = 'back'
        fig, ax = plt.subplots()
        self.plot_background(files, ax)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        save_plot(files, save_name)
        return

    def plot_comp_paths(self, files):
        save_name = 'comp_paths'
        fig, ax = plt.subplots()
        self.plot_background(files, ax)
        data = np.load(files.comp_folder + 'comp_paths.npy')
        range_v = np.arange(-0.08, 0.081, 0.002)
        range_v = np.round(range_v, 3)
        comp_plot = ax.contourf(data, levels=range_v, cmap=self.error_cmap)
        self.plot_standard_lims()
        add_cbar(comp_plot, fig, ax, axis_title='%')
        save_comp_plot(files, save_name)
        return

    def plot_comp_retention(self, files):
        save_name = 'comp_retention'
        fig, ax = plt.subplots()
        self.plot_background(files, ax)
        data = np.load(files.comp_folder + 'comp_retention.npy')
        x_sites = data[0, :]
        y_sites = data[1, :]
        c_list = data[2, :]
        ax.scatter(x_sites, y_sites, s=10, facecolor='none', edgecolors='gray', alpha=0.3, linewidth=0.2)
        comp_plot = ax.scatter(x_sites, y_sites, s=10, c=c_list, cmap=self.error_cmap,
                               edgecolors='gray', vmin=np.nanmean(c_list[c_list < 0]) * 2.4,
                               vmax=np.nanmean(c_list[c_list > 0]) * 2.4, linewidth=0.2)
        self.plot_site_lims(files, x_sites, y_sites)
        add_cbar(comp_plot, fig, ax, axis_title='hours')
        save_comp_plot(files, save_name)
        return

    def plot_background(self, files, ax):
        # Assumes figure is already created
        self.get_land(files)
        ax.contourf(self.land, self.land_par, cmap=self.land_cmap)
        ax.contour(self.land, self.depth_contours, colors='k', alpha=0.15)
        return

    def get_land(self, files):
        depth_file = files.save + files.sim + '/depth.npy'
        depth = np.load(depth_file)
        self.land = np.array(depth[:])
        self.land[np.isnan(self.land)] = -35000
        return

    def plot_site_lims(self, files, x_sites, y_sites):
        depth_file = files.save + files.sim + '/depth.npy'
        depth = np.load(depth_file)
        self.xlim_standard = [np.nanmax([np.nanmin(x_sites) - 20, 0]),
                              np.nanmin([np.nanmax(x_sites) + 20, np.shape(depth)[1]])]
        self.ylim_standard = [np.nanmax([np.nanmin(y_sites) - 20, 0]),
                              np.nanmin([np.nanmax(y_sites) + 20, np.shape(depth)[0]])]
        add_latlon(self.xlim_standard, self.ylim_standard)
        plt.ylim(self.ylim_standard[0], self.ylim_standard[1])
        plt.xlim(self.xlim_standard[0], self.xlim_standard[1])
        plt.grid(alpha=0.45)
        return

    def plot_standard_lims(self):
        add_latlon(self.xlim_standard, self.ylim_standard)
        plt.ylim(self.ylim_standard[0], self.ylim_standard[1])
        plt.xlim(self.xlim_standard[0], self.xlim_standard[1])
        plt.grid(alpha=0.45)
        return

    def plot_depth(self, files, data):
        save_name = 'depth'
        fig, ax = plt.subplots()
        depth_1 = data.depth[:]
        depth_plot = ax.contourf(depth_1, levels=self.depth_colors, extend='both', cmap=self.depth_cmap, alpha=1)
        self.plot_background(files, ax)
        cbar_title = 'Depth ' + '(m)'
        add_cbar(depth_plot, fig, ax, cbar_title)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        save_plot(files, save_name)
        return


def plot_depth_profile(f):
    filepath = f.save + f.sim + '/depth_profile.npy'
    savepath = f.save + f.sim + '/depth_profile.png'
    z = np.load(filepath)
    fig, ax = plt.subplots()
    f1 = ax.contourf(z)

    ax1 = f1.axes
    ax1.invert_yaxis()
    f1.set_cmap('hot')
    ax.set_xlabel('Time step', fontsize=13)
    ax.set_ylabel('Depth (m)', fontsize=13)
    plt.colorbar(f1)
    plt.savefig(savepath, dpi=400)
    plt.close(fig)
    return


def plot_connectivity(f, k_r):
    key_list = ['SO', 'SG']
    for target_area in key_list:
        filepath = f.save + f.sim + '/' + target_area + '_transit.npy'
        save_path = f.save + f.sim + '/' + target_area + '_' + 'connect.png'
        if not os.path.exists(filepath):
            print('Error: Directory ' + filepath + ' with intermediate matrix does not exist, exiting')
            sys.exit()
        else:
            df = np.load(filepath)
            x = df[:, 0]
            y = df[:, 1]
            x1 = x[x > 0]
            y1 = y[y > 0]
            df2 = np.array([x1, y1]).T
            df3 = np.unique(df2, axis=0)
            xlist = df3[:, 0]
            ylist = df3[:, 1]
            tot_part = np.zeros(np.shape(df3)[0])
            if np.shape(tot_part)[0] == 0:
                print('No arrivals for connectivitiy plot, skipping')
            else:
                x = k_r.x[:, 0]
                y = k_r.y[:, 0]
                dfn = np.array([x, y]).T
                dfn2 = np.unique(dfn, axis=0)
                xl1 = dfn2[:, 0]
                yl1 = dfn2[:, 1]
                xl1 = xl1[xl1 > 0]
                yl1 = yl1[yl1 > 0]
                for i in range(0, np.shape(df3)[0]):
                    xv = xlist[i]
                    yv = ylist[i]
                    b1 = np.isin(x1, xv)
                    b2 = np.isin(y1, yv)
                    sumb = b1 * b2
                    tot_part[i] = np.sum(sumb)
                depth = k_r.depth
                grid_lims = dict()
                grid_lims['idlimx'] = [np.nanmax([np.nanmin(xl1) - 20, 0]),
                                       np.nanmin([np.nanmax(xl1) + 20, np.shape(depth)[1]])]
                grid_lims['idlimy'] = [np.nanmax([np.nanmin(yl1) - 20, 0]),
                                       np.nanmin([np.nanmax(yl1) + 20, np.shape(depth)[0]])]
                fig, ax = plot_background(f, grid_lims)
                ax.scatter(xl1, yl1, s=10, c='w', edgecolors='gray', alpha=0.3, linewidth=0.2)
                p2 = ax.scatter(xlist, ylist, s=10, c=tot_part, cmap='jet', edgecolors='gray', linewidth=0.2,
                                vmin=np.nanmin(tot_part), vmax=np.nanmax(tot_part))
                fig.colorbar(p2, label='Number of particles')
                # plt.legend(loc='upper center', bbox_to_anchor=(1.23, 1), ncol=1, fancybox=True, shadow=True, title='s2')
                plt.savefig(save_path, dpi=400)
                plt.close(fig)
    return


def plot_transit(f):
    # Arguments:
    # Pointer to subset of individuals (sub_key)
    # Use that pointer to identify folder
    # Use target function to identify target area (ssmu_target in this case)- loop over target areas (subs dictionary)
    key_list = ['SO', 'SG']
    for target_area in key_list:
        filepath = f.save + f.sim + '/' + target_area + '_transit.npy'
        save_path = f.save + f.sim + '/' + target_area + '_' + 'transit.png'
        save_path2 = f.save + f.sim + '/' + target_area + '_' + 'transit_map.png'
        if not os.path.exists(filepath):
            print('Error: Directory ' + filepath + ' with intermediate matrix does not exist, exiting')
            sys.exit()
        else:
            df = np.load(filepath)
            x = df[:, 0]
            y = df[:, 1]
            t = df[:, 2]
            # Calculate the fraction that reach SOI, the mean time, std, etc., plot histogram;
            x_arrive = x[x > 0.1]
            y_arrive = y[y > 0.1]
            t_arrive = t[t > 0.1]
            t_arrive = t_arrive / 24
            if np.shape(t_arrive)[0] == 0:
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
    return


def plot_retention(f):
    filepath = f.save + f.sim + '/retention.npy'
    save_path = f.save + f.sim + '/retention.png'
    if not os.path.exists(filepath):
        print('Error: Directory ' + filepath + ' with intermediate matrix does not exist, exiting')
        sys.exit()
    else:
        df = np.load(filepath)
        xlist = df[:, 0]
        ylist = df[:, 1]
        clist = df[:, 2]
        depth = np.load(f.depth_file)
        grid_lims = dict()
        grid_lims['idlimx'] = [np.nanmax([np.nanmin(xlist) - 20, 0]),
                               np.nanmin([np.nanmax(xlist) + 20, np.shape(depth)[1]])]
        grid_lims['idlimy'] = [np.nanmax([np.nanmin(ylist) - 20, 0]),
                               np.nanmin([np.nanmax(ylist) + 20, np.shape(depth)[0]])]
        fig, ax = plot_background(f, grid_lims)
        ax.scatter(xlist, ylist, s=10, facecolor='none', edgecolors='gray', alpha=0.3, linewidth=0.2)
        lim = 0.25 * np.nanmean(clist)
        p2 = ax.scatter(xlist[clist > lim], ylist[clist > lim], s=10, c=clist[clist > lim], cmap='YlOrBr',
                        edgecolors='gray', vmin=np.nanmin(clist),
                        vmax=np.nanmean(clist) * 1.75, linewidth=0.2)
        fig.colorbar(p2, label='hours')
        # plt.legend(loc='upper center', bbox_to_anchor=(1.23, 1), ncol=1, fancybox=True, shadow=True, title='s2')
        plt.savefig(save_path, dpi=400)
        plt.close(fig)
    return


def animate_transit(list_dir, sub_idx):
    key_list = list(sub_idx.keys())
    # for sub_key in key_list:
    for sub_key in ['WAP']:
        # subs = ssmu_target(sub_key)  # Target regions stored in dictionary (there may be multiple target areas);
        # for target_area in subs:
        for target_area in ['SO']:
            filepath = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit.npy'
            save_path = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit.gif'
            if not os.path.exists(filepath):
                print('Error: Directory ' + filepath + ' with intermediate matrix does not exist, exiting')
                sys.exit()
            else:
                df = np.load(filepath)
                x = df[:, 0]
                y = df[:, 1]
                t = df[:, 2]
                # Calculate the fraction that reach SOI, the mean time, std, etc., plot histogram;
                t_arrive = t[t > 0.1]
                t_arrive = t_arrive / 24
                if np.shape(t_arrive)[0] == 0:
                    print('No particles arrived, skipping')
                else:
                    filepath2 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormx.npy'
                    filepath3 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormy.npy'
                    wormx = np.load(filepath2)
                    wormy = np.load(filepath3)
                    sub_id = np.sum(wormx, axis=1) > 0
                    x_end = df[:, 3]
                    y_end = df[:, 4]
                    wx = wormx  #
                    wy = wormy
                    wxp = wx
                    wyp = wy
                    wxp[wxp == 0] = np.nan
                    wyp[wyp == 0] = np.nan

                    grid_lims = dict()
                    depth = np.load(list_dir['depth_file'])
                    grid_lims['idlimx'] = [np.nanmax([np.nanmin(wxp) - 50, 0]),
                                           np.nanmin([np.nanmax(wxp) + 20, np.shape(depth)[0]])]
                    grid_lims['idlimy'] = [np.nanmax([np.nanmin(wyp) - 100, 0]),
                                           np.nanmin([np.nanmax(wyp) + 50, np.shape(depth)[1]])]
                    fig, ax = plot_depth(list_dir, grid_lims)

                    artists = []
                    max_steps = np.floor(np.shape(wxp)[1] / 5).astype(int)
                    for i in range(0, max_steps):
                        container = ax.plot(wxp[:, i * 5], wyp[:, i * 5], 'r.', markersize=1.5, linewidth=0.8)
                        artists.append(container)

                    ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=2)
                    ani.save(filename=save_path, dpi=400, writer=PillowWriter(fps=30))
    return


def animate_dom_paths(list_dir, sub_idx):
    depth = np.load(list_dir['depth_file'])
    cmap1 = plt.get_cmap('OrRd')  # Oranges, Reds- sequential coolwarm= divergent, jet, seismic
    cmax = 100
    crange = 0.01
    sub_key = 'WAP'
    idv = (sub_idx[sub_key])
    save_path = list_dir[sub_key + '_folder'] + 'dom_paths.gif'
    nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
    df = np.zeros(np.shape(depth))
    df[np.isnan(depth)] = np.nan
    list_ids = np.where([depth > 0])

    grid_lims = dict()
    grid_lims['idlimx'] = [np.nanmax([np.nanmin(list_ids[2]) - 50, 0]),
                           np.nanmin([np.nanmax(list_ids[2]) + 20, np.shape(depth)[1]])]
    grid_lims['idlimy'] = [np.nanmax([np.nanmin(list_ids[1]) - 100, 0]),
                           np.nanmin([np.nanmax(list_ids[1]) + 50, np.shape(depth)[0]])]
    fig, ax = plot_background(list_dir, grid_lims)
    x = nc_file['xp'][idv, :]
    y = nc_file['yp'][idv, :]
    inc = 80
    artists = []
    tot_v = np.floor(np.shape(x)[1] / inc).astype(int)
    for k in range(0, tot_v):
        df = np.zeros(np.shape(depth))
        df[np.isnan(depth)] = np.nan
        # v1 = np.min([k + 10, np.shape(x)[1]])
        for i in range(0, np.shape(x)[0]):
            yi = y[i, (k * inc) + 1:((k + 1) * inc)]
            xi = x[i, (k * inc) + 1:((k + 1) * inc)]
            df[yi, xi] = df[yi, xi] + 1
        df[df > 0] = ((df[df > 0]) / np.nanmax(df[df > 0])) * 100
        cs = ax.contourf(df, levels=np.arange(0, cmax, crange), extend='both', cmap=cmap1)
        # divider = make_axes_locatable(ax)
        # ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
        # fig.add_axes(ax_cb)
        # cbar = plt.colorbar(cs, cax=ax_cb)
        # cbar.ax.set_ylabel('Probability (%)', loc='center', size=12.5, weight='bold')
        # cbar.ax.tick_params(labelsize=9, rotation=0)
        # ax.tick_params(axis='x', labelsize=12)
        # ax.tick_params(axis='y', labelsize=12)
        artists.append(cs.collections)
        cs.remove()
        plt.show()
    ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=2)
    ani.save(filename=save_path, dpi=400, writer=PillowWriter(fps=30))
    return


def plot_dom_paths(f, k):
    depth = k.depth
    dom_file = f.save + f.sim + '/dominant_paths.npy'
    df = np.load(dom_file)
    df[np.isnan(depth)] = np.nan
    list_ids = np.where([df > 0])

    grid_lims = dict()
    grid_lims['idlimx'] = [np.nanmax([np.nanmin(list_ids[2]) - 50, 0]),
                           np.nanmin([np.nanmax(list_ids[2]) + 20, np.shape(depth)[1]])]
    grid_lims['idlimy'] = [np.nanmax([np.nanmin(list_ids[1]) - 100, 0]),
                           np.nanmin([np.nanmax(list_ids[1]) + 50, np.shape(depth)[0]])]
    fig, ax = plot_background(f, grid_lims)
    # breakpoint()
    # NOTE: Do a check for the existence of file
    # Modify for the limits of the plot- hard coding elsewhere- create a dictionary with options for plotting
    # tot_file = sv_dir + tr_folder + '/' + area_name + '_total_visits.npy'

    cmax = 0.3
    crange = 0.02

    # Dominant pathways
    cmap1 = plt.get_cmap('OrRd')  # Oranges, Reds- sequential coolwarm= divergent, jet, seismic
    data_plot = ax.contourf(df, levels=np.arange(0, cmax, crange), extend='both', cmap=cmap1)
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(data_plot, cax=ax_cb)
    cbar.ax.set_ylabel('Probability (%)', loc='center', size=12.5, weight='bold')
    cbar.ax.tick_params(labelsize=9, rotation=0)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    file = f.save + f.sim + '/dominant_fig.png'
    plt.savefig(file, dpi=400)
    plt.close(fig)
    return


def plot_background(f, grid_lims):
    depth_file = f.save + f.sim + '/depth.npy'
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


def add_cbar(plot, fig, ax, axis_title):
    # plt.grid(alpha=0.45)
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(plot, cax=ax_cb)
    cbar.ax.set_ylabel(axis_title, loc='center', size=7, weight='bold')
    cbar.ax.tick_params(labelsize=7, rotation=0)
    return


def save_plot(files, save_name):
    savefile = files.save + files.sim + '/' + save_name + '.png'
    print('Saving file: ' + savefile)
    plt.savefig(savefile, dpi=400)
    plt.close()
    return


def save_comp_plot(files, save_name):
    savefile = files.comp_folder + save_name + '.png'
    print('Saving file: ' + savefile)
    plt.savefig(savefile, dpi=400)
    plt.close()
    return

# def plot_regions(list_dir):
#     import matplotlib.pyplot as plt
#     from plot_trajectory import plot_background
#     import numpy as np
#
#     filepath = list_dir['shape_folder'] + str('poly2grid.npy')
#     depth = np.load(list_dir['depth_file'])
#     df = np.load(filepath)
#     df[np.isnan(depth)] = np.nan
#     df[df < 0] = -100
#     list_ids = np.where([df > 0])
#     grid_lims = dict()
#     grid_lims['idlimx'] = [np.nanmax([np.nanmin(list_ids[2]) - 50, 0]),
#                            np.nanmin([np.nanmax(list_ids[2]) + 20, np.shape(depth)[1]])]
#     grid_lims['idlimy'] = [np.nanmax([np.nanmin(list_ids[1]) - 100, 0]),
#                            np.nanmin([np.nanmax(list_ids[1]) + 50, np.shape(depth)[0]])]
#     fig, ax = plot_background(list_dir, grid_lims)
#     cmap1 = plt.get_cmap('coolwarm')  # Oranges, Reds- sequential coolwarm= divergent, jet, seismic
#     ax.contourf(df, levels=np.arange(0, 1, 0.1), extend='both', cmap=cmap1)
#     plt.savefig(list_dir['shape_folder'] + 'regions.png', dpi=400)
#     plt.close(fig)
#     return
