import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
#from Krillmod.get_trajectory import geo2grid, store_traj
#from Krillmod.import_settings import tr_file
from io import BytesIO

def plot_geo(file):
    #store = store_traj(file)

    vals = store_traj(file, 0)
    # plot_grid(vals)

    x1 = vals['xi']
    y1 = vals['yi']
    idx = (x1 > 0)
    x2 = x1[idx]
    y2 = y1[idx]
    lat1, lon1 = geo2grid(x2, y2, 'get_bl')
    #fig = plt.figure(figsize=(8, 6))
    lonW = -85
    lonE = -20 #lonE=-29
    latS = -70 #latS=-74
    latN = -50 #latN = -54
    coordinates = (lonW, lonE, latS, latN)
    res = 'f'
    m = Basemap(projection='merc',
                llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                resolution=res)
    m.drawcoastlines(linewidth=.5)
    m.drawmeridians(np.arange(-180., 180., 10.), labels=[False, False, False, True])
    m.drawparallels(np.arange(-90., 90., 4.), labels=[True, False, False, False])
    #m.etopo()
    m.shadedrelief()
    x, y = m(lon1, lat1)
    m.scatter(x, y, 3, marker='o', color='r')
    file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/init_geo.png'
    plt.savefig(file)
    # ax.set_extent([125, 150, 35, 63])
    #     ax.stock_img()
    #
    #     ax.add_feature(cfeature.LAND)  # If I comment this => all ok, but I need
    #     ax.add_feature(cfeature.LAKES)
    #     ax.add_feature(cfeature.RIVERS)
    #ax.coastlines(resolution='10m')
    # coast = cfeature.GSHHSFeature(scale='coarse', levels=[1])
    # ax.add_feature(coast)
    #     ax.coastlines()
    return


def plot_dom_path(df, store):
    shp = np.shape(store['xp'])
    dv = shp[1] #Note, this equals total particles for the entire simulation- results are only valid for entire frame
    df2 = (df/dv)*100
    plt.pcolor(df2)
    y2, x2 = np.where(df > 0)
    fv1 = 1
    fv2 = 1
    idxlimx = [min(x2[:]) - fv1, max(x2[:]) + fv2]
    idxlimy = [min(y2[:]) - fv1, max(y2[:]) + fv2]
    cmap = plt.get_cmap('coolwarm')  # coolwarm, jet, seismic
    #plt.clim(0, 2)
    plt.set_cmap(cmap)
    plt.colorbar()
    add_latlon(idxlimx, idxlimy)
    plt.ylim([idxlimy[0], idxlimy[1]])
    plt.xlim([idxlimx[0], idxlimx[1]])
    file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.png'
    plt.savefig(file)
    #plt.show()
    return


def plot_grid(file):


    store = store_traj(tr_file, 0)
    vals = store_traj(store, 0)

    x1 = vals['xi']
    y1 = vals['yi']
    idx = (x1 > 0)
    x2 = x1[idx]
    y2 = y1[idx]
    fv1 = 50
    fv2 = 100
    idxlimx = [min(x2[:])-fv1, max(x2[:])+fv2]
    idxlimy = [min(y2[:])-fv1, max(y2[:])+fv2]
    dp = store['depth']

    plt.pcolor(dp)
    plt.colorbar()
    cmap = plt.get_cmap('coolwarm')  # coolwarm, jet
    plt.clim(0, 6000)
    plt.set_cmap(cmap)
    plt.scatter(x2, y2, s=1, c='k')
    plt.scatter(x2[:, 0], y2[:, 0], s=5, c='w')
    #plt.scatter(x2[:, -1], y2[:, -1], s=2, c='b')
    add_latlon(idxlimx, idxlimy)
    plt.ylim([idxlimy[0], idxlimy[1]])
    plt.xlim([idxlimx[0], idxlimx[1]])
    file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/init_grid.png'
    plt.savefig(file)
    # open(file)
    # plt.ioff()
    # plt.show()
    return


def init_grid(store, stepv):
    from Krillmod.get_trajectory import geo2grid
    field = np.array(store['depth'])
    plt.figure()
    plt.pcolor(field)
    plt.colorbar()
    cmap = plt.get_cmap('seismic')  # coolwarm, jet
    plt.clim(0, 6000)
    plt.set_cmap(cmap)

    c = 0
    xp = np.empty(200000)
    yp = np.empty(200000)
    # stepv = 5
    for i in range(200, 550, stepv):
        for j in range(50, 450, stepv):
            if 1 < field[i, j] < 7000:
                c = c + 1
                xp[c] = i
                yp[c] = j

    print('Number of sites: ' + str(c))
    x1 = xp[xp > -1]
    y1 = yp[yp > -1]
    plt.scatter(y1, x1, s=1, c='k')
    plt.ylim([200, 600])
    plt.xlim([50, 500])
    file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/init_grid.png'
    plt.savefig(file)

    # geo_plot
    # lat, lon = geo2grid(y1, x1, case='get_bl')
    # cLon = -44.0
    # cLat = -45.0
    # lonW = -96
    # lonE = -28
    # latS = -78
    # latN = -57
    # projPC = ccrs.PlateCarree()
    # projps = ccrs.SouthPolarStereo()
    # fig = plt.figure(figsize=(11, 8.5))
    # ax = plt.subplot(1, 1, 1, projection=projps)
    # # ax.set_title('Lambert Conformal')
    #
    # # ax.set_extent([lonW, lonE, latS, latN], crs=projLcc)
    # ax.coastlines(resolution='110m', color='black')
    # gl = ax.gridlines(
    #     draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--'
    # )
    # ax.set_extent([-90, -20, -90, -50], crs=projps)
    # ax.coastlines(resolution='50m')
    # ax.scatter(lon, lat, c='b', marker='.', transform=projps)
    # file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/init_geo.png'
    # plt.savefig(file)
    return


def add_latlon(idxlimx, idxlimy):
    from Krillmod.get_trajectory import geo2grid
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


def plot_dom_paths(tr_folder, reg_file, sv_dir):
    depth_file = sv_dir + tr_folder + '/depth.npy'
    depth = np.load(depth_file)


    from Krillmod.get_trajectory import region_part
    area_idx = region_part(reg_file)
    key_list = list(area_idx.keys())
    for area_name in key_list:
        dom_file = sv_dir + tr_folder + '/' + area_name + '_dom_paths.npy'
        tot_file = sv_dir + tr_folder + '/' + area_name + '_total_visits.npy'
        for i in range(0, 2):
            if i == 1:
                cmax = 2
                crange = 0.1
            else:
                cmax = 0.01
                crange = 0.001

            if area_name == 'WAP':
                idxlimx = [100, 500]
                idxlimy = [250, 650]
                so_nparts = 45621
            if area_name == 'SO':
                idxlimx = [300, 600]
                idxlimy = [400, 700]
                so_nparts = 12221
            if area_name == 'APP':
                idxlimx = [100, 500]
                idxlimy = [250, 650]
                so_nparts = 122452
            if area_name == 'EAP':
                idxlimx = [100, 500]
                idxlimy = [250, 650]
                so_nparts = 17061
            if area_name == 'SOP':
                idxlimx = [300, 600]
                idxlimy = [400, 700]
                so_nparts = 66792
            #     cmax_dom = np.max(df[df > 0])/3
            #     cmax_vis = np.max(df2[df2 > 0]) / 6
            #Load dataframe from file
            if i == 1:
                df = np.load(dom_file)
                depth[np.isnan(depth)] = -35000
                df[depth < 0] = np.nan
                df[df>0] = (df[df>0]/so_nparts)*100

                fig, ax = plt.subplots()
                params = {'axes.labelsize': 12, 'axes.titlesize': 12, 'legend.fontsize': 12,
                          'xtick.labelsize': 12, 'ytick.labelsize': 12}
                plt.rcParams['font.sans-serif'] = "Times New Roman"
                #plt.rc('font', size=12)

                #Dominant pathways
                cmap1 = plt.get_cmap('OrRd')  # Oranges, Reds- sequential coolwarm= divergent, jet, seismic
                cmap2 = plt.get_cmap('gray')
                land_plot = ax.contourf(depth,levels= [-40000,-20000],extend='both',cmap=cmap2)
                depth_plot = ax.contour(depth, np.linspace(0, 1500, 4), extend='both', colors='k', alpha=0.15)
                data_plot = ax.contourf(df,levels=np.arange(0, cmax, crange),extend='both',cmap=cmap1)
                #
            #levels=np.linspace(0,cmax,10)
                add_latlon(idxlimx, idxlimy)
                plt.ylim([idxlimy[0], idxlimy[1]])
                plt.xlim([idxlimx[0], idxlimx[1]])
                plt.grid(alpha=0.45)

                divider = make_axes_locatable(ax)
                ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
                fig.add_axes(ax_cb)
                cbar = plt.colorbar(data_plot, cax=ax_cb)
                cbar.ax.set_ylabel('Probability (%)', loc='center', size=12.5, weight='bold')
                cbar.ax.tick_params(labelsize=9,rotation=0)
                ax.tick_params(axis='x', labelsize=12)
                ax.tick_params(axis='y', labelsize=12)
                #plt.rcParams.update(params)
                #plt.show()
                file = sv_dir + tr_folder + '/' + area_name + '_dom_paths.svg'
                plt.savefig(file)
            else:
                df = np.load(tot_file)
                depth[np.isnan(depth)] = -35000
                df[depth < 0] = np.nan
                df[df > 0] = (df[df > 0] / np.sum(df[df > 0])) * 100

                fig, ax = plt.subplots()
                params = {'axes.labelsize': 12, 'axes.titlesize': 12, 'legend.fontsize': 12,
                          'xtick.labelsize': 12, 'ytick.labelsize': 12}
                plt.rcParams['font.sans-serif'] = "Times New Roman"
                # plt.rc('font', size=12)

                # Dominant pathways
                cmap1 = plt.get_cmap('Blues')  # Oranges, Reds- sequential coolwarm= divergent, jet, seismic
                cmap2 = plt.get_cmap('gray')
                land_plot = ax.contourf(depth, levels=[-40000, -20000], extend='both', cmap=cmap2)
                depth_plot = ax.contour(depth, np.linspace(0, 1500, 4), extend='both', colors='k', alpha=0.15)
                data_plot = ax.contourf(df, levels=np.arange(0, cmax, crange), extend='both', cmap=cmap1)
                #
                # levels=np.linspace(0,cmax,10)
                add_latlon(idxlimx, idxlimy)
                plt.ylim([idxlimy[0], idxlimy[1]])
                plt.xlim([idxlimx[0], idxlimx[1]])
                plt.grid(alpha=0.45)

                divider = make_axes_locatable(ax)
                ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
                fig.add_axes(ax_cb)
                cbar = plt.colorbar(data_plot, cax=ax_cb)
                cbar.ax.set_ylabel('Probability (%)', loc='center', size=12.5, weight='bold')
                cbar.ax.tick_params(labelsize=9, rotation=0)
                ax.tick_params(axis='x', labelsize=12)
                ax.tick_params(axis='y', labelsize=12)
                # plt.rcParams.update(params)
                # plt.show()
                file = sv_dir + tr_folder + '/' + area_name + '_tot_visits.svg'
                plt.savefig(file)
            plt.close(fig)

    return fig


def get_frame_as_image(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    return plt.imread(buf)


def animate_frame(image):
    fig, ax = plt.subplots()
    ax.imshow(image)
    plt.show()


def plot_retention(sv_dir, tr_folder, time_file, reg_file):
    from Krillmod.get_trajectory import region_part
    area_idx = region_part(reg_file)
    key_list = list(area_idx.keys())
    for area_name in key_list:
        rfile = sv_dir + tr_folder + '/' + area_name + '_retention_days.npy'
        data = np.load(rfile)

        fig, ax = plt.subplots()
        list_id = np.unique(data[:, 1]).astype(int)
        v = np.zeros([len(list_id), 1])
        v2 = np.zeros([len(list_id), 1])
        time = np.load(time_file, allow_pickle=True)
        c = -1
        for i in range(0, np.shape(list_id)[0]):
            id1 = data[:, 1] == list_id[i]
            vals = data[id1, 0]
            v1 = vals[vals > -1]
            v[i, 0] = np.median(v1)
            v2[i, 0] = np.mean(v1)

        date_time = time[list_id]
        list_labels = []
        for i in range(0, np.shape(date_time)[0]):
            if i == 0:
                list_labels = [date_time[i].strftime("%d/%m/%Y")]
            else:
                list_labels.append(date_time[i].strftime("%d/%m/%Y"))

        ax.plot(list_labels, v)
        ax.plot(list_labels, v2)
        ax.tick_params(axis='x', labelsize=8, rotation=20)
        ax.tick_params(axis='y', labelsize=10)
        ax.xaxis.set_major_locator(plt.MaxNLocator(8))
        ax.legend(['median', 'mean'])
        ax.set_ylabel('Days in ' + area_name + ' region')
        file = sv_dir + tr_folder + '/' + area_name + '_retention.svg'
        plt.savefig(file)
    return

# start_p = 0
# stop_p = np.shape(x)[0]
# iter_p = np.shape(start_b[(start_b == uniq_b[0]) & (idv == True)])[0]
# xp1 = x[start_p:stop_p:iter_p, :]
# yp1 = y[start_p:stop_p:iter_p, :]
# plt.figure()
# # xp1[xp1 == 0] = np.empty
# # yp1[yp1 == 0] = np.empty
# plt.pcolor(depth)
# x1 = np.array(x)
# y1 = np.array(y)
# plt.scatter(x1, y1)
# plt.ylim([idxlimy[0], idxlimy[1]])
# plt.xlim([idxlimx[0], idxlimx[1]])
# start_id = nc_file['start'][:]
#     start_b = nc_file['act_part'][:]
#     uniq_b = np.unique(start_b)
#     shp_b = np.shape(uniq_b)[0]
#
# #Worm plot: Show the most important paths somehow, also could display start points on colourmaps;
# #Dom. pathways
# depth[np.isnan(depth)] = -35000
# df[depth < 0] = np.nan
# df2[depth<0] = np.nan
# if area_name == 'WAP':
#     corr_max = 30
#     cmax_dom = np.max(df[df > 0])/1.2
#     cmax_vis = np.max(df2[df2 > 0]) / 10
# if area_name == 'SO':
#     corr_max = 100
#     cmax_dom = np.max(df[df > 0])/3
#     cmax_vis = np.max(df2[df2 > 0]) / 6
#

# idxlimx = [np.min(x[x>0]), np.max(x[x>0]) - corr_max]
# idxlimy = [np.min(y[y>0]), np.max(y[y>0]) - corr_max]

#
#
#
# #clim_dom = [0, cmax_dom]
# # plt.clim(clim_dom[0], clim_dom[1])
# # plt.set_cmap(cmap)
# # plt.colorbar()
#
# file = sv_dir + keysList[ar_idv] + '_' + 'dom_paths.png'
# plt.savefig(file)
#
#
#
# #Output = number of unique individuals in each cell & number of visits;
# # To do: use poly2grid.npy indices for SO to count number of visits to each cell;
# # and put the plotting stuff below into a functions for 3 separate plots;
# fig, ax = plt.subplots()
# #Total number of visits
# depth[np.isnan(depth)] = -35000
# cmap1 = plt.get_cmap('Blues')  # coolwarm= divergent, jet, seismic
# cmap2 = plt.get_cmap('gray')
# landplot = ax.contourf(depth,levels= [-35000,-20000],extend='both',cmap=cmap2)
# #landplot = ax.pcolor(depth,cmap=cmap2,clim=[-35000,-20000])
#
# #testplot = ax.pcolor(df,cmap=cmap1,clim=[0, 3000])
#
# testplot = ax.contourf(df2,levels=np.linspace(0,cmax_vis,15),extend='both',cmap=cmap1)
#
# add_latlon(idxlimx, idxlimy)
# plt.ylim([idxlimy[0], idxlimy[1]])
# plt.xlim([idxlimx[0], idxlimx[1]])
#
#
# divider = make_axes_locatable(ax)
# ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
# fig.add_axes(ax_cb)
# cbar = plt.colorbar(testplot,cax=ax_cb)
# cbar.ax.tick_params(labelsize=8)
#
#
# file = sv_dir + keysList[ar_idv] + '_' + 'total_visits.png'
# plt.savefig(file)
#
#
#
#
# #Occupancy rate plots:
#
# # Worm plots
# plt.scatter(st_x, st_y, 1, 'b')