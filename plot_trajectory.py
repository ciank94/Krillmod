import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from Krillmod.get_trajectory import geo2grid


def plot_geo(lat1, lon1):
    ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent([125, 150, 35, 63])
    #     ax.stock_img()
    #
    #     ax.add_feature(cfeature.LAND)  # If I comment this => all ok, but I need
    #     ax.add_feature(cfeature.LAKES)
    #     ax.add_feature(cfeature.RIVERS)
    ax.coastlines()
    #     ax.coastlines()
    ax.scatter(lon1, lat1, transform=ccrs.PlateCarree())
    return


def plot_dom_path(df):
    plt.pcolor(df)
    plt.show()
    return


def plot_grid(vals):

    x1 = vals['xi']
    y1 = vals['yi']
    idx = (x1 > 0)
    x2 = x1[idx]
    y2 = y1[idx]
    idxlimx = [min(x2[:]), max(x2[:])]
    idxlimy = [min(y2[:]), max(y2[:])]
    dp = vals['dp']

    plt.pcolor(dp)
    plt.colorbar()
    cmap = plt.get_cmap('seismic')  # coolwarm, jet
    plt.clim(0, 6000)
    plt.set_cmap(cmap)
    plt.scatter(x2, y2, s=1, c='k')
    add_latlon(idxlimx, idxlimy)
    plt.ylim([idxlimy[0], idxlimy[1]])
    plt.xlim([idxlimx[0], idxlimx[1]])
    file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/init_grid.png'
    plt.savefig(file)
    open(file)
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
                vs = str(-1 * np.floor(f1[i]))
                tick_labels.append(vs + 'S')
            else:
                vs = str(-1 * np.floor(f2[i]))
                tick_labels.append(vs + 'W')
            # Clear ticks and labels

        # Set ticks and labels based on the axis
        if axfig == 1:
            plt.yticks([])
            plt.yticks(a2, tick_labels)
        else:
            plt.xticks([])
            plt.xticks(a1, tick_labels)
    return
