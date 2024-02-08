# master file for krill model analysis
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Krillmod.import_settings import *

# Trajectory data reformatted:
nc_file = Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
#
start_id = nc_file['start'][:]
start_b = nc_file['act_part'][:]
uniq_b = np.unique(start_b)
shp_b = np.shape(uniq_b)[0]
#
subs_so = (start_id == 9) | (start_id == 10) | (start_id == 11)
subs_wap = (start_id == 2) | (start_id == 3) | (start_id == 4) | (start_id == 5) | (start_id == 6) | (start_id == 7)
area_idx = dict([('SO', subs_so), ('WAP', subs_wap)])

keysList = list(area_idx.keys())

ar_idv = 0
area_name = keysList[ar_idv]

idv = (area_idx[area_name])
#
#
#
df = np.zeros(np.shape(depth))
df[np.isnan(depth)] = np.nan
df2 = np.zeros(np.shape(depth))
df2[np.isnan(depth)] = np.nan

x = nc_file['xp'][idv, :].astype(int)
y = nc_file['yp'][idv, :].astype(int)
array_shape = np.shape(depth)
for i in range(0, np.shape(x)[0]):
    yi = y[i, :]
    xi = x[i, :]
    df[yi, xi] = df[yi, xi] + 1
    for j in range(0, np.shape(xi)[0]):
        df2[yi[j], xi[j]] = df2[yi[j],xi[j]] + 1

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

#Worm plot: Show the most important paths somehow, also could display start points on colourmaps;
#Dom. pathways
depth[np.isnan(depth)] = -35000
df[depth < 0] = np.nan
df2[depth<0] = np.nan
if area_name == 'WAP':
    corr_max = 30
    cmax_dom = np.max(df[df > 0])/1.2
    cmax_vis = np.max(df2[df2 > 0]) / 10
if area_name == 'SO':
    corr_max = 100
    cmax_dom = np.max(df[df > 0])/3
    cmax_vis = np.max(df2[df2 > 0]) / 6


idxlimx = [np.min(x[x>0]), np.max(x[x>0]) - corr_max]
idxlimy = [np.min(y[y>0]), np.max(y[y>0]) - corr_max]
#plt.figure()
fig, ax = plt.subplots()
#Total number of visits

cmap1 = plt.get_cmap('Oranges')  # coolwarm= divergent, jet, seismic
cmap2 = plt.get_cmap('gray')
landplot = ax.contourf(depth,levels= [-35000,-20000],extend='both',cmap=cmap2)
#landplot = ax.pcolor(depth,cmap=cmap2,clim=[-35000,-20000])

#testplot = ax.contourf(df,cmap=cmap1,clim=[0,3000])

testplot = ax.contourf(df,levels=np.linspace(0,cmax_dom,15),extend='both',cmap=cmap1)

add_latlon(idxlimx, idxlimy)
plt.ylim([idxlimy[0], idxlimy[1]])
plt.xlim([idxlimx[0], idxlimx[1]])


divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cbar = plt.colorbar(testplot,cax=ax_cb)
cbar.ax.tick_params(labelsize=8)
#plt.set_cmap('Oranges')



#clim_dom = [0, cmax_dom]
# plt.clim(clim_dom[0], clim_dom[1])
# plt.set_cmap(cmap)
# plt.colorbar()

file = sv_dir + keysList[ar_idv] + '_' + 'dom_paths.png'
plt.savefig(file)



#Output = number of unique individuals in each cell & number of visits;
# To do: use poly2grid.npy indices for SO to count number of visits to each cell;
# and put the plotting stuff below into a functions for 3 separate plots;
fig, ax = plt.subplots()
#Total number of visits
depth[np.isnan(depth)] = -35000
cmap1 = plt.get_cmap('Blues')  # coolwarm= divergent, jet, seismic
cmap2 = plt.get_cmap('gray')
landplot = ax.contourf(depth,levels= [-35000,-20000],extend='both',cmap=cmap2)
#landplot = ax.pcolor(depth,cmap=cmap2,clim=[-35000,-20000])

#testplot = ax.pcolor(df,cmap=cmap1,clim=[0, 3000])

testplot = ax.contourf(df2,levels=np.linspace(0,cmax_vis,15),extend='both',cmap=cmap1)

add_latlon(idxlimx, idxlimy)
plt.ylim([idxlimy[0], idxlimy[1]])
plt.xlim([idxlimx[0], idxlimx[1]])


divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="3%", pad=0.05, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cbar = plt.colorbar(testplot,cax=ax_cb)
cbar.ax.tick_params(labelsize=8)


file = sv_dir + keysList[ar_idv] + '_' + 'total_visits.png'
plt.savefig(file)

# Visits per particle- Retention regions;
df2[df2==0] = 0.000001
df[df==0] = 0.000001
df3 = np.divide(df2, df)

plt.figure()

plt.pcolor(df3)
cmap = plt.get_cmap('coolwarm')  # coolwarm, jet, seismic
if area_name == 'WAP':
    cmax_ret = np.max(df3[df3 > 0]) / 50
elif area_name == 'SO':
    cmax_ret = np.max(df3[df3 > 0]) / 20

clim_ret = [0, cmax_ret]
plt.clim(clim_ret[0], clim_ret[1])
plt.set_cmap(cmap)
plt.colorbar()
add_latlon(idxlimx, idxlimy)
plt.ylim([idxlimy[0], idxlimy[1]])
plt.xlim([idxlimx[0], idxlimx[1]])
file = sv_dir + keysList[ar_idv] + '_' + 'retention.png'
plt.savefig(file)


lonW=-70
lonE=-20
latS=-70
latN=-50
coordinates = (lonW, lonE, latS, latN)
coordv = np.where(df3>0)
x2 = coordv[1]
y2 = coordv[0]
lat1, lon1 = geo2grid(x2, y2, 'get_bl')
lon_bins_2d, lat_bins_2d = np.meshgrid(lon1, lat1)

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
plt.pcolormesh(x, y, df3[df3>0])
plt.colorbar(orientation='horizontal')

m.scatter(x, y, 3, marker='o', color='r')


#Occupancy rate plots:

# Worm plots
plt.scatter(st_x, st_y, 1, 'b')
