# master file for krill model analysis
from Krillmod.import_settings import *

# Trajectory data reformatted:
nc_file = Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
#
start_id = nc_file['start'][:]
start_b = nc_file['act_part'][:]
uniq_b = np.unique(start_b)
shp_b = np.shape(uniq_b)[0]
#
subs = (start_id == 9) | (start_id == 10) | (start_id == 11)
area_idx = dict([('SO', subs)])
# subs = (start_id == 2) | (start_id == 3) | (start_id == 4) | (start_id == 5) | (start_id == 6) | (start_id == 7)
#
#
df = np.zeros(np.shape(depth))
df[np.isnan(depth)] = np.nan
idv = (area_idx['SO'])
x = nc_file['xp'][idv, :].astype(int)
y = nc_file['yp'][idv, :].astype(int)
array_shape = np.shape(depth)
for i in range(0, np.shape(x)[0]):
    yi = y[i, :]
    xi = x[i, :]
    df[yi, xi] = df[yi, xi] + 1

#Output = number of unique individuals in each cell & number of visits;
# To do: use poly2grid.npy indices for SO to count number of visits to each cell;
# and put the plotting stuff below into a functions for 3 separate plots;
idxlimx = [250, 550]
idxlimy = [400, 650]
plt.figure()

#Dom. pathways
plt.pcolor(df)
cmap = plt.get_cmap('coolwarm')  # coolwarm, jet, seismic
plt.clim(0,7)
plt.set_cmap(cmap)
plt.colorbar()
plt.ylim([idxlimy[0], idxlimy[1]])
plt.xlim([idxlimx[0], idxlimx[1]])

#Occupancy rate plots:

# Worm plots
plt.scatter(st_x, st_y, 1, 'b')
