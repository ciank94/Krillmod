# master file for krill model analysis
from Krillmod.get_trajectory import *
from Krillmod.analyse_trajectory import *
from Krillmod.plot_trajectory import *

#Function that outputs trajectory.nc file into a dictionary storing arrays, cftime and integers
# file = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'
shp_file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/ssmu/ssmusPolygon.shp'
file = 'A:/Cian_sinmod/meeso_sim/sim_2016/trajectory.nc'
#file = 'A:/Cian_sinmod/sim_2017/trajectory.nc'
#plot_geo(file)
#plot_geo(file)
store = store_traj(file)
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#df = np.load(file)
df = dom_path(store)
plot_dom_path(df, store)
sim_account(store)
# Save the array to a binary file
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#np.save(file, np.array(df))

# Load the saved array
#df = np.load(file)
#
#plot_grid(file)
import shapely.geometry
shape_v = read_ssmu(shp_file)
pol_idx = 0
poly1 = shape_v[pol_idx:pol_idx +1].geometry
idx = 0
store_t = single_traj(store, idx)
x = store_t['xp']
y = store_t['yp']
lat, lon = geo2grid(x, y, 'get_bl')
shp = np.shape(lat)
part_in = np.zeros(shp[0])
for i in range(0, shp[0]):
    lon_i = lon[i]
    lat_i = lat[i]
    point1 = shapely.geometry.Point(lon_i, lat_i)
    in_poly = point1.within(poly1)
    if in_poly[1]:
        part_in[i] = 1


lat1 = lat[part_in > 0]
lon1 = lon[part_in > 0]

# First index should be for
list_v = [{}, {}, {}, {}, {}]
list_v[0]['Polygon ' + str(pol_idx)] = [lon, lat]
list_v = ["R1", "R2"]
foo = {}
foo['R1'][0] = lat1

idx = 0







point1 = lon[0:10], lat[0:10]





# field = dom_path(store)


