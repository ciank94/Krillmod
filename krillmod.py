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

# Save the array to a binary file
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#np.save(file, np.array(df))

# Load the saved array
#df = np.load(file)
#
#plot_grid(file)







# field = dom_path(store)


