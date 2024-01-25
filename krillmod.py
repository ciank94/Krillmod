# master file for krill model analysis
from Krillmod.get_trajectory import *
from Krillmod.analyse_trajectory import *
from Krillmod.plot_trajectory import *

#Function that outputs trajectory.nc file into a dictionary storing arrays, cftime and integers
# file = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'

file = 'A:/Cian_sinmod/antkrill_sim/trajectory.nc'
plot_geo(file)
#store = store_traj(file)
# df = dom_path(store)
# plot_dom_path(df)

#
plot_grid(file)







# field = dom_path(store)


