# master file for krill model analysis
from Krillmod.get_trajectory import *
from Krillmod.analyse_trajectory import *
from Krillmod.plot_trajectory import *

#Function that outputs trajectory.nc file into a dictionary storing arrays, cftime and integers
# file = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'

file = 'A:/Cian_sinmod/antkrill_sim/trajectory.nc'
store = store_traj(file)
# field = dom_path(store)

#vals = get_traj(store, 150)
#plot_grid(vals)

# x1 = vals['xi']
# y1 = vals['yi']
# idx = (x1 > 0)
# x2 = x1[idx]
# y2 = y1[idx]
# lat1, lon1 = geo2grid(x2, y2, 'get_bl')
#
# plot_geo(lat1, lon1)
#
#
# x2, y2 = geo2grid(lat1, lon1, 'get_xy')

