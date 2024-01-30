from Krillmod.get_trajectory import *
from Krillmod.analyse_trajectory import *
from Krillmod.plot_trajectory import *

shp_dir = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/ssmu/'
tr_dir = 'A:/Cian_sinmod/meeso_sim/'
sv_dir = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/'


#Consider using functions with input: Number of files, index of files etc.
shp_file = shp_dir + 'ssmusPolygon.shp'
tr_file = tr_dir + 'sim_2016/trajectory.nc'

# Consider using conditional statements for using either the output from a function or the data saved already
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#df = np.load(file)


#Function that outputs trajectory.nc file into a dictionary storing arrays, cftime and integers
# file = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'