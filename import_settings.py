from Krillmod.get_trajectory import *
from Krillmod.analyse_trajectory import *
from Krillmod.plot_trajectory import *
import sys

shp_dir = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/ssmu/'
tr_dir = 'A:/Cian_sinmod/meeso_sim/'
sv_dir = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/'


#Consider using functions with input: Number of files, index of files etc.
shp_file = shp_dir + 'ssmusPolygon.shp'
tr_file = tr_dir + 'sim_2016/trajectory.nc'
reg_file = sv_dir + 'regions.nc'
time_file = sv_dir + 'times.npy'
depth_file = sv_dir + 'depth.npy'


if not os.path.exists(shp_file):
    print('Error: Directory "' + shp_file + '" does not exist')
    sys.exit()


if not os.path.exists(tr_file):
    print('Error: Directory "' + tr_file + '" does not exist')
    sys.exit()


if not os.path.exists(time_file):
    get_times(tr_file, time_file)

time = np.load(time_file, allow_pickle=True)


if not os.path.exists(depth_file):
    get_depth(tr_file, depth_file)

depth = np.load(depth_file)


if not os.path.exists(reg_file):
    print('Note: Reformatting trajectory data for analysis')
    store_regions(shp_file, tr_file, reg_file)

# Consider using conditional statements for using either the output from a function or the data saved already
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#df = np.load(file)


#Function that outputs trajectory.nc file into a dictionary storing arrays, cftime and integers
# file = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'