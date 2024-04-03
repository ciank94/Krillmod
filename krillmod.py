# master file for krill model analysis
import matplotlib.pyplot as plt

from set_directory import Folder, FolderComp
from get_trajectory import Trajectory, Regional
from analyse_trajectory import Analyse, AnalyseComp
#from get_trajectory import store_regions
#from analyse_trajectory import lagrangian_analysis, ssmu_start, sim_account
from plot_trajectory import (Plots, plot_connectivity, plot_retention, plot_transit, plot_dom_paths,
                             plot_depth_profile, animate_transit, animate_dom_paths)

#todo: put different comparisons in their own module- get_sim_analysis, compare_sim, and so on- which take folders and
#todo: files as input and do basic analysis; also put network drives in a list for easier extraction e.g. ssmu_shape,
#todo: betzy folder etc.

# Comparisons between simulations- use another comparison class:
sim_folder1 = 'DVM'  # Simulation identifier
sim_folder2 = '2017'  # Second simulation identifier
#trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
t_folder1 = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
t_folder2 = 'A:/Cian_sinmod/sim_'  # trajectory folder- saga

# set data source 1
# set data source 2
# retrieve relevant data for dimensions to be compared
# calculated RMSE
# plot RMSE for each data point with plotting class (can do this later)

f = FolderComp(t_folder1, sim_folder1, t_folder2, sim_folder2)
f.set_folder('local')
f.exist()

df = AnalyseComp(f)
df.compare_pathways()
df.compare_retention()

p = Plots()
#p.plot_back(f)
p.plot_comp_paths(files=f)
p.plot_comp_retention(files=f)
breakpoint()

sim = ['DVM']  # Simulation identifier
#trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
for s in sim:
    # Define names of main folders
    sim_folder = s  # Name of simulation instance

    # Setup directory information for analysis
    f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
    f.set_folder('local')  # Sets folders based on remote or local node
    f.exist()  # Checks existence of relevant folders and files

    # Two classes for dealing with data:
    k_r = Regional(f)
    k_t = Trajectory(f)

    # Now the analysis
    df = Analyse(f)
    df.dom_paths(k_r)
    df.depth_profile(k_t)  # Note
    df.retention_times(k_r)
    df.transit_times('ALL', k_r)

    plot_depth_profile(f)
    plot_connectivity(f, k_r)
    plot_transit(f)
    plot_dom_paths(f, k_r)
    plot_retention(f)





##### Extra;
#import netCDF4 as nc
# import matplotlib.pyplot as plt
# import numpy as np
# dl_file = 'A:/Cian_sinmod/meeso_sim/sim_generic/DL_save.nc'
# dl = nc.Dataset(dl_file, "r", format="NETCDF4")
# l_vals = dl.variables['DL']
# d = 12
# plt.figure()
# light_sample = l_vals[d,:,:]
# light_sample[light_sample<0]=np.nan
# t = dl.variables['time'][d]
# title_x = ('time= ' + str(t[2]) + '.' + str(t[1]) + '.' + str(t[0]) + ' ' + str(t[3]) + ':' + str(t[4]) + ':' + str(t[5]))# +
#             #  '     ' + str('layer') + '= ' + str(l))
# plt.contourf(light_sample)
# plt.title(title_x)
# plt.colorbar()
# plt.savefig('C:/Users/ciank/PycharmProjects/sinmod/Krillmod/dl.png',dpi=400)
# plt.show()
#
# # #
# #
# d = np.arange(24,48,1).astype(int)
# light_sample = l_vals[d,100,300]
# t = np.array(dl.variables['time'][d,:])
# plt.plot(t[:,3], light_sample, '-')
# #
# #
# # plt.show()
# # plt.savefig('C:/Users/ciank/PycharmProjects/sinmod/Krillmod/dl.png',dpi=400)
# light_sample =np.zeros([24])
# for i in range(0,24):
#     light_sample[i] = l_vals[i+1, l, 100, 300]


#import netCDF4 as nc
#tr_file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/sim_generic/trajectory.nc'
#tr_file= 'A:/Cian_sinmod/meeso_sim/sim_generic/trajectory.nc'
#tr = nc.Dataset(tr_file, "r", format="NETCDF4")

  # import matplotlib.pyplot as plt
    # import numpy as np
    # shp = np.shape(k_t.light)[0]
    # k_l = np.zeros(shp)
    # for i in range(0, shp):
    #     l_range = k_t.light[i, 0:4500]
    #     k_l[i] = np.nanmean(l_range[l_range>0])
    #
    # plt.scatter(k_l)
    # breakpoint()


    # Some analysis