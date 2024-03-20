# master file for krill model analysis
from set_directory import Folder
from get_trajectory import Trajectory, Regional
from analyse_trajectory import Analyse
#from get_trajectory import store_regions
#from analyse_trajectory import lagrangian_analysis, ssmu_start, sim_account
from plot_trajectory import (plot_connectivity, plot_retention, plot_transit, plot_dom_paths,
                             plot_depth_profile, animate_transit, animate_dom_paths)

name = ['samp_part.nc']
for n in name:
    # Define names of main folders
    trj_folder = 'B:/'  # trajectory folder
    sim_name = n # Name of simulation instance
    tr_file = trj_folder + sim_name

    # Setup directory information for analysis
    f = Folder(tr_file)  # Initialize folders
    f.set_folder('local')  # Sets folders based on remote or local node
    f.exist()  # Checks existence of relevant folders and files

    # Two classes for dealing with data:
    k_r = Regional(f)
    k_t = Trajectory(f)

    # Now the analysis
    df = Analyse(f)
    df.dom_paths(k_r)
    df.depth_profile(k_t) #Note
    df.retention_times(k_r)
    # df.transit_times('ALL', k_r)
    #
    plot_depth_profile(f)
    # plot_connectivity(f, k_r)
    # plot_transit(f)
    plot_dom_paths(f, k_r)
    plot_retention(f)


#store_regions(f)

#
# # Analyse trajectories:
# sub_idx = ssmu_start(list_dir['reg_file'])
# list_dir = lagrangian_analysis(comp_node, list_dir, sub_idx)

# plot trajectories
#plot_connectivity(f, k_r)
#plot_transit(f)
#plot_dom_paths(f, k_r)
#plot_retention(list_dir, sub_idx)

#animate_transit(list_dir, sub_idx)
#animate_dom_paths(list_dir, sub_idx)





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