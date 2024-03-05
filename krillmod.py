# master file for krill model analysis
from import_settings import locate_folders
from analyse_trajectory import lagrangian_analysis, ssmu_start, sim_account
from plot_trajectory import (plot_connectivity, plot_retention, plot_transit, plot_dom_paths,
                             animate_transit, animate_dom_paths)

# # Define directories
#shape_name = 'mpasPolygon.shp'

#for sim_folder in ['2016', '2017', '2018', '2019']:
shape_name = 'ssmusPolygon.shp'
sim_folder = '2016'
comp_node = 'local'
list_dir = locate_folders(comp_node, sim_folder, shape_name)

# Simulation summary:
sim_account(list_dir)

# Analyse trajectories:
sub_idx = ssmu_start(list_dir['reg_file'])
list_dir = lagrangian_analysis(comp_node, list_dir, sub_idx)

# plot trajectories
plot_connectivity(list_dir, sub_idx)
#plot_transit(list_dir, sub_idx)
#plot_dom_paths(list_dir, sub_idx)
plot_retention(list_dir, sub_idx)

#animate_transit(list_dir, sub_idx)
#animate_dom_paths(list_dir, sub_idx)





##### Extra;
# import netCDF4 as nc
# import matplotlib.pyplot as plt
# import numpy as np
# dl_file = 'A:/Cian_sinmod/meeso_sim/sim_generic/DL_save.nc'
# dl = nc.Dataset(dl_file)
# l_vals = dl.variables['DL']
# d = 6
# l = 0
# light_sample = l_vals[d,l,:,:]
# t = dl.variables['time'][d]
# title_x = ('time= ' + str(t[2]) + '.' + str(t[1]) + '.' + str(t[0]) + ' ' + str(t[3]) + ':' + str(t[4]) + ':' + str(t[5]) +
#            '     ' + str('layer') + '= ' + str(l))
# plt.contourf(light_sample)
# plt.title(title_x)
# plt.colorbar()
# plt.savefig('C:/Users/ciank/PycharmProjects/sinmod/Krillmod/dl.png',dpi=400)
# #plt.show()
#

