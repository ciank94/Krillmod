# master file for krill model analysis
from Krillmod.import_settings import init_folders
from Krillmod.analyse_trajectory import particle_visits

# Reformat trajectory data and setup target folders:
list_dir = init_folders()

# Analyse trajectories
particle_visits(list_dir)



from Krillmod.get_trajectory import store_traj
from Krillmod.analyse_trajectory import retention_part, particle_visits
from Krillmod.plot_trajectory import plot_retention, plot_dom_paths



# for sim_folder in ['2016', '2017', '2018', '2019']:

    #particle_visits(reg_file, sv_dir, tr_folder, depth_file)
    #plot_dom_paths(tr_folder, reg_file, sv_dir)
    #retention_part(reg_file, time_file, sv_dir, tr_folder)
    #plot_retention(sv_dir, tr_folder, time_file, reg_file)

# To do: Adapt this code for analysing the vertical depth of particles;
#zv = np.zeros([np.shape(range(0, 700))[0]])
#for idv in range(0,700):
    #store = store_traj(tr_file, idv)
    #if idv == 1:
        #act = store['active']
        #zv[idv] = np.mean(store['zp'][act == 1])
    #else:
        #zv[idv] = np.mean(store['zp'][act == 1])
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


from mpl_toolkits.axes_grid1 import make_axes_locatable
from Krillmod.get_trajectory import *
from Krillmod.analyse_trajectory import *
from Krillmod.plot_trajectory import *

from Krillmod.import_settings import *


#fig = plt.figure()


#for i in ['2016', '2017', '2018', '2019']:
 #   dom_file = sv_dir + tr_folder + '/WAP_dom_paths.npy'
    # import_files(tr_folder)
  #  plot_dom_paths(dom_file, tr_folder)
   # tr_folder = i
    #from Krillmod.import_settings import *
        # import_files(tr_folder)




    #writer.save(sv_dir + 'writer_test.mp4', writer = FFwriter)
    #writer.grab_frame()
#particle_visits(reg_file)







