# master file for krill model analysis
from set_directory import Folder
from get_trajectory import Trajectory
from analyse_trajectory import Analyse
from plot_trajectory import Plots
import numpy as np
import matplotlib.pyplot as plt

# #trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
sim_list = ['2017']  # Simulation identifier
trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
node = 'local'
shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant

for s in sim_list:
    # Define names of main folders
    sim_folder = s  # Name of simulation instance

    # Setup directory information for analysis
    f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
    f.set_folder(node)  # Sets folders based on remote or local node
    f.exist()  # Checks existence of relevant folders and files

    # Two classes for dealing with data:
    k_t = Trajectory(f)

    p1 = Plots()
    p1.plot_dom_paths(f, k_t)
    breakpoint()





    shp_t = np.shape(k_t.x)[0]
    shp_i = np.shape(k_t.x)[1]

    time_id = np.zeros(shp_i)
    in_region = np.zeros(shp_i)
    df = np.zeros(np.shape(k_t.depth))
    df[np.isnan(k_t.depth)] = np.nan

    for i in range(0, shp_t, 1):
        print(i)
        [x1, y1] = k_t.xy_slice(i)

        near_SG = k_t.SG_dist(x1, y1)
        time_id[(near_SG) & (in_region == 0)] = i
        in_region[near_SG] = np.nan

        x2 = x1[x1 > 0]
        y2 = y1[x1 > 0]
        df[y2, x2] = df[y2, x2] + 1



    save_name = 'dom_paths.npy'

    df = df / np.nanmax(df)
    #todo: make save intermediate file function;
    save_file = f.save + f.sim + '/' + save_name
    print('Saving: ' + save_file)
    np.save(save_file, df)



    #np.save()




    breakpoint()
    id = np.where(time_id > 0)
    iv = id[0]
    n_ent = np.shape(iv)[0]
    ti = time_id[time_id > 0].astype(int)
    x_samp = np.zeros(n_ent)
    y_samp = np.zeros(n_ent)
    for i in range(0, 100):
        print(i)
        x_samp[i] = k_t.x[ti[i], iv[i]]
        y_samp[i] = k_t.y[ti[i], iv[i]]

    breakpoint()



    #todo: call functions that contain vectors for processing different conditions e.g. visit_region-;
    #todo: make regions
    #todo: fix the trajectory code;

    # p = Plots()
    # p.plot_depth(f, k_t)
    # p.plot_back(f)

    #
    #
    #     # Now the analysis
    # df = Analyse(f)
    # df.dom_paths(k_r)
    # #df.depth_profile(k_t)  # Note
    # df.retention_times(k_r)
    # df.transit_times('ALL', k_r)
    #
    # p = Plots()
    # p.plot_connectivity(f, k_r)
    # p.plot_transit(f)
    # p.plot_retention(f)
    # p.plot_dom_paths(f, k_r)
    #
    # #k_t.nc_file.close()
    # k_r.nc_file.close()

#
# def split_analysis(sim_list, trj_folder, shp_name, node):
#     from set_directory import Folder
#     from get_trajectory import Trajectory, Regional
#     from analyse_trajectory import Analyse
#     from plot_trajectory import Plots
#     for s in sim_list:
#         # Define names of main folders
#         sim_folder = s  # Name of simulation instance
#
#         # Setup directory information for analysis
#         f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
#         f.set_folder(node)  # Sets folders based on remote or local node
#         f.exist()  # Checks existence of relevant folders and files
#
#         # Two classes for dealing with data:
#         k_r = Regional(f)
#         breakpoint()
#         # k_t = Trajectory(f)
#
#         # Now the analysis
#         # df = Analyse(f)
#         # df.dom_paths(k_r)
#         # # df.depth_profile(k_t)  # Note
#         # df.retention_times(k_r)
#         # df.transit_times('ALL', k_r)
#         #
#         # p = Plots()
#         # p.plot_connectivity(f, k_r)
#         # p.plot_transit(f)
#         # p.plot_retention(f)
#         # p.plot_dom_paths(f, k_r)
#
#         # k_t.nc_file.close()
#         k_r.nc_file.close()
#
#
#
# def get_quant(sim_list, trj_folder, shp_name, node):
#     from set_directory import Folder
#     from plot_trajectory import Plots
#     for s in sim_list:
#         # Define names of main folders
#         sim_folder = s  # Name of simulation instance
#
#         # Setup directory information for analysis
#         f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
#         f.set_folder(node)  # Sets folders based on remote or local node
#
#         df = Plots()
#         df.save_data(f)
#
#
#
#
# def compare_2_years(t_folder1, sim_folder1, t_folder2, sim_folder2, node):
#     from set_directory import FolderComp
#     from analyse_trajectory import AnalyseComp
#     from plot_trajectory import Plots
#
#     files = FolderComp(t_folder1, sim_folder1, t_folder2, sim_folder2)
#     files.set_folder(node)
#     files.exist()
#
#     df = AnalyseComp(files)
#     df.compare_pathways()
#     df.compare_retention()
#
#     p = Plots()
#     p.plot_comp_paths(files)
#     p.plot_comp_retention(files)
#     return
#
#
#
# from sim_scenarios import main_analysis, compare_2_years, split_analysis, get_quant  # simulation scenarios
#
# # sim_list = ['2017']  # Simulation identifier
# # trj_folder = 'D:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
# # sim_list = ['DVM']
# # trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
# # shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
# # get_quant(sim_list, trj_folder, shp_name, node='local')
# #
# # breakpoint()
# # 1: Main analysis for a list of simulations;
# sim_list = ['2017']  # Simulation identifier
# trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
# #trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
# shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
# main_analysis(sim_list, trj_folder, shp_name, node='local')
#
#
#
# split_analysis(sim_list, trj_folder, shp_name, node='local')
# breakpoint()
#
#
#
# sim_list = ['2017']  # Simulation identifier
# trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
# #trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
# shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
# main_analysis(sim_list, trj_folder, shp_name, node='local')
#
# breakpoint()
#
# # 2: Comparisons between 2 simulations
# sim_folder1 = 'DVM'  # Simulation identifier
# sim_folder2 = '2017'  # Second simulation identifier
# t_folder1 = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
# t_folder2 = 'A:/Cian_sinmod/sim_'  # trajectory folder- saga
# compare_2_years(t_folder1, sim_folder1, t_folder2, sim_folder2, node='local')
#
#
