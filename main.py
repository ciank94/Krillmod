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

    p = Plots()
    p.plot_depth(f, k_t)
    p.plot_back(f)
    breakpoint()

    shp_t = np.shape(k_t.x)[0]
    shp_i = np.shape(k_t.x)[1]
    SG_x = 520
    SG_y = 750
    site_matrix = np.zeros(k_t.n_sites)
    time_id = np.zeros(shp_i)
    in_region = np.zeros(shp_i)

    for i in range(0, shp_t, 10):





        print(i)
        [x1, y1] = k_t.xy_slice(i)
        e_dist = np.sqrt(((x1 - SG_x) ** 2) + (y1 - SG_y) ** 2)
        time_id[(e_dist<50) & (in_region == 0)] = i
        in_region[e_dist<50] = np.nan
    breakpoint()



    #todo: call functions that contain vectors for processing different conditions e.g. visit_region-;
    #todo: make regions


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
