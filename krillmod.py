# master file for krill model analysis
from import_settings import locate_folders
from analyse_trajectory import lagrangian_analysis, ssmu_start, sim_account
from plot_trajectory import plot_transit, plot_dom_paths, animate_transit, animate_dom_paths

# # Define directories
shape_name = 'ssmusPolygon.shp'
#shape_name = 'mpasPolygon.shp'
sim_folder = '2016'
comp_node = 'local'
list_dir = locate_folders(comp_node, sim_folder, shape_name)

# Simulation summary:
sim_account(list_dir)

# Analyse trajectories:
sub_idx = ssmu_start(list_dir['reg_file'])
list_dir = lagrangian_analysis(comp_node, list_dir, sub_idx)

# plot trajectories
#plot_transit(list_dir, sub_idx)
#plot_dom_paths(list_dir, sub_idx)

animate_transit(list_dir, sub_idx)
#animate_dom_paths(list_dir, sub_idx)

