# master file for krill model analysis
from import_settings import locate_folders
from analyse_trajectory import particle_visits, ssmu_start, sim_account
from plot_trajectory import plot_dom_paths

# Reformat trajectory data and setup target folders:
shape_name = 'ssmusPolygon.shp'
sim_folder = '2018'
comp_node = 'local'
list_dir = locate_folders(comp_node, sim_folder, shape_name)

# Simulation summary:
sim_account(list_dir)

# Analyse trajectories
sub_idx = ssmu_start(list_dir['reg_file'])
particle_visits(list_dir, sub_idx)

# plot trajectories
#plot_dom_paths(list_dir, sub_idx)








