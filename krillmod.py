# master file for krill model analysis
from import_settings import locate_folders
#from analyse_trajectory import particle_visits

# Reformat trajectory data and setup target folders:
shape_name = 'ssmusPolygon.shp'
sim_folder = '2016'
comp_node = 'local'
list_dir = locate_folders(comp_node, sim_folder, shape_name)

# Analyse trajectories
# particle_visits(list_dir)









