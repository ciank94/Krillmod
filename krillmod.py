# master file for krill model analysis
from Krillmod.import_settings import *

# Trajectory data in dictionary format
store = store_traj(tr_file)    # Read trajectory file
shape_v = read_ssmu(shp_file)  # Read ssmu shapefiles


list_v = get_regions(store, shape_v)
np.save('myfile.npy', list_v)
np.load('myfile.npy',allow_pickle='TRUE')



df = dom_path(store)
plot_dom_path(df, store)
sim_account(store)
# Save the array to a binary file
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#np.save(file, np.array(df))

# Load the saved array
#df = np.load(file)
#
#plot_grid(file)


# First index should be for


list_v = ["R1", "R2"]
foo = {}
foo['R1'][0] = lat1

idx = 0







point1 = lon[0:10], lat[0:10]





# field = dom_path(store)


