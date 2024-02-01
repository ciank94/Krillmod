# master file for krill model analysis
from Krillmod.import_settings import *

# Trajectory data in dictionary format
# store = store_traj(tr_file)    # Read trajectory file
# shape_v = read_ssmu(shp_file)  # Read ssmu shapefiles

poly_id = read_regions()
new_mat = np.zeros(np.shape(poly_id))


new_mat[poly_id == 1] = 1
new_mat[poly_id == 2] = 1


ps = poly_id[0, :]
poly_np = [1, 2]
shp_p = np.shape(poly_np)
shp_t = np.shape(poly_id)



for t in range(0, shp_t[0]):
    for i in range(0, shp_p[0]):
        region_id[region_id == 1] = 1


#get_regions(store, shape_v)



# np.save('myfile.npy', list_v)
# np.load('myfile.npy',allow_pickle='TRUE')



#df = dom_path(store)
#plot_dom_path(df, store)
#sim_account(store)
# Save the array to a binary file
#file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/dom_paths.npy'
#np.save(file, np.array(df))

# Load the saved array
#df = np.load(file)
#
#plot_grid(file)


# First index should be for





# field = dom_path(store)


