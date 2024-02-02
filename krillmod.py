# master file for krill model analysis
from Krillmod.import_settings import *

# Trajectory data in dictionary format
# store = store_traj(tr_file)    # Read trajectory file
# shape_v = read_ssmu(shp_file)  # Read ssmu shapefiles


# Initialization:
ret = retention(tr_file, reg_file)

reg_vals = ret['start_region']

ids = np.where(reg_vals == 1)
ids[0]

st_t = ret['start_time'][ids[0], :]
et_t = ret['exit_time'][ids[0], :]


id_leave = np.where(et_t[:, 1] > 0)



np.shape(ids)

et_t[et_t[:,1]>0,:]


idv = 0
store_reg = store_regions(reg_file, idv)
store_tra = store_traj(tr_file, idv)
part_act = store_reg[store_tra['active'] == 1]



