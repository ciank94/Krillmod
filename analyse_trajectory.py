import os
import numpy as np
import xarray as xr


# Retention in regions
def retention(tr_file, reg_file):
    from Krillmod.get_trajectory import store_regions, store_traj
    from Krillmod.import_settings import sv_dir
    file_inter = sv_dir + str('retention.npy')
    if not os.path.exists(file_inter):
        print('Note: Creating intermediate .npy file with retention dictionary')
        idv = 0
        store = store_traj(tr_file, idv)
        store_reg = store_regions(reg_file, idv)

        nc_file = nc.Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
        # nc_file.variables.keys()

        act_ind = np.array(nc_file.variables["act_part"][:])
        n_time = nc_file.variables["in_region"][0, :]
        list_start = np.unique(act_ind)

        id1 = list_start[0]
        log_id1 = act_ind == id1
        in_polt = np.array(nc_file.variables["in_region"][log_id1, id1:np.shape(n_time)[0]])
        idp1 = in_polt[:, 0] == 1






        # store_pol = dict([("in_reg", in_polt)])
        nc_file.close()



        Init_reg = np.zeros(np.shape(store_reg))
        Exit_reg = Init_reg
        count_id = Exit_reg
        Init_time = np.zeros([np.shape(store_reg)[0], 2])
        Exit_time = Init_time
        shp_t = store['t_id']
        store_mat = np.zeros([np.shape(store_reg)[0], shp_t[0]])


        for idv in range(0, shp_t[0]):

            print('Percent complete = ' + str(np.ceil((idv / shp_t[0]) * 100)))
            # store_reg = store_regions(reg_file, idv)
            store
            store_mat[:, idv] = store_reg




            print('Percent complete = ' + str(np.ceil((idv / shp_t[0]) * 100)))
            store = store_traj(tr_file, idv)
            store_reg = store_regions(reg_file, idv)
            if idv == 0:
                idx = store_reg > 0
                Init_reg[idx] = store_reg[idx]  # Initial regions
                Init_time[idx, 0] = store['time'].day
                Init_time[idx, 1] = store['time'].hour
                act_id = store['active']
            else:
                # Check for new individuals using the active flag
                act_id2 = store['active']
                new_id = act_id2 - act_id
                act_id = act_id2
                idx2 = new_id[new_id == 1]
                Init_reg[idx2] = store_reg[idx2]  # Initial regions
                Init_time[idx2, 0] = store['time'].day
                Init_time[idx2, 1] = store['time'].hour

                #Check for exit of region:
                idx3 = ((store_reg - Init_reg) != 0) & (Init_reg > 0) & (count_id != 1)
                Exit_reg[idx3] = store_reg[idx3]
                Exit_time[idx3, 0] = store_reg[idx3]
                Exit_time[idx3, 1] = store_reg[idx3]
                count_id[idx3] = 1  # Mark individual as accounted for
    else:
        store_dict = np.load(file_inter, allow_pickle=True).item()
    return store_dict



def sim_account(store):

    shp = np.shape(store['xp'])
    # x = np.array(store['xp'][0:shp[1], :])
    # y = np.array(store['yp'][0:shp[1], :])

    # Table with information:
    # Particle information:
    n_parts = shp[1]
    date_0 = store['time'][0]
    date_1 = store['time'][shp[0]-1]
    timeframe = date_1 - date_0
    start_date = 'Start datetime (dd.mm.yyyy) = ' + date_0.strftime("%d.%m.%Y, %H:%M:%S")
    end_date = 'End datetime (dd.mm.yyyy) = ' + date_1.strftime("%d.%m.%Y, %H:%M:%S")
    time_frame = 'Time frame = ' + str(timeframe.days) + ' days and ' + str(timeframe.seconds*1/(60*60)) + ' hours '
    # Add time steps  here;
    # Write to text file
    file = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/sim_summary.txt'
    f = open(file, 'w')
    top_line = ['Simulation summary: \n']
    time_lines = [start_date + '\n', end_date + '\n', time_frame + '\n']
    f.writelines(top_line)
    f.writelines(time_lines)
    f.close()
    return


def dom_path(store):
    shp = np.shape(store['xp'])
    array_shape = np.shape(store['depth'])
    # dom_path = np.empty(array_shape)
    df = store['depth']
    df[df > 0] = 0
    x = np.array(store['xp'])
    y = np.array(store['yp'])
    av = np.array(store['active'])

    for i in range(0, shp[1]):
        act_id = av[:, i]
        xp = np.floor(x[act_id == 1, i])
        yp = np.floor(y[act_id == 1, i])
        idx = sub2ind(array_shape, xp, yp)  # Outputs positional indexes for all times for individual i
        idc = np.unique(np.floor(idx), return_index=True)  # First index where there are unique x and y
        xi = np.floor(xp[idc[1]]).astype(int)  # Subset unique visits
        yi = np.floor(yp[idc[1]]).astype(int)
        # xi = xi[xi > 0]
        # yi = yi[yi > 0]
        if np.shape(xi) != np.shape(yi):
            breakpoint()

        df[yi, xi] = df[yi, xi] + 1
    return df


def sub2ind(array_shape, rows, cols):
    ind = rows * array_shape[1] + cols
    return ind


