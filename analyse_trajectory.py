import os
import numpy as np


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
        start_vec = np.zeros(np.shape(store_reg))
        store_dict = {'start_time': np.zeros([np.shape(store_reg)[0], 2]),
                  'exit_time': np.zeros([np.shape(store_reg)[0], 2]),
                  'start_region': np.zeros([np.shape(store_reg)[0], 1]),
                  'exit_region': np.zeros([np.shape(store_reg)[0], 1])}
        shp_t = np.shape(store['time'])
        for idv in range(0, shp_t[0]):
            print('Percent complete = ' + str(np.ceil((idv / shp_t[0]) * 100)))
            store = store_traj(tr_file, idv)
            store_reg = store_regions(reg_file, idv)
            act_ind = np.sum(store['active'] == 1)
            if idv == 0:
                store_dict['start_time'][0:act_ind, 0] = store['time'][idv].day
                store_dict['start_time'][0:act_ind, 1] = store['time'][idv].hour
                store_dict['start_region'][0:act_ind, 0] = store_reg[0:act_ind]
                start_vec[store_reg > 0] = 1  # Assign all individuals an active number
                t_act = act_ind
            else:
                current_vec = np.zeros(np.shape(start_vec))
                current_vec[store_reg > 0] = 1
                diff_v = start_vec - current_vec
                if len(diff_v[diff_v == 1]) > 0:
                    store_dict['exit_region'][diff_v == 1, 0] = store_reg[diff_v == 1][:]
                    store_dict['exit_time'][diff_v == 1, 0] = store['time'][idv].day
                    store_dict['exit_time'][diff_v == 1, 1] = store['time'][idv].hour
            if act_ind > t_act:
                new_reg = store_reg[t_act:act_ind]
                new_reg[new_reg < 1] = 0
                store_dict['start_time'][t_act:act_ind, 0] = store['time'][idv].day
                store_dict['start_time'][t_act:act_ind, 1] = store['time'][idv].hour
                store_dict['start_region'][t_act:act_ind, 0] = store_reg[t_act:act_ind]
                start_vec[t_act:act_ind] = new_reg
                t_act = act_ind

        np.save(file_inter, store_dict)
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


