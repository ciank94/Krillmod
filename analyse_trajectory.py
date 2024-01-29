import numpy as np


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


