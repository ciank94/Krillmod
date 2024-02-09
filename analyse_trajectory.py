import os
import numpy as np
import netCDF4 as nc
import xarray as xr


def sim_account(store):
    # Note: Adapt this file for  new data storage
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


def particle_visits(reg_file):
    from Krillmod.get_trajectory import region_part
    from Krillmod.import_settings import depth, sv_dir

    # Trajectory data reformatted:
    nc_file = nc.Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
    area_idx = region_part(reg_file)
    key_list = list(area_idx.keys())

    for area_name in key_list:
        idv = (area_idx[area_name])
        df = np.zeros(np.shape(depth))
        df[np.isnan(depth)] = np.nan
        df2 = np.zeros(np.shape(depth))
        df2[np.isnan(depth)] = np.nan
        x = nc_file['xp'][idv, :].astype(int)
        y = nc_file['yp'][idv, :].astype(int)
        for i in range(0, np.shape(x)[0]):
            yi = y[i, :]
            xi = x[i, :]
            df[yi, xi] = df[yi, xi] + 1
            for j in range(0, np.shape(xi)[0]):
                df2[yi[j], xi[j]] = df2[yi[j], xi[j]] + 1

        file1 = sv_dir + area_name + '_dom_paths.npy'
        file2 = sv_dir + area_name + '_total_visits.npy'
        file3 = sv_dir + area_name + '_n_parts.txt'

        # Save matrices to intermediate file;
        np.save(file1, df)
        np.save(file2, df2)

        # Save unique particle number in region:
        line1 = 'Number of particles in ' + area_name + ' region: \n'
        line2 = str(np.shape(x)[0])
        f = open(file3, 'w')
        f.writelines(line1)
        f.writelines(line2)
        f.close()

    nc_file.close()  # close nc file
    return



