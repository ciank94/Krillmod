import numpy as np
import netCDF4 as nc


def particle_visits(list_dir):
    depth = np.load(list_dir['depth_file'])
    # Trajectory data reformatted:
    nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
    area_idx = region_part(list_dir['reg_file'])

    #Modify function to accept a subset of particles generally- rather than just area indices;
    key_list = list(area_idx.keys())
    print(key_list)

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
        #    for j in range(0, np.shape(xi)[0]):
        #        df2[yi[j], xi[j]] = df2[yi[j], xi[j]] + 1

        df[df > 0] = ((df[df > 0])/np.shape(x)[0])*100

        file1 = list_dir['save_folder'] + list_dir['sim_folder'] + '/' + area_name + '_dom_paths.npy'
        #file2 = list_dir['save_folder'] + list_dir['sim_folder'] + '/' + area_name + '_total_visits.npy'
        #file3 = list_dir['save_folder'] + list_dir['sim_folder'] + '/' + area_name + '_n_parts.txt'

        print('Saving ' + area_name + ' visits')

        # Save matrices to intermediate file;

        # For any particle visiting a region; did it previously visit others?
        # Did it pass through different areas? Test on saga; Copy files to saga
        #load python packages; use module load netcdf4 python; scipy; numpy; module purge; then module load;
        np.save(file1, df)
        # np.save(file2, df2)
        #
        # # Save unique particle number in region:
        # line1 = 'Number of particles in ' + area_name + ' region: \n'
        # line2 = str(np.shape(x)[0])
        # f = open(file3, 'w')
        # f.writelines(line1)
        # f.writelines(line2)
        # f.close()

    nc_file.close()  # close nc file
    return


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


def retention_part(reg_file, time_file, sv_dir, tr_folder):
    # Trajectory data reformatted:
    nc_file = nc.Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
    area_idx = region_part(reg_file)
    key_list = list(area_idx.keys())
    time = np.load(time_file, allow_pickle=True)
    for area_name in key_list:
        idv = (area_idx[area_name])
        x = nc_file['xp'][idv, :].astype(int)
        reg_v = nc_file['in_region'][idv, :]
        start_p = nc_file['act_part'][idv]
        uniq_b = np.unique(start_p)
        shp_b = np.shape(uniq_b)[0]
        shp_t = np.shape(x)[1]
        store_t = np.zeros([np.shape(start_p)[0], 2])
        c = -1
        for j in range(0, shp_b):
            print(area_name + ': Percent complete = ' + str(np.ceil((j / shp_b) * 100)))
            id_b = start_p == uniq_b[j]
            reg_sim = reg_v[id_b, uniq_b[j]:shp_t]
            b_it = np.shape(reg_sim)[0]
            if j == 0:
                reg_uniq = np.unique(reg_sim[:, 0])

            for i in range(0, b_it):
                c = c + 1
                id_in = np.isin(reg_sim[i, :], reg_uniq)
                id_f = np.where(id_in == False)
                if len(id_f[0]) == 0:
                    store_t[c, 0] = -1
                else:
                    dt = time[id_f[0][0]] - time[0]
                    store_t[c, 0] = dt.days

        store_t[:, 1] = start_p
        file1 = sv_dir + tr_folder + '/' + area_name + '_retention_days.npy'
        np.save(file1, store_t)
    nc_file.close()
    return


def region_part(reg_file):
    nc_file = nc.Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
    start_id = nc_file['start'][:]
    #subs_so = (start_id == 9) | (start_id == 10) | (start_id == 11)
    subs_wap = (start_id == 2) | (start_id == 3) | (start_id == 4) | (start_id == 5) | (start_id == 6) | (start_id == 7)
    #subs_eap = (start_id == 17)
    #subs_apa = (start_id == 1)
    #subs_sopa = (start_id == 8)
    #area_idx = dict([('SO', subs_so), ('WAP', subs_wap), ('EAP', subs_eap), ('APP', subs_apa), ('SOP', subs_sopa)])
    area_idx = dict([('WAP', subs_wap)])
    nc_file.close()
    return area_idx
#     #% !(AP: 1:APPA; 2: APW; 3: DPW; 4: DPE; 5: BSW; 6:BSE; 7: EI; 17: APE)
#     #% !(SOI: 8: SOPA; 9: SOW; 10:SONE; 11: SOSE)
#     #% !(SG: 12: SGPA; 13: SGW; 14:SGE)
#     #% !(SSI: 15: SSPA; 16: SSI)
#     return




