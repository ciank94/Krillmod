import numpy as np
import netCDF4 as nc
import sys
import os


def particle_visits(list_dir, sub_idx):
    # Load depth file and region file:
    depth = np.load(list_dir['depth_file'])
    nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')

    # NOTE: Modify function to save keys in separate folders;
    # Also use subsets for all analysis at once to increase efficiency
    key_list = list(sub_idx.keys())
    print('Key list for analysis: ')
    print(key_list)
    for sub_key in key_list:
        dom_file = list_dir['save_folder'] + list_dir['sim_folder'] + '/' + sub_key + '_dom_paths.npy'
        if not os.path.exists(dom_file):
            idv = (sub_idx[sub_key])
            df = np.zeros(np.shape(depth))
            df[np.isnan(depth)] = np.nan
            x = nc_file['xp'][idv, :].astype(int)
            y = nc_file['yp'][idv, :].astype(int)
            for i in range(0, np.shape(x)[0]):
                yi = y[i, :]
                xi = x[i, :]
                df[yi, xi] = df[yi, xi] + 1
            df[df > 0] = ((df[df > 0])/np.shape(x)[0])*100
            print('Saving ' + sub_key + ' visits')
            np.save(dom_file, df)
        else:
            print('Directory: ' + dom_file + ' already exists, skipping')
    nc_file.close()  # close nc file
    return


def sim_account(list_dir):
    # Note: Adapt this file for  new data storage
    # Access regional file and look at number of particles, number of batches, individuals in each batch
    summary_file = list_dir['save_folder'] + list_dir['sim_folder'] + '/sim_summary.txt'
    if not os.path.exists(summary_file):
        nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
        n_parts = np.shape(nc_file['xp'])[0]
        uniq_id = np.unique(nc_file['act_part'][:])
        n_releases = np.shape(uniq_id)[0]
        n_per_batch = np.floor(n_parts/n_releases).astype(int)
        time = np.load(list_dir['time_file'], allow_pickle=True)
        shp = np.shape(time)
        date_0 = time[0]
        date_1 = time[shp[0]-1]
        timeframe = date_1 - date_0
        start_date = 'Start datetime (dd.mm.yyyy) = ' + date_0.strftime("%d.%m.%Y, %H:%M:%S")
        end_date = 'End datetime (dd.mm.yyyy) = ' + date_1.strftime("%d.%m.%Y, %H:%M:%S")
        time_frame = 'Time frame = ' + str(timeframe.days) + ' days and ' + str(timeframe.seconds*1/(60*60)) + ' hours '
        f = open(summary_file, 'w')
        top_line = ['Simulation summary: \n']
        time_lines = [start_date + '\n', end_date + '\n', time_frame + '\n']
        n_parts_line1 = ['Total number of particles: ' + str(n_parts) + '\n']
        n_parts_line2 = ['Number of releases: ' + str(n_releases) + '\n']
        n_parts_line3 = ['Number per release: ' + str(n_per_batch) + '\n']
        f.writelines(top_line)
        f.writelines(time_lines)
        f.writelines(n_parts_line1)
        f.writelines(n_parts_line2)
        f.writelines(n_parts_line3)
        f.close()
        nc_file.close()  # close nc file
    else:
        print('Directory: ' + summary_file + ' already exists, skipping')
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


def ssmu_start(reg_file):
    # Hard-coded indices for particles starting in ssmu's. This can be adapted for both time and region;
    # Below is the id number for each of the ssmu regions:
    # !(AP: 1:APPA; 2: APW; 3: DPW; 4: DPE; 5: BSW; 6:BSE; 7: EI; 17: APE)
    # !(SOI: 8: SOPA; 9: SOW; 10:SONE; 11: SOSE)
    # !(SG: 12: SGPA; 13: SGW; 14:SGE)
    # !(SSI: 15: SSPA; 16: SSI)

    nc_file = nc.Dataset(reg_file, mode='r', format='NETCDF4_CLASSIC')
    start_id = nc_file['start'][:]
    sub_idx = dict()
    area_list = ['ALL']
    for reg_id in area_list:
        if reg_id == 'WAP':
            subs = np.arange(2, 7 + 1)
        elif reg_id == 'SO':
            subs = np.arange(9, 11 + 1)
        elif reg_id == 'SG':
            subs = np.arange(12, 14 + 1)
        elif reg_id == 'ALL':
            subs = np.unique(start_id)
        else:
            print('Error: SSMU start area not specified')
            sys.exit()
        idx_in = np.isin(start_id, subs)  # Boolean vector of individuals in set
        sub_idx[reg_id] = idx_in  # Set key to dictionary
    nc_file.close()
    return sub_idx





