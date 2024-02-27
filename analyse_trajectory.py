import numpy as np
import netCDF4 as nc
import sys
import os
from import_settings import make_directory


def lagrangian_analysis(comp_node, list_dir, sub_idx):
    # Arguments:
    # list_dir: target directories for saving intermediate output files
    # sub_idx: Dictionary with n keys each containing p x 1 boolean vectors with True for subset of individuals
    # NOTE: Modify function to save keys in separate folders;
    # Also use subsets for all analysis at once to increase efficiency
    key_list = list(sub_idx.keys())
    print('Key list for analysis: ')
    print(key_list)
    for sub_key in key_list:
        idv = (sub_idx[sub_key])  # Boolean vector for subset of individuals
        list_dir = sub_folder(comp_node, list_dir, sub_key)  # Creates sub_folder for subset of individuals

        # Dominant pathways algorithm for subset of individuals (Van Sebille paper)
        list_dir['dom_file'] = list_dir[sub_key + '_folder'] + 'dominant_paths.npy'
        if not os.path.exists(list_dir['dom_file']):
            dominant_paths(idv, list_dir)
        else:
            print('Directory: ' + list_dir['dom_file'] + ' already exists, skipping')

        # Transit time distribution (Van Sebille mainly): os.path.exists() exception handled inside function
        transit_times(idv, list_dir, sub_key)

        # Finite-size Lyapunov exponent (FSLE)-
        # (Bettencourt mainly: https://www.nature.com/articles/ngeo2570 & check email)
        # Connectivity estimates;
    return list_dir


def dominant_paths(idv, list_dir):
    # Load depth file and region file:
    depth = np.load(list_dir['depth_file'])
    nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
    df = np.zeros(np.shape(depth))
    df[np.isnan(depth)] = np.nan
    x = nc_file['xp'][idv, :].astype(int)
    y = nc_file['yp'][idv, :].astype(int)
    for i in range(0, np.shape(x)[0]):
        yi = y[i, :]
        xi = x[i, :]
        df[yi, xi] = df[yi, xi] + 1
    df[df > 0] = ((df[df > 0]) / np.shape(x)[0]) * 100
    print('Saving: ' + list_dir['dom_file'])
    np.save(list_dir['dom_file'], df)
    nc_file.close()  # close nc file
    return


def transit_times(idv, list_dir, sub_key):
    subs = ssmu_target(sub_key)  # Target regions stored in dictionary (there may be multiple target areas);
    for target_area in subs:
        filepath = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit.npy'
        filepath2 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormx.npy'
        filepath3 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormy.npy'
        filepath4 = list_dir[sub_key + '_folder'] + target_area + '_' + 'act_part.npy'
        if not os.path.exists(filepath):
            # Load depth file and region file:
            depth = np.load(list_dir['depth_file'])
            nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
            df = np.zeros(np.shape(depth))
            df[np.isnan(depth)] = np.nan
            x = nc_file['xp'][idv, :]
            y = nc_file['yp'][idv, :]
            in_reg = nc_file['in_region'][idv, :]
            act_part = nc_file['act_part'][idv]
            time = np.load(list_dir['time_file'], allow_pickle=True)
            transit_mat = np.zeros([np.shape(x)[0], 5])
            n_saves = 2500
            worm_matx = np.empty([n_saves, np.shape(x)[1]]) # Too much information for plotting, put a limit on the first thousand, also create a function for saving the first batches of particles for plotting
            worm_maty = np.empty([n_saves, np.shape(x)[1]])
            sub_ids = subs[target_area]
            print(subs[target_area])
            c = -1
            for i in range(0, np.shape(x)[0]):
                act_id = act_part[i]  # when particle becomes active
                visit_reg = in_reg[i, :]  # All the regions particle visits
                if not np.isin(in_reg[i, act_id], sub_ids): #Make sure it doesn't start in region
                    ids = np.where(np.isin(visit_reg, sub_ids)) # where there are overlaps
                    if not np.shape(ids)[0]*np.shape(ids)[1] == 0:
                        transit_mat[i, 0] = x[i, act_id]
                        transit_mat[i, 1] = y[i, act_id]
                        date_0 = time[act_id]  # time particle became active
                        date_1 = time[ids[0][0]]  # time particle reached the target destination;
                        timeframe = date_1 - date_0
                        transit_hours = (timeframe.days*24) + np.floor(timeframe.seconds*1/(60*60))
                        transit_mat[i, 2] = transit_hours
                        transit_mat[i, 3] = x[i, ids[0][0]]
                        transit_mat[i, 4] = y[i, ids[0][0]]
                        #Worm matrix
                        if not c == n_saves:
                            worm_matx[c, act_id:np.shape(x)[1]] = x[i, act_id:np.shape(x)[1]]
                            worm_maty[c, act_id:np.shape(x)[1]] = y[i, act_id:np.shape(x)[1]]
                            c = c + 1

                    # Add column with index of arrival perhaps
            print('Saving: ' + filepath)
            np.save(filepath, transit_mat)
            np.save(filepath2, worm_matx)
            np.save(filepath3, worm_maty)
            # when particles are active:
            nc_file.close()  # close nc file
        else:
            print('Directory: ' + filepath + ' already exists, skipping')
    return


def sim_account(list_dir):
    # Note: Adapt this file for  new data storage
    # Access regional file and look at number of particles, number of batches, individuals in each batch
    summary_file = list_dir['save_folder'] + list_dir['sim_folder'] + '/sim_summary.txt'
    if not os.path.exists(summary_file):
        print('')
        print('Saving: ' + summary_file)
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
        print('')
        print('Directory: ' + summary_file + ' already exists, skipping')
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
    area_list = ['WAP', 'SO', 'ALL']
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


def ssmu_target(sub_key):
    subs = dict()
    if sub_key == 'WAP':
        subs['SO'] = np.arange(9, 11 + 1)  # SO neritic zone target
        subs['SG'] = np.arange(12, 14 + 1)  # SG neritic zone target
    elif sub_key == 'SO':
        subs['SG'] = np.arange(13, 14 + 1)
    elif sub_key == 'ALL':
        subs['SO'] = np.arange(9, 11 + 1)
        subs['SG'] = np.arange(13, 14 + 1)
    else:
        print('Error: SSMU end area not specified for transit')
        sys.exit()
    return subs


def sub_folder(comp_node, list_dir, sub_key):
    folder1 = list_dir['save_folder'] + list_dir['sim_folder']
    folder2 = folder1 + '/' + sub_key + '/'
    list_dir[sub_key + '_folder'] = folder2
    if not os.path.exists(folder2):
        make_directory(comp_node, folder1, folder2)
    return list_dir


