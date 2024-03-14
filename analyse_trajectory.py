import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import sys
import os


class Analyse:

    def __init__(self, f):

        self.t_mat = None
        self.begin_key = dict()  # specify subset of individuals based on beginning key
        self.end_key = dict()  # specify the target end point of individuals

        # Save files:
        self.transit_file = f.save + f.sim + '/'
        self.dom_file = f.save + f.sim + '/dominant_paths.npy'  # save dominant paths
        self.depth_file = f.save + f.sim + '/depth_profile.npy'  # save dominant paths
        return

    def dom_paths(self, k):
        depth = k.depth
        df = np.zeros(np.shape(depth))
        df[np.isnan(depth)] = np.nan
        for t_id in range(0, k.t_max):
            x_act, y_act, bool_act = k.get_slice(t_id)
            df[y_act.astype(int), x_act.astype(int)] = df[y_act.astype(int), x_act.astype(int)] + 1
        df[df > 0] = ((df[df > 0]) / k.n_parts) * 100
        print('Saving: ' + self.dom_file)
        np.save(self.dom_file, df)
        return

    def depth_profile(self, k):
        max_depth = 200
        z_bins = np.zeros([max_depth, k.t_max])
        #  Find and loop over batches of individuals
        for t_id in range(0, np.shape(k.x)[0]):
            z = k.get_slice(t_id)
            z_ik = np.floor(z[:])
            for i in range(1, max_depth):
                z_bins[i, t_id] = np.sum(i == z_ik)

        np.save(self.depth_file, z_bins)
        return

        # import matplotlib.pyplot as plt
        # f1 = plt.contourf(z_bins)
        # ax = f1.axes
        # ax.invert_yaxis()
        # breakpoint()

    def transit_times(self, key, k):
        self.ssmu_target(key)
        # depth = k.depth
        # df = np.zeros(np.shape(depth))
        # df[np.isnan(depth)] = np.nan
        f_name = self.transit_file
        for keys in self.end_key:
            key_vals = self.end_key[keys]
            self.transit_file = f_name + keys + '_transit.npy'

            self.t_mat = np.zeros([np.shape(k.x)[1]])
            td_mat = np.zeros([np.shape(k.x)[1]])
            ret_mat = np.zeros([np.shape(k.x)[1]])

            for t_id in range(0, np.shape(k.x)[0]):
                x_act, y_act, bool_act = k.get_slice(t_id)
                in_reg = k.in_region(x_act.astype(int), y_act.astype(int))

                if t_id == 0:
                    start_region = in_reg
                    start_idx = np.isin(start_region, key_vals)
                    self.t_mat[bool_act] = start_idx
                    self.t_mat[self.t_mat == 1] = np.nan
                else:
                    idx = np.isin(in_reg, key_vals)
                    temp_vec = self.t_mat[bool_act]
                    idv = (temp_vec == 0) * idx
                    self.t_mat[bool_act] = idv*t_id

            print('Saving ' + self.transit_file)
            np.save(self.transit_file, self.t_mat)




                #else

                #todo: save sites instead of all particles;




    #     subs = ssmu_target(sub_key)  # Target regions stored in dictionary (there may be multiple target areas);
    #     for target_area in subs:
    #         filepath = list_dir[sub_key + '_folder'] + target_area + '_' + 'transit.npy'
    #         filepath2 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormx.npy'
    #         filepath3 = list_dir[sub_key + '_folder'] + target_area + '_' + 'wormy.npy'
    #         filepath4 = list_dir[sub_key + '_folder'] + target_area + '_' + 'act_part.npy'
    #         if not os.path.exists(filepath):
    #             # Load depth file and region file:
    #             depth = np.load(list_dir['depth_file'])
    #             nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
    #             df = np.zeros(np.shape(depth))
    #             df[np.isnan(depth)] = np.nan
    #             x = nc_file['xp'][idv, :]
    #             y = nc_file['yp'][idv, :]
    #             in_reg = nc_file['in_region'][idv, :]
    #             act_part = nc_file['act_part'][idv]
    #             time = np.load(list_dir['time_file'], allow_pickle=True)
    #             transit_mat = np.zeros([np.shape(x)[0], 5])
    #             n_saves = 2500
    #             worm_matx = np.empty([n_saves, np.shape(x)[1]])
    #             worm_maty = np.empty([n_saves, np.shape(x)[1]])
    #             sub_ids = subs[target_area]
    #             print(subs[target_area])
    #             c = -1
    #             for i in range(0, np.shape(x)[0]):
    #                 act_id = act_part[i]  # when particle becomes active
    #                 visit_reg = in_reg[i, :]  # All the regions particle visits
    #                 if not np.isin(in_reg[i, act_id], sub_ids):  # Make sure it doesn't start in region
    #                     ids = np.where(np.isin(visit_reg, sub_ids))  # where there are overlaps
    #                     if not np.shape(ids)[0]*np.shape(ids)[1] == 0:
    #                         transit_mat[i, 0] = x[i, act_id]
    #                         transit_mat[i, 1] = y[i, act_id]
    #                         date_0 = time[act_id]  # time particle became active
    #                         date_1 = time[ids[0][0]]  # time particle reached the target destination;
    #                         timeframe = date_1 - date_0
    #                         transit_hours = (timeframe.days*24) + np.floor(timeframe.seconds*1/(60*60))
    #                         transit_mat[i, 2] = transit_hours
    #                         transit_mat[i, 3] = x[i, ids[0][0]]
    #                         transit_mat[i, 4] = y[i, ids[0][0]]
    #                         #Worm matrix
    #                         if not c == n_saves:
    #                             worm_matx[c, act_id:np.shape(x)[1]] = x[i, act_id:np.shape(x)[1]]
    #                             worm_maty[c, act_id:np.shape(x)[1]] = y[i, act_id:np.shape(x)[1]]
    #                             c = c + 1
    #
    #                     # Add column with index of arrival perhaps
    #             print('Saving: ' + filepath)
    #             np.save(filepath, transit_mat)
    #             np.save(filepath2, worm_matx)
    #             np.save(filepath3, worm_maty)
    #             # when particles are active:
    #             nc_file.close()  # close nc file
    #         else:
    #             print('Directory: ' + filepath + ' already exists, skipping')
    #     return
    #

    def ssmu_target(self, key):
        # Hard-coded indices for particles starting in ssmu's. This can be adapted for both time and region;
        # Below is the id number for each of the ssmu regions:
        # !(AP: 1:APPA; 2: APW; 3: DPW; 4: DPE; 5: BSW; 6:BSE; 7: EI; 17: APE)
        # !(SOI: 8: SOPA; 9: SOW; 10:SONE; 11: SOSE)
        # !(SG: 12: SGPA; 13: SGW; 14:SGE)
        # !(SSI: 15: SSPA; 16: SSI)
        self.begin_key = key
        if key == 'WAP':
            self.end_key['SO'] = np.arange(9, 11 + 1)  # SO neritic zone target
            self.end_key['SG'] = np.arange(12, 14 + 1)  # SG neritic zone target
        elif key == 'SO':
            self.end_key['SG'] = np.arange(13, 14 + 1)
        elif key == 'ALL':
            self.end_key['SO'] = np.arange(9, 11 + 1)
            self.end_key['SG'] = np.arange(13, 14 + 1)
        else:
            print('Error: SSMU end area not specified for transit')
            sys.exit()
        return




#
# def lagrangian_analysis(comp_node, list_dir, sub_idx):
#     # Arguments:
#     # list_dir: target directories for saving intermediate output files
#     # sub_idx: Dictionary with n keys each containing p x 1 boolean vectors with True for subset of individuals
#     # NOTE: Modify function to save keys in separate folders;
#     # Also use subsets for all analysis at once to increase efficiency
#     key_list = list(sub_idx.keys())
#     print('Key list for analysis: ')
#     print(key_list)
#     for sub_key in key_list:
#         idv = (sub_idx[sub_key])  # Boolean vector for subset of individuals
#         list_dir = sub_folder(comp_node, list_dir, sub_key)  # Creates sub_folder for subset of individuals
#
#         # Dominant pathways algorithm for subset of individuals (Van Sebille paper)
#         list_dir['dom_file'] = list_dir[sub_key + '_folder'] + 'dominant_paths.npy'
#         if not os.path.exists(list_dir['dom_file']):
#             dominant_paths(idv, list_dir)
#         else:
#             print('Directory: ' + list_dir['dom_file'] + ' already exists, skipping')
#
#         # Transit time distribution (Van Sebille mainly): os.path.exists() exception handled inside function
#         #transit_times(idv, list_dir, sub_key)
#
#         # Retention:
#         list_dir['ret_file'] = list_dir[sub_key + '_folder'] + 'retention.npy'
#         if not os.path.exists(list_dir['ret_file']):
#             retention_times(idv, list_dir, sub_key)
#         else:
#             print('Directory: ' + list_dir['ret_file'] + ' already exists, skipping')
#
#         #depth_times(idv, list_dir)
#
#
#         # Finite-size Lyapunov exponent (FSLE)-
#         # (Bettencourt mainly: https://www.nature.com/articles/ngeo2570 & check email)
#         # Connectivity estimates;
#     return list_dir
#
#

#
#
#

#     breakpoint()
#
#     #plt.gca().invert_yaxis()
#     #print('Saving: ' + list_dir['dom_file'])
#     plt.close()
#     nc_file.close()
#     return
#
# def retention_times(idv, list_dir, sub_key):
#     depth = np.load(list_dir['depth_file'])
#     nc_file = nc.Dataset(list_dir['reg_file'], mode='r', format='NETCDF4_CLASSIC')
#     df = np.zeros(np.shape(depth))
#     df[np.isnan(depth)] = np.nan
#     x = nc_file['xp'][idv, :].astype(int)
#     y = nc_file['yp'][idv, :].astype(int)
#     act_part = nc_file['act_part'][idv]
#     store_mat = np.array([x[0:np.sum(act_part==0),0], y[0:np.sum(act_part==0), 0], np.zeros(np.shape(y[0:np.sum(act_part==0), 0]))]).T
#     time = np.load(list_dir['time_file'], allow_pickle=True)
#     d_lim = 5
#     temp_vec = np.zeros(np.shape(x)[0])
#     for i in range(0, np.shape(x)[0]):
#         act_id = act_part[i]
#         xt = x[i, act_id:np.shape(x)[1]]
#         yt = y[i, act_id:np.shape(x)[1]]
#         dist_t = np.sqrt(((xt-xt[0])**2) + ((yt-yt[0])**2))
#         cond = np.where(dist_t > d_lim)
#         if np.shape(cond)[0]*np.shape(cond)[1] > 0:
#             lim_v = cond[0][0] + act_id
#             date_0 = time[act_id]  # time particle became active
#             date_1 = time[lim_v]  # time particle reached the target destination;
#             timeframe = date_1 - date_0
#             transit_hours = (timeframe.days * 24) + np.floor(timeframe.seconds * 1 / (60 * 60))
#             temp_vec[i] = transit_hours
#         else:
#             temp_vec[i] = np.nan
#
#     for i in range(0, np.shape(store_mat)[0]):
#         ids = np.arange(i,np.shape(x)[0], np.shape(store_mat)[0])
#         store_mat[i, 2] = np.nanmean(temp_vec[ids])
#
#     save_path = list_dir[sub_key + '_folder'] + 'retention.npy'
#     np.save(save_path, store_mat)
#     nc_file.close()
#     return
#
#









