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
        self.retention_file = f.save + f.sim + '/retention.npy'  # save dominant paths
        return

    def dom_paths(self, k):
        if not os.path.exists(self.dom_file):
            df = np.zeros(np.shape(k.depth))
            df[np.isnan(k.depth)] = np.nan
            x = k.x[:]
            y = k.y[:]
            for p_id in range(0, k.p_max):
                yi = y[p_id, :]
                xi = x[p_id, :]
                df[yi, xi] = df[yi, xi] + 1
            df[df > 0] = ((df[df > 0]) / k.p_max) * 100

            print('Saving: ' + self.dom_file)
            np.save(self.dom_file, df)
        else:
            print('Directory: ' + self.dom_file + ' already exists, skipping')
            return

    def depth_profile(self, k):
        if not os.path.exists(self.depth_file):
            max_depth = 200
            z_bins = np.zeros([max_depth, k.t_max])
            #  Find and loop over batches of individuals
            for t_id in range(0, np.shape(k.x)[0]):
                z = k.get_slice_z(t_id)
                z_ik = np.floor(z[:])
                for i in range(1, max_depth):
                    z_bins[i, t_id] = np.sum(i == z_ik)

            np.save(self.depth_file, z_bins)
        else:
            print('Directory: ' + self.depth_file + ' already exists, skipping')
        return

    def transit_times(self, key, k):
        self.ssmu_target(key)
        f_name = self.transit_file
        for keys in self.end_key:
            key_vals = self.end_key[keys]
            self.transit_file = f_name + keys + '_transit.npy'
            if not os.path.exists(self.transit_file):
                transit_mat = np.zeros([k.p_max, 5])
                for i in range(0, k.p_max):
                    act_id = k.act_part[i]  # when particle becomes active
                    visit_reg = k.in_reg[i, :]  # All the regions particle visits
                    if not np.isin(k.in_reg[i, act_id], key_vals):  # Make sure it doesn't start in region
                        ids = np.where(np.isin(visit_reg, key_vals))  # where there are overlaps
                        if not np.shape(ids)[0] * np.shape(ids)[1] == 0:
                            transit_mat[i, 0] = k.x[i, act_id]
                            transit_mat[i, 1] = k.y[i, act_id]
                            date_0 = k.time[act_id]  # time particle became active
                            date_1 = k.time[ids[0][0]]  # time particle reached the target destination;
                            timeframe = date_1 - date_0
                            transit_hours = (timeframe.days * 24) + np.floor(timeframe.seconds * 1 / (60 * 60))
                            transit_mat[i, 2] = transit_hours
                            transit_mat[i, 3] = k.x[i, ids[0][0]]
                            transit_mat[i, 4] = k.y[i, ids[0][0]]

                print('Saving ' + self.transit_file)
                np.save(self.transit_file, transit_mat)
            else:
                print('Directory: ' + self.transit_file + ' already exists, skipping')
        return

    def retention_times(self, k):
        if not os.path.exists(self.retention_file):
            df = np.zeros(np.shape(k.depth))
            df[np.isnan(k.depth)] = np.nan
            act_part = k.act_part[:]
            time = k.time
            d_lim = 5
            temp_vec = np.zeros(np.shape(k.x)[0])
            x = np.array(k.x[:])
            y = np.array(k.y[:])
            store_mat = np.array([x[0:np.sum(act_part==0),0], y[0:np.sum(act_part==0),0], np.zeros(np.shape(y[0:np.sum(act_part==0), 0]))]).T
            for i in range(0, np.shape(k.x)[0]):
                act_id = act_part[i]
                xt = x[i, act_id:np.shape(x)[1]].astype(float)
                yt = y[i, act_id:np.shape(y)[1]].astype(float)
                y_dist = (yt - yt[0])**2
                x_dist = (xt - xt[0])**2
                dist_t = np.sqrt(x_dist + y_dist)
                cond = np.where(dist_t > d_lim)
                if np.shape(cond)[0]*np.shape(cond)[1] > 0:
                    lim_v = cond[0][0] + act_id
                    date_0 = time[act_id]  # time particle became active
                    date_1 = time[lim_v]  # time particle reached the target destination;
                    timeframe = date_1 - date_0
                    transit_hours = (timeframe.days * 24) + np.floor(timeframe.seconds * 1 / (60 * 60))
                    temp_vec[i] = transit_hours
                else:
                    temp_vec[i] = np.nan

            for i in range(0, np.shape(store_mat)[0]):
                ids = np.arange(i, np.shape(x)[0], np.shape(store_mat)[0])
                store_mat[i, 2] = np.nanmean(temp_vec[ids])
            np.save(self.retention_file, store_mat)
        else:
            print('Directory: ' + self.retention_file + ' already exists, skipping')
        return

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


class AnalyseComp:
    def __init__(self, f):

        # Intermediate files:
        self.transit_file1 = f.save + f.sim + '/'
        self.dom_file1 = f.save + f.sim + '/dominant_paths.npy'  # save dominant paths
        self.depth_file1 = f.save + f.sim + '/depth_profile.npy'  # save dominant paths
        self.retention_file1 = f.save + f.sim + '/retention.npy'  # save dominant paths

        self.transit_file2 = f.save + f.sim2 + '/'
        self.dom_file2 = f.save + f.sim2 + '/dominant_paths.npy'  # save dominant paths
        self.depth_file2 = f.save + f.sim2 + '/depth_profile.npy'  # save dominant paths
        self.retention_file2 = f.save + f.sim2 + '/retention.npy'  # save dominant paths

        self.compare_path_file = f.comp_folder + 'comp_paths.npy'
        self.compare_retention_file = f.comp_folder + 'comp_retention.npy'
        return

    def compare_pathways(self):
        if not os.path.exists(self.compare_path_file):
            df1 = np.load(self.dom_file1)
            df2 = np.load(self.dom_file2)
            err = (df1-df2)
            np.save(self.compare_path_file, err)
        else:
            print('Directory: ' + self.compare_path_file + ' already exists, skipping')
        return

    def compare_retention(self):
        if not os.path.exists(self.compare_retention_file):
            df1 = np.load(self.retention_file1)
            df2 = np.load(self.retention_file2)
            idx1 = np.shape(df1)[0]
            idx2 = np.shape(df2)[0]
            err = (df1[0:(np.min([idx1, idx2])), 2] - df2[0:(np.min([idx1, idx2])), 2])
            ret_err = [df1[0:(np.min([idx1, idx2])), 0], df1[0:(np.min([idx1, idx2])), 1], err]
            np.save(self.compare_retention_file, ret_err)
        else:
            print('Directory: ' + self.compare_retention_file + ' already exists, skipping')
        return










