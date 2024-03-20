from get_trajectory import geo2grid
import os
import sys
import numpy as np
import geopandas as gpd
import netCDF4 as nc
from netCDF4 import num2date


class Folder:

    def __init__(self, tr_file):
        self.comp_node = None
        self.depth_file = None
        self.time_file = None
        self.reg_file = None
        self.save = None
        self.trj_file = tr_file

    def set_folder(self, comp_node):
        # Sets folder and file names based on whether the files are stored locally or on a remote server, where
        # comp node can be a number of switches
        nodes = ['local', 'saga']
        if comp_node not in nodes:
            print('Expected one of the following nodes for comp_node:  ')
            [print(i) for i in nodes]
            sys.exit()

        if comp_node == 'local':
            self.comp_node = 'local'
            self.save = 'C:/Users/ciank/PycharmProjects/sinmod/sewagemod/'
        elif comp_node == 'saga':
            self.comp_node = 'saga'
            self.save = '/cluster/projects/nn9828k/Cian_sinmod/python/Krillmod/results/'
        else:
            print('Error: computation node not specified')
            sys.exit()

        # Define files based on folder setup
        # 1) Input files for analysis
        self.time_file = self.save + '/times.npy'
        self.depth_file = self.save + '/depth.npy'

        return

    def exist(self):

        if not os.path.exists(self.trj_file):
            print('Error: Directory "' + self.trj_file + '" does not exist\nCheck connection to directory')
            sys.exit()
        else:
            print(self.trj_file + ' exists')

        if not os.path.exists(self.save):
            self.make_directory()
        else:
            print(self.save + ' exists')

        if not os.path.exists(self.time_file):
            self.get_times()  # Creates time file in folder
        else:
            print(self.time_file + ' exists')

        if not os.path.exists(self.depth_file):
            self.get_depth()  # Creates depth file in folder
        else:
            print(self.depth_file + ' exists')

        print('### NOTE:')
        print('All intermediate files exist for analysis')
        print('###')
        print('')
        return

    def make_directory(self):
        folder1 = self.save
        # Work around for remote server when we need to create two folders at once;
        if self.comp_node == 'local':
            os.mkdir(folder1)
        else:
            sys.exit('Need to use create remote directory for saving analysis output')
            # N
            # cmd = 'mkdir ' + folder1
            # cmd2 = 'mkdir ' + folder2
            # os.system(cmd)
            # os.system(cmd2)
        return

    def get_times(self):
        nc_file = nc.Dataset(self.trj_file)
        times = nc_file.variables['time']
        date_save = num2date(times, times.units)
        np.save(self.time_file, date_save.data)
        nc_file.close()
        return

    def get_depth(self):
        nc_file = nc.Dataset(self.trj_file)
        depth = nc_file.variables['depth'][:]
        dp = depth[:]  # get depth matrix
        idx = dp < 0  # Fill values
        dp[idx] = np.nan  # Set fill values to invalid;
        dp2 = np.array(dp)
        np.save(self.depth_file, dp2)
        nc_file.close()
        return



