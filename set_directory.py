from get_trajectory import geo2grid
import os
import sys
import numpy as np
import geopandas as gpd
import netCDF4 as nc
from netCDF4 import num2date


class Folder:

    def __init__(self, t_folder, sim_folder, shape_name):
        self.poly_file = None
        self.comp_node = None
        self.depth_file = None
        self.time_file = None
        self.reg_file = None
        self.shp_file = None
        self.trj_file = None
        self.save = None
        self.shapes = None
        self.traj = t_folder
        self.sim = sim_folder
        self.shp_name = shape_name

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
            self.shapes = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/ssmu/'
            self.save = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/'
        elif comp_node == 'saga':
            self.comp_node = 'saga'
            self.shapes = '/cluster/projects/nn9828k/Cian_sinmod/python/Krillmod/ssmu/'
            self.save = '/cluster/projects/nn9828k/Cian_sinmod/python/Krillmod/results/'
        else:
            print('Error: computation node not specified')
            sys.exit()

        # Define files based on folder setup
        # 1) Input files for analysis
        self.shp_file = self.shapes + self.shp_name
        self.trj_file = self.traj + self.sim + '/trajectory.nc'

        # 2) Intermediate save files
        self.poly_file = self.shapes + 'poly2grid.npy'
        self.time_file = self.save + self.sim + '/times.npy'
        self.depth_file = self.save + self.sim + '/depth.npy'

        return

    def exist(self):

        if not os.path.exists(self.shp_file):
            print('Error: Directory "' + self.shp_file + '" does not exist\nCheck connection to directory')
            sys.exit()
        else:
            print(self.shp_file + ' exists')

        if not os.path.exists(self.trj_file):
            print('Error: Directory "' + self.trj_file + '" does not exist\nCheck connection to directory')
            sys.exit()
        else:
            print(self.trj_file + ' exists')

        if not os.path.exists(self.save + self.sim):
            self.make_directory()
        else:
            print(self.save + self.sim + ' exists')

        if not os.path.exists(self.time_file):
            self.get_times()  # Creates time file in folder
        else:
            print(self.time_file + ' exists')

        if not os.path.exists(self.depth_file):
            self.get_depth()  # Creates depth file in folder
        else:
            print(self.depth_file + ' exists')

        if not os.path.exists(self.poly_file):
            self.get_poly()
        else:
            print(self.poly_file + ' exists')

        print('### NOTE:')
        print('All intermediate files exist for analysis')
        print('###')
        print('')
        return

    def make_directory(self):
        folder1 = self.save
        folder2 = folder1 + self.sim
        # Work around for remote server when we need to create two folders at once;
        if self.comp_node == 'local':
            os.mkdir(folder2)
        else:
            cmd = 'mkdir ' + folder1
            cmd2 = 'mkdir ' + folder2
            os.system(cmd)
            os.system(cmd2)
        return

    def get_poly(self):
        shape_v = self.read_shape()  # Extract shape data from file
        shp_p = np.shape(shape_v)  # Number of polygons
        dp = np.load(self.depth_file)
        shp_d = np.shape(dp)  # Shape of depth file
        file_inter = self.poly_file
        if not os.path.exists(file_inter):
            print('Note: Creating intermediate .npy file with polygon indices mapped to grid coordinates')
            store_ids = np.zeros([shp_d[0] * shp_d[1], 2])
            store_poly = np.zeros([shp_d[0] * shp_d[1], 1])
            c = 0
            for i in range(0, shp_d[0]):
                for j in range(0, shp_d[1]):
                    store_ids[c, 0] = j
                    store_ids[c, 1] = i
                    c = c + 1

            for p in range(0, shp_p[0]):
                print('poly = ' + str(p + 1) + ' of ' + str(shp_p[0]))
                poly1 = shape_v[p:p + 1].geometry  # Extract polygon geometry (gpd object)
                p_id = in_poly(store_ids[:, 0], store_ids[:, 1], poly1)
                store_poly[p_id == 1] = p + 1
            poly_ids = np.reshape(store_poly, [shp_d[0], shp_d[1]])
            np.save(file_inter, poly_ids)

    def read_shape(self):
        shape_p = gpd.read_file(self.shp_file)
        return shape_p

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



def in_poly(x, y, poly1):
    lat, lon = geo2grid(x, y, 'get_bl')  # Convert to geographic coordinates
    gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_xy(lon, lat))  # Create Geo dataframe from positional data
    polya = gpd.GeoDataFrame(poly1)  # Geo dataframe from polygon data
    gdf.crs = polya.crs  # Make sure we use the same projections for both
    in_pol = gpd.tools.sjoin(gdf, polya, predicate="within", how='left')  # Check points are in polygon
    in_idx = np.where(in_pol.index_right > -1)  # Find all values in the polygon
    p_id = np.zeros(np.shape(x))  # Initialize vector for storage
    p_id[in_idx[0]] = 1  # Store as ones
    return p_id



