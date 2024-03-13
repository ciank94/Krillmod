import netCDF4 as nc
from netCDF4 import num2date
import numpy as np
import geopandas as gpd
import math
import os


class Trajectory:

    def __init__(self, f):
        # From files initialized in f
        self.poly = np.load(f.poly_file)
        self.nc_file = nc.Dataset(f.trj_file)
        self.depth = np.load(f.depth_file)
        self.time = np.load(f.time_file, allow_pickle=True)

        # From the trajectory file
        self.x = self.nc_file['x']
        self.y = self.nc_file['y']
        self.z = self.nc_file['z']
        self.t_max = np.shape(self.x)[0]
        self.i_max = np.shape(self.depth)[0]
        self.j_max = np.shape(self.depth)[1]
        self.active = self.nc_file['active']

        # Derived data, batches and number of batches, when batches arrive
        self.n_sites = np.sum(self.active[0, :])
        self.n_parts = np.shape(self.x)[1]
        self.n_batches = np.floor(self.n_parts/self.n_sites).astype(int)
        self.sim_account(f)
        # self.nc_file.close()
        return

    def get_batch(self, i):
        # Method to subset individuals
        start_i = (i-0)*self.n_sites
        end_i = (i+1)*self.n_sites
        slice_i = slice(start_i, end_i, 1)
        act_t = self.t_max - np.sum(self.active[:, start_i])
        slice_t = slice(act_t, self.t_max, 1)
        return slice_i, slice_t

    def in_region(self, x, y):
        # Function for finding polygon of individuals from polygrid.npy
        p_id = self.poly[x, y]
        return p_id

    def subset_data(self, id_start, id_end):
        self.z = self.z[id_start:id_end, :]
        self.x = self.x[id_start:id_end, :]
        self.y = self.y[id_start:id_end, :]
        return




    def sim_account(self, f):
        # Note: Adapt this file for  new data storage
        # Access regional file and look at number of particles, number of batches, individuals in each batch
        summary_file = f.save + f.sim + '/sim_summary.txt'
        if not os.path.exists(summary_file):
            print('')
            print('Saving: ' + summary_file)
            time = self.time
            shp = np.shape(time)
            date_0 = time[0]
            date_1 = time[shp[0] - 1]
            timeframe = date_1 - date_0
            start_date = 'Start datetime (dd.mm.yyyy) = ' + date_0.strftime("%d.%m.%Y, %H:%M:%S")
            end_date = 'End datetime (dd.mm.yyyy) = ' + date_1.strftime("%d.%m.%Y, %H:%M:%S")
            time_frame = 'Time frame = ' + str(timeframe.days) + ' days and ' + str(
                timeframe.seconds * 1 / (60 * 60)) + ' hours '
            f = open(summary_file, 'w')
            top_line = ['Simulation summary: \n']
            time_lines = [start_date + '\n', end_date + '\n', time_frame + '\n']
            n_parts_line0 = ['Saved every: ' + str(np.round(timeframe.days/self.t_max, 2)) + ' days ' +  '\n']
            n_parts_line1 = ['Total number of particles: ' + str(self.n_parts) + '\n']
            n_parts_line2 = ['Number of batches: ' + str(self.n_batches) + '\n']
            n_parts_line3 = ['Number of sites: ' + str(self.n_sites) + '\n']
            n_parts_line4 = ['Time steps: ' + str(self.t_max) + '\n']
            n_parts_line5 = ['i_max: ' + str(self.i_max) + '\n']
            n_parts_line6 = ['j_max: ' + str(self.j_max) + '\n']
            f.writelines(top_line)
            f.writelines(time_lines)
            f.writelines(n_parts_line0)
            f.writelines(n_parts_line1)
            f.writelines(n_parts_line2)
            f.writelines(n_parts_line3)
            f.writelines(n_parts_line4)
            f.writelines(n_parts_line5)
            f.writelines(n_parts_line6)

            f.close()
        else:
            print('')
            print('Directory: ' + summary_file + ' already exists, skipping')
        return


def geo2grid(lat, lon, case):
    case_types = ['get_xy', 'get_bl']
    if case not in case_types:
        raise ValueError("Invalid case type in 3rd argument. Expected one of: %s" % case_types)

    a = 6378206.4  # Earth Radius
    fm = 294.97870  # Inverse flattening
    f = 1 / fm  # Flattening
    e = math.sqrt((2 * f) - (f ** 2))  # Eccentricity
    lon_0 = -45.0  # False origin longitude
    lat_0 = -44.0  # False origin latitude
    lat_1 = -40  # First parallel
    lat_2 = -68  # Second parallel
    x_0 = 376  # Easting at false origin
    y_0 = 971  # Northing at false origin
    dx = 4
    imax = 825  # Weddell Sea domain

    rval = len(lat)
    xs = np.empty(rval)
    ys = np.empty(rval)
    for i in range(rval):
        FiF = lat_0 * math.pi / 180
        Fi1 = lat_1 * math.pi / 180
        Fi2 = lat_2 * math.pi / 180
        LamdaF = lon_0 * math.pi / 180

        if case == 'get_xy':
            Fi = lat[i] * math.pi / 180
            Lamda = lon[i] * math.pi / 180

        EF = x_0 * dx * 1000
        NF = y_0 * dx * 1000

        if case == 'get_bl':
            E = lat[i] * dx * 1000
            N = lon[i] * dx * 1000

        m1 = math.cos(Fi1) / math.sqrt(1 - e ** 2 * (math.sin(Fi1)) ** 2)
        m2 = math.cos(Fi2) / math.sqrt(1 - e ** 2 * (math.sin(Fi2)) ** 2)

        t1 = math.tan(math.pi / 4 - Fi1 / 2) / (((1 - e * math.sin(Fi1)) / (1 + e * math.sin(Fi1))) ** (e / 2))
        t2 = math.tan(math.pi / 4 - Fi2 / 2) / (((1 - e * math.sin(Fi2)) / (1 + e * math.sin(Fi2))) ** (e / 2))
        tF = math.tan(math.pi / 4 - FiF / 2) / (((1 - e * math.sin(FiF)) / (1 + e * math.sin(FiF))) ** (e / 2))

        if case == 'get_xy':
            t = math.tan(math.pi / 4 - Fi / 2) / (((1 - e * math.sin(Fi)) / (1 + e * math.sin(Fi))) ** (e / 2))

        n = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
        F = m1 / (n * t1 ** n)
        rF = a * F * tF ** n

        if case == 'get_xy':
            r = a * F * t ** n

        if case == 'get_bl':
            rm = np.sign(n) * np.sqrt((E - EF) ** 2 + (rF - (N - NF)) ** 2)
            tm = (rm / (a * F)) ** (1 / n)
            Tetam = np.arctan((E - EF) / (rF - (N - NF)))
            Fim = np.pi / 2 - 2 * np.arctan(Tetam)
            for j in range(1, 9):

                Fi = np.pi / 2 - 2 * np.arctan(tm * ((1 - e * np.sin(Fim)) / (1 + e * np.sin(Fim))) ** (e / 2))
                if np.abs(Fi - Fim) < 1e-7:
                    break
                else:
                    Fim = Fi
            Lamda = Tetam / n + LamdaF
            xs[i] = Fi * 180 / math.pi
            ys[i] = Lamda * 180 / math.pi

        if case == 'get_xy':
            Teta = n * (Lamda - LamdaF)
            x = EF + r * math.sin(Teta)
            y = NF + rF - r * math.cos(Teta)
            xs[i] = x / (dx * 1000)
            ys[i] = y / (dx * 1000)

    return xs, ys