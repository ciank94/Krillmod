import netCDF4 as nc
from netCDF4 import num2date
import numpy as np
import math
import geopandas as gpd


def store_traj(file):
    # Importing libraries
    # get dataset
    traj = nc.Dataset(file)
    # print(traj) #similar to ncdump -h
    # print(traj.variables.keys()) #prints list of variables
    # Subset variables:
    depth = traj.variables['depth']
    dp = depth[:]  # get depth matrix
    dp_size = np.shape(dp)  # Size of domain
    idx = dp < 0  # Fill values
    dp[idx] = np.nan  # Set fill values to invalid;
    imax = dp_size[0]
    jmax = dp_size[1]
    times = traj.variables['time']
    dates = num2date(times[:], times.units)
    x, y, z = traj.variables['x'], traj.variables['y'], traj.variables['z']
    act = traj.variables['active']

    # Store in dictionary
    store = dict([('imax', imax), ('jmax', jmax), ('depth', dp), ('active', act),
                  ('xp', x), ('yp', y), ('zp', z), ('time', dates)])
    # traj.close()
    return store


def single_traj(store, idx):
    # Importing libraries
    # get dataset
    x = np.array(store['xp'][idx:idx + 1][0])
    y = np.array(store['yp'][idx:idx + 1][0])
    active = np.array(store['active'][idx:idx + 1][0])
    xi = x[active == 1]
    yi = y[active == 1]
    idx = x <= 0  # Fill values
    # x[idx] = np.nan  # Set fill values to invalid;
    # y[idx] = np.nan
    slice_t = dict([('xp', xi), ('yp', yi)])
    return slice_t


def read_ssmu(shp_file):

    shape_p = gpd.read_file(shp_file)
    # print(shape.boundary)
    # pl = shape.plot()
    # pl.imshow()
    #
    # print(shape[:17])
    return shape_p


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

# def region_part()
#     #% !(AP: 1:APPA; 2: APW; 3: DPW; 4: DPE; 5: BSW; 6:BSE; 7: EI; 17: APE)
#     #% !(SOI: 8: SOPA; 9: SOW; 10:SONE; 11: SOSE)
#     #% !(SG: 12: SGPA; 13: SGW; 14:SGE)
#     #% !(SSI: 15: SSPA; 16: SSI)
#     return
#
