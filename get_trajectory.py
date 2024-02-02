import netCDF4 as nc
from netCDF4 import num2date
from netCDF4 import Dataset    # Note: python is case-sensitive!
import numpy as np
import math
import shapely.geometry
import pandas as pd
import geopandas as gpd
import pickle
import os

# from Krillmod.import_settings import tr_file


# from Krillmod.import_settings import sv_dir


def store_traj(file, idv):
    # Importing libraries: Problem with code- not closing netcdf files here;
    # get dataset
    traj = nc.Dataset(file)
    # print(traj) #similar to ncdump -h
    # print(traj.variables.keys()) #prints list of variables
    # Subset variables:
    depth = traj.variables['depth'][:]
    dp = depth[:]  # get depth matrix
    dp_size = np.shape(dp)  # Size of domain
    idx = dp < 0  # Fill values
    dp[idx] = np.nan  # Set fill values to invalid;
    imax = dp_size[0]
    jmax = dp_size[1]
    times = traj.variables['time']
    dates = num2date(times[:], times.units)
    x, y, z = (np.array(traj.variables['x'][idv, :]), np.array(traj.variables['y'][idv, :]),
               np.array(traj.variables['z'][idv, :]))
    act = np.array(traj.variables['active'][idv, :])

    #x_act = x[act == 1]
    #y_act = y[act == 1]
    #z_act = z[act == 1]

    # Store in dictionary
    store = dict([('imax', imax), ('jmax', jmax), ('depth', dp), ('active', act), ('xp', x), ('yp', y),
                  ('zp', z), ('time', dates)])
    traj.close()
    return store


def store_regions(file, idv):
    nc_file = nc.Dataset(file, mode='r', format='NETCDF4_CLASSIC')
    # nc_file.variables.keys()
    in_polt = np.array(nc_file.variables["in_region"][idv, :])
    # store_pol = dict([("in_reg", in_polt)])
    nc_file.close()
    # Make a check for the year I am trying to read;
    return in_polt


def get_regions(shape_v):
    store = store_traj(tr_file, 0)
    # Shapes from trajectory file for creating nc dimensions
    shp_p = np.shape(shape_v)  # Number of polygons
    shp_t = np.shape(store['time'])  # Time steps
    shp_i = np.shape(store['active'])  # Number of individuals
    iv = np.shape(store['depth'])

    from Krillmod.import_settings import shp_dir
    file_inter = shp_dir + str('poly2grid.npy')
    if not os.path.exists(file_inter):
        print('Note: Creating intermediate .npy file with polygon indices mapped to grid coordinates')
        store_ids = np.zeros([iv[0]*iv[1], 2])
        store_poly = np.zeros([iv[0]*iv[1], 1])
        c = 0
        for i in range(0, iv[0]):
            for j in range(0, iv[1]):
                store_ids[c, 0] = j
                store_ids[c, 1] = i
                c = c + 1

        for p in range(0, shp_p[0]):
            print('poly = ' + str(p+1) + ' of ' + str(shp_p[0]))
            poly1 = shape_v[p:p + 1].geometry  # Extract polygon geometry (gpd object)
            p_id = inpoly(store_ids[:, 0], store_ids[:, 1], poly1)
            store_poly[p_id == 1] = p + 1
        poly_ids = np.reshape(store_poly, [iv[0], iv[1]])
        np.save(file_inter, poly_ids)

    poly_ids = np.load(file_inter)
    temp_poly = np.empty([shp_t[0], shp_i[1]])
    for t in range(0, shp_t[0]):
        store_t = store_traj(tr_file, t)  # Time slice of trajectory
        x = store_t['xp'].astype(int)
        y = store_t['yp'].astype(int)
        poly_fid = poly_ids[y, x]
        idx_pid = (np.where(poly_fid > 0))
        temp_poly[t, idx_pid] = poly_fid[idx_pid]
        print('t = ' + str(t) + ' of ' + str(shp_t[0]) + ' steps')
        print('Percent complete = ' + str(np.ceil((t/shp_t[0])*100)))

    # Check that file isn't already open
    try:
        nc_file.close()  # just to be safe, make sure dataset is not already open.
    except:
        pass

    # Create filepath
    from Krillmod.import_settings import sv_dir, tr_dir
    # r_save = sv_dir + 'regions.nc'
    r_save = tr_dir + 'regions.nc'
    nc_file = Dataset(r_save, mode='w', format='NETCDF4_CLASSIC')

    # Specify nc dimensions
    part_dim = nc_file.createDimension('particle', shp_i[1])
    time_dim = nc_file.createDimension('time', shp_t[0])
    area_dim = nc_file.createDimension('polygon', shp_p[0])

    # for dim in nc_file.dimensions.items():
    #     print(dim)

    # Create variable for storing presence/ absence in each region at each time:
    in_region = nc_file.createVariable('in_region', np.int32, ('time', 'particle'))  # , 'polygon'))
    # consider adding variables: active part; time; polygon names;

    in_region[:] = temp_poly
    print(nc_file)
    nc_file.close()
    print('Dataset is closed!')
    return


def inpoly(x, y, poly1):
    lat, lon = geo2grid(x, y, 'get_bl')  # Convert to geographic coordinates
    # Create geodataframe from numpy arrays
    gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_xy(lon, lat))
    polya = gpd.GeoDataFrame(poly1)
    gdf.crs = polya.crs
    in_pol = gpd.tools.sjoin(gdf, polya, predicate="within", how='left')
    in_idx = np.where(in_pol.index_right > -1)
    p_id = np.zeros(np.shape(x))
    p_id[in_idx[0]] = 1
    return p_id


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
