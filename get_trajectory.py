import netCDF4 as nc
from netCDF4 import num2date
from netCDF4 import Dataset    # Note: python is case-sensitive!
import numpy as np
import math
import geopandas as gpd
import os


def store_regions(list_dir):
    # Shapes from trajectory file for creating nc dimensions
    store = store_traj(list_dir['traj_file'], 0)
    shape_v = read_shape(list_dir['shape_file'])  # Extract shape data from file
    shp_p = np.shape(shape_v)  # Number of polygons
    shp_t = store['t_id']  # Time steps
    shp_i = np.shape(store['active'])  # Number of individuals
    iv = np.shape(store['depth'])

    # Creates an intermediate matrix that stores polygon indices in each grid cell for fast identification of polygon
    # values for individuals
    file_inter = list_dir['shape_folder'] + str('poly2grid.npy')
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
            p_id = in_poly(store_ids[:, 0], store_ids[:, 1], poly1)
            store_poly[p_id == 1] = p + 1
        poly_ids = np.reshape(store_poly, [iv[0], iv[1]])
        np.save(file_inter, poly_ids)

    # Create nc file with the transpose of the columns
    poly_ids = np.load(file_inter)
    act_poly = np.empty(shp_i[0], dtype=np.int16)
    in_area = np.empty([shp_i[0], shp_t[0]], dtype=np.int16)
    x_vals = np.empty([shp_i[0], shp_t[0]], dtype=np.int16)
    y_vals = np.empty([shp_i[0], shp_t[0]], dtype=np.int16)
    z_vals = np.empty([shp_i[0], shp_t[0]], dtype=np.int16)

    for t in range(0, shp_t[0]):
        store_t = store_traj(list_dir['traj_file'], t)  # Time slice of trajectory
        x = store_t['xp'].astype(int)
        y = store_t['yp'].astype(int)
        z = store_t['zp'].astype(int)
        act = store_t['active'].astype(int)  # Note: I don't deactivate individuals, Ingrid suggests to do this;
        act_poly = act_poly + act
        poly_fid = poly_ids[y, x]
        idx_pid = (np.where(poly_fid > 0))
        in_area[idx_pid, t] = poly_fid[idx_pid]
        x_vals[:, t] = x
        y_vals[:, t] = y
        z_vals[:, t] = z
        print('t = ' + str(t) + ' of ' + str(shp_t[0]) + ' steps')
        print('Percent complete = ' + str(np.ceil((t/shp_t[0])*100)))

    act_poly = shp_t[0] - act_poly  # Index of activity
    list_start = np.unique(act_poly).astype(int)  # find starting points of each individual
    start_point = np.zeros(shp_i[0])
    for i in range(0, len(list_start)):
        id1 = list_start[i]
        log_id1 = act_poly == id1
        in_polt = in_area[log_id1, id1:shp_t[0]]
        if np.shape(in_polt)[0] <= 0:
            print('Empty start areas ' + str(np.shape(in_polt)))
        else:
            start_point[log_id1] = in_polt[:, 0]

    nc_file = Dataset(list_dir['reg_file'], mode='w', format='NETCDF4_CLASSIC')

    # Specify nc dimensions
    part_dim = nc_file.createDimension('particle', shp_i[0])
    time_dim = nc_file.createDimension('time', shp_t[0])

    # Create variable for storing presence/ absence in each region at each time:
    act_ind = nc_file.createVariable('act_part', np.int16, 'particle')
    start_area = nc_file.createVariable('start', np.int16, 'particle')
    in_region = nc_file.createVariable('in_region', np.int16, ('particle', 'time'))
    xv = nc_file.createVariable('xp', np.int16, ('particle', 'time'))
    yv = nc_file.createVariable('yp', np.int16, ('particle', 'time'))
    zv = nc_file.createVariable('zp', np.int16, ('particle', 'time'))

    # Store data in nc variables
    act_ind[:] = act_poly
    xv[:] = x_vals
    yv[:] = y_vals
    zv[:] = z_vals
    in_region[:] = in_area
    start_area[:] = start_point

    print(nc_file)
    nc_file.close()
    print('Dataset1 is closed!')
    return


def get_times(list_dir):
    traj = nc.Dataset(list_dir['traj_file'])
    times = traj.variables['time']
    date_save = num2date(times, times.units)
    np.save(list_dir['time_file'], date_save.data)
    return


def get_depth(list_dir):
    traj = nc.Dataset(list_dir['traj_file'])
    depth = traj.variables['depth'][:]
    dp = depth[:]  # get depth matrix
    idx = dp < 0  # Fill values
    dp[idx] = np.nan  # Set fill values to invalid;
    dp2 = np.array(dp)
    np.save(list_dir['depth_file'], dp2)
    return


def store_traj(tr_file, idv):
    # netcdf4 library for extracting dataset
    traj = nc.Dataset(tr_file)
    depth = traj.variables['depth'][:]
    dp = depth[:]  # get depth matrix
    dp_size = np.shape(dp)  # Size of domain
    idx = dp < 0  # Fill values
    dp[idx] = np.nan  # Set fill values to invalid;
    imax = dp_size[0]
    jmax = dp_size[1]
    times = traj.variables['time']
    t_ids = np.shape(times)
    dates = num2date(times[idv], times.units)
    x, y, z = (np.array(traj.variables['x'][idv, :]), np.array(traj.variables['y'][idv, :]),
               np.array(traj.variables['z'][idv, :]))
    act = np.array(traj.variables['active'][idv, :])

    # Store outputs in a dictionary (similar to a structure)
    store = dict([('imax', imax), ('jmax', jmax), ('depth', dp), ('active', act), ('xp', x), ('yp', y),
                  ('zp', z), ('time', dates), ('t_id', t_ids)])
    traj.close()
    return store


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


def read_shape(shp_file):
    shape_p = gpd.read_file(shp_file)
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