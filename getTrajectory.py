def getTrajectory(file):
    # Importing libraries
    import pandas as pd
    import netCDF4 as nc
    import matplotlib.pyplot as plt
    import numpy as np
    from netCDF4 import num2date, date2num, date2index

    # get dataset
    traj = nc.Dataset(file)
    # print(traj) #similar to ncdump -h
    # print(traj.variables.keys()) #prints list of variables

    # Subset variables:
    depth = traj.variables['depth']
    x, y = traj.variables['x'], traj.variables['y']
    dp = depth[:] # get depth matrix
    dp_size = np.shape(dp) # Size of domain
    idx = dp < 0 # Fill values
    dp[idx] = np.nan # Set fill values to invalid;
    imax = dp_size[0]
    jmax = dp_size[1]
    times = traj.variables['time']
    dates = num2date(times[:], times.units)

    # Store in dictionary
    store = dict([('imax', imax), (jmax, 'jmax'), ('depth', dp),
                  ('xp', x), ('yp', y), ('time', dates)])
    return store



