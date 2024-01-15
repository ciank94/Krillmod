#Importing libraries
import pandas as pd
import netCDF4 as nc
from netCDF4 import num2date, date2num, date2index

#Specify directories
filepath = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'

#get dataset
traj = nc.Dataset(filepath)
#print(traj) #similar to ncdump -h
#print(traj.variables.keys())

#Subset variables:
times = traj.variables['time']
dates = num2date(times[:], times.units)
print(dates)
print([date.strftime('%Y-%m-%d %H:%M:%S') for date in dates[:10]]) # print only first ten...
#strftime('%d-%m-%Y %H:%M:%S')
#t = traj.variables['time']
#time = t[:]
#print(time)
#xp = traj.variables['x']
#print(xp)
#print(xp.head(10))
#print(ds.variables.keys())