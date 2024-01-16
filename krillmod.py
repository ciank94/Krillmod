# master file for krill model analysis
from Krillmod.getTrajectory import getTrajectory

#Function that outputs trajectory.nc file into a dictionary storing arrays, cftime and integers
file = 'C:/Users/ciank/OneDrive - NTNU/PostDoc/d_Data visualization/trajectory.nc'
store = getTrajectory(file)
