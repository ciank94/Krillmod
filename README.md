# Krillmod
## Modules:
``krillmod.py``: Main python file for calling modules.\
``get_trajectory.py``: Module containing functions for accessing and pre-processing netcdf data stored based on CF conventions.\
``analyse_trajectory.py``: Module with functions for analysing trajectory data accessed with ``get_trajectory.py``.\
``plot_trajectory.py``: Module with functionality for plotting output from ``analyse_trajectory.py``.

### Module ``get_trajectory.py``:
The function ``store_traj`` takes in file path as input and stores
trajectory data in netcdf format ``store_traj(file)``
```python
file = 'A:/sim_2016/trajectory.nc' #Example filepath
store = store_traj(file) #Store netcdf data using dictionary keys
print(store.keys()) #Print keys defined in store_traj
```
_Output:_
```python
dict_keys(['imax', 'jmax', 'depth', 'xp', 'yp', 'zp', 'time'])
```
Each key can be accessed using ``store[key]``, for example, 
``store['xp']`` accesses the x positions.  The ``['time']`` key accesses times in cftime format, where the ``trajectory.nc`` file is
stored with units 'days since', for example: ``units: days since 2016-01-01 00:00:00``. For more information on accessing data from netcdf files using netcdf4, see examples: 
https://github.com/Unidata/netcdf4-python/blob/master/examples/reading_netCDF.ipynb.

We create one large netcdf file regions.nc which transposes original file to particles x time (p x t) and 
adds vector with indexes for time and region they begin

Function ``single_traj`` accesses a single time slice of the trajectories using the store dictionary defined in store_traj above.
extracts a single trajectory slice in time based on the index ``idx`` [TODO: Intermediate step: find idx based on time range]:
```python
slice = single_traj(store, idx)
print(slice['xi'])
```
_Output:_
```python
[103.23898 121.04438 123.98057 ... 445.      445.      445.]
```

Function ``read_ssmu`` is used for accessing polygon data from shapefiles using geopandas. It stores information on subunits file  in the Southern Ocean defining important regions for krill predators. It assumes the 
``.shp`` file is in the same folder as the ``.shx`` file required. Shapefiles courtesy of https://github.com/ccamlr/data.
```python
shp_file = 'A:/ssmu/ssmusPolygon.shp'
shape = read_ssmu(shp_file)
print(shape[:5])
```
_Output:_
```
GAR_ID  ...                                           geometry
0   95561  ...  POLYGON ((-50.00000 -64.00000, -50.00000 -65.0...
1   95562  ...  POLYGON ((-61.11353 -64.32835, -61.11216 -64.3...
2   95563  ...  POLYGON ((-60.53607 -61.73011, -60.53600 -61.7...
3   95564  ...  POLYGON ((-56.37107 -61.54253, -56.38300 -61.5...
4   95588  ...  POLYGON ((-58.24770 -63.45529, -58.25221 -63.4...

[5 rows x 15 columns]
```




