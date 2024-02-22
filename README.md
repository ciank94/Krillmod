# Krillmod
#### Short description of project: 
Description of moduels used to analyse trajectory data stored in netcdf format (``trajectory.nc``). Firstly, the modules
are briefly summarised. Secondly, the main ``krillmod.py`` file is described with example output and plots. Thirdly, I 
go into more detail on how the individual modules work in the background.
## Modules:
``krillmod.py``: Main python file for calling modules used in the analysis of trajectory.nc data.\
``import_settings.py``: Stores directories and file names in a dictionary, with settings for local and remote directories.\
``get_trajectory.py``: Functions for accessing and pre-processing trajectory netcdf data stored with CF conventions.\
``analyse_trajectory.py``: Functions for analysing trajectory data accessed with ``get_trajectory.py``.\
``plot_trajectory.py``: Module with functionality for plotting output from ``analyse_trajectory.py``.

### Module ``krillmod.py``:
#todo: Insert description here with example output for each step; \
#1: import settings: print keys and what it points to as in sinmodules read_average_file; \
#2: get_trajectory: regions.nc file key for which polygon an individual was in obs x time; with printout for a 
sample(text file as in import settings.)\
#3: analyse_trajectory: sample output.\
#4: plot_trajectory: nicest figures.

### Module ``import_settings.py``:
* Function ``locate_folders`` accepts arguments for specifying whether the computation node is local or 
remote (``comp_node``), the name of the simulation folder (``sim_folder``) and the ``.shp`` shapefile used for 
subsets of individuals (``shape_name``). These inputs as well as directory specifications in locate_folders function, 
are used to create a dictionary (``list_dir``) of directories for ease of usage:

    ```
    list_dir = locate_folders(comp_node, sim_folder, shape_name)
    ```
* Within ``locate_folders``, the directory for saving analysis output (``save_folder``), containing the trajectory.nc 
file (``traj_folder``) and folder with shape files (``shape_file``) should be specified by the user, depending on 
whether the script is being run locally or on a remote server. Note that the shape folder should also contain a ``.shx`` 
file necessary for reading many ``.shp`` files.
* Finally, the trajectory file (``trajectory.nc``), regional file (intermediate netcdf file named ``regions.nc`` 
described below), a file with simulation timestamps (``time.npy``) and a file with bottom depth data (``depth.npy``) are
all added to ``list_dir`` or first, created with ``os.mkdir()`` remotely or ``os.system('mkdir ')`` on a linux server.
* The directories above are stored in the  for ease of usage:

### Module ``get_trajectory.py``:
* The function ``store_traj`` takes in file path as input and stores
trajectory data in netcdf format ``store_traj(file)``:
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




