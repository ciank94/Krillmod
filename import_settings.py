from Krillmod.get_trajectory import store_regions, get_depth, get_times
import os
import sys


def store_folders(sim_folder, shape_name, shape_folder, save_folder, traj_folder):
    # Stores relevant folders in dictionary
    list_dir = dict([('shape_folder', shape_folder), ('traj_folder', traj_folder), ('save_folder', save_folder),
                    ('sim_folder', sim_folder), ('shape_name', shape_name)])
    return list_dir


def store_files(list_dir):
    # Adds relevant files to dictionary. The regions.nc file is derived from trajectory.nc and shape files,
    # storing information on polygon index of individual for each time step
    list_dir['shape_file'] = list_dir['shape_folder'] + list_dir['shape_name']
    list_dir['traj_file'] = list_dir['traj_folder'] + list_dir['sim_folder'] + '/trajectory.nc'
    list_dir['reg_file'] = list_dir['save_folder'] + list_dir['sim_folder'] + '/regions.nc'
    list_dir['time_file'] = list_dir['save_folder'] + list_dir['sim_folder'] + '/times.npy'
    list_dir['depth_file '] = list_dir['save_folder'] + list_dir['sim_folder'] + '/depth.npy'

    if not os.path.exists(list_dir['shape_file']):
        print('Error: Directory "' + list_dir['shape_file'] + '" does not exist')
        sys.exit()

    if not os.path.exists(list_dir['traj_file']):
        print('Error: Directory "' + list_dir['traj_file'] + '" does not exist')
        sys.exit()

    if not os.path.exists(list_dir['save_folder'] + list_dir['sim_folder']):
        os.mkdir(list_dir['save_folder'] + list_dir['sim_folder'])

    if not os.path.exists(list_dir['time_file']):
        get_times(list_dir)

    if not os.path.exists(list_dir['depth_file']):
        get_depth(list_dir)

    if not os.path.exists(list_dir['reg_file']):
        print('Note: Reformatting trajectory data for analysis')
        store_regions(list_dir)

    return list_dir
