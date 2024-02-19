from get_trajectory import store_regions, get_depth, get_times
import os
import sys


def locate_folders(comp_node, sim_folder, shape_name):
    # Name directories here
    if comp_node == 'local':
        shape_folder = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/ssmu/'  # Directory where the shape files are stored
        traj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # Directory where the trajectory file is stored
        save_folder = 'C:/Users/ciank/PycharmProjects/sinmod/Krillmod/results/'  # Directory where data will be saved
        list_dir = init_folders(comp_node, sim_folder, shape_name, shape_folder, save_folder, traj_folder)
    elif comp_node == 'saga':
        shape_folder = '/cluster/projects/nn9828k/Cian_sinmod/python/Krillmod/ssmu/'  # Directory where the shape files are stored
        traj_folder = '/cluster/projects/nn9828k/Cian_sinmod/meeso_sim/sim_'  # Directory where the trajectory file is stored
        save_folder = '/cluster/projects/nn9828k/Cian_sinmod/python/Krillmod/results/'  # Directory where data will be saved
        list_dir = init_folders(comp_node, sim_folder, shape_name, shape_folder, save_folder, traj_folder)
    else:
        print('Error: computation node not specified')
        sys.exit()
    return list_dir


def init_folders(comp_node, sim_folder, shape_name, shape_folder, save_folder, traj_folder):
    # Stores directories in a dictionary
    list_dir = dict([('shape_folder', shape_folder), ('traj_folder', traj_folder), ('save_folder', save_folder),
                     ('sim_folder', sim_folder), ('shape_name', shape_name)])

    # Adds relevant files to dictionary. The regions.nc file is derived from trajectory.nc and shape files,
    # storing information on polygon index of individual for each time step
    list_dir['shape_file'] = list_dir['shape_folder'] + list_dir['shape_name']
    list_dir['traj_file'] = list_dir['traj_folder'] + list_dir['sim_folder'] + '/trajectory.nc'
    list_dir['reg_file'] = list_dir['save_folder'] + list_dir['sim_folder'] + '/regions.nc'
    list_dir['time_file'] = list_dir['save_folder'] + list_dir['sim_folder'] + '/times.npy'
    list_dir['depth_file'] = list_dir['save_folder'] + list_dir['sim_folder'] + '/depth.npy'

    if not os.path.exists(list_dir['shape_file']):
        print('Error: Directory "' + list_dir['shape_file'] + '" does not exist\nCheck connection to directory')
        sys.exit()

    if not os.path.exists(list_dir['traj_file']):
        print('Error: Directory "' + list_dir['traj_file'] + '" does not exist\nCheck connection to directory')
        sys.exit()

    if not os.path.exists(list_dir['save_folder'] + list_dir['sim_folder']):
        if comp_node == 'local':
            os.mkdir(list_dir['save_folder'] + list_dir['sim_folder'])
        else:
            os.system('mkdir' + list_dir['save_folder'])
            os.system('mkdir' + list_dir['save_folder'] + list_dir['sim_folder'])

    if not os.path.exists(list_dir['time_file']):
        get_times(list_dir)  # Creates time file in folder

    if not os.path.exists(list_dir['depth_file']):
        get_depth(list_dir)  # Creates depth file in folder

    if not (os.path.exists(list_dir['reg_file'])):
        # Use try os.system(echo) or something of the sort  here;
        print('Note: Reformatting trajectory data for analysis')
        store_regions(list_dir)  # Creates intermediate file with trajectory and regional data

    return list_dir


