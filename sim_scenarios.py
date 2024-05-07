# File that contains different simulation comparisons for specific projects e.g. comparing two years;
def main_analysis(sim_list, trj_folder, shp_name, node):
    from set_directory import Folder
    from get_trajectory import Trajectory, Regional
    from analyse_trajectory import Analyse
    from plot_trajectory import Plots
    for s in sim_list:
        # Define names of main folders
        sim_folder = s  # Name of simulation instance

        # Setup directory information for analysis
        f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
        f.set_folder(node)  # Sets folders based on remote or local node
        f.exist()  # Checks existence of relevant folders and files

        # Two classes for dealing with data:
        #k_r = Regional(f)
        k_t = Trajectory(f)

        import numpy as np
        import matplotlib.pyplot as plt
        shp_t = np.shape(k_t.x)[0]
        batch_x = np.zeros([30, k_t.n_sites])
        batch_y = np.zeros([30, k_t.n_sites])
        SG_x = 520
        SG_y = 650
        test_mat = np.zeros(k_t.n_sites)
        for j in range(0, k_t.n_sites):
            print(str(j))
            for i in range(0, shp_t, 100):
                x_v = k_t.x[i, j]
                y_v = k_t.y[i, j]
                e_dist = np.sqrt(((x_v - SG_x)**2) + (y_v - SG_y)**2)
                if e_dist < 100:
                    print('close to SG')
                    test_mat[j] = i
                    break
        breakpoint()
        #[a1,b1,c1] = k_t.get_slice(0)
        breakpoint()


        # Now the analysis
        df = Analyse(f)
        df.dom_paths(k_r)
        #df.depth_profile(k_t)  # Note
        df.retention_times(k_r)
        df.transit_times('ALL', k_r)

        p = Plots()
        p.plot_connectivity(f, k_r)
        p.plot_transit(f)
        p.plot_retention(f)
        p.plot_dom_paths(f, k_r)

        #k_t.nc_file.close()
        k_r.nc_file.close()
    return

def split_analysis(sim_list, trj_folder, shp_name, node):
    from set_directory import Folder
    from get_trajectory import Trajectory, Regional
    from analyse_trajectory import Analyse
    from plot_trajectory import Plots
    for s in sim_list:
        # Define names of main folders
        sim_folder = s  # Name of simulation instance

        # Setup directory information for analysis
        f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
        f.set_folder(node)  # Sets folders based on remote or local node
        f.exist()  # Checks existence of relevant folders and files

        # Two classes for dealing with data:
        k_r = Regional(f)
        breakpoint()
        # k_t = Trajectory(f)

        # Now the analysis
        # df = Analyse(f)
        # df.dom_paths(k_r)
        # # df.depth_profile(k_t)  # Note
        # df.retention_times(k_r)
        # df.transit_times('ALL', k_r)
        #
        # p = Plots()
        # p.plot_connectivity(f, k_r)
        # p.plot_transit(f)
        # p.plot_retention(f)
        # p.plot_dom_paths(f, k_r)

        # k_t.nc_file.close()
        k_r.nc_file.close()



def get_quant(sim_list, trj_folder, shp_name, node):
    from set_directory import Folder
    from plot_trajectory import Plots
    for s in sim_list:
        # Define names of main folders
        sim_folder = s  # Name of simulation instance

        # Setup directory information for analysis
        f = Folder(trj_folder, sim_folder, shp_name)  # Initialize folders
        f.set_folder(node)  # Sets folders based on remote or local node

        df = Plots()
        df.save_data(f)




def compare_2_years(t_folder1, sim_folder1, t_folder2, sim_folder2, node):
    from set_directory import FolderComp
    from analyse_trajectory import AnalyseComp
    from plot_trajectory import Plots

    files = FolderComp(t_folder1, sim_folder1, t_folder2, sim_folder2)
    files.set_folder(node)
    files.exist()

    df = AnalyseComp(files)
    df.compare_pathways()
    df.compare_retention()

    p = Plots()
    p.plot_comp_paths(files)
    p.plot_comp_retention(files)
    return
