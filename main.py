# master file for krill model analysis
from sim_scenarios import main_analysis, compare_2_years, split_analysis, get_quant  # simulation scenarios

# sim_list = ['2017']  # Simulation identifier
# trj_folder = 'D:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
# sim_list = ['DVM']
# trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
# shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
# get_quant(sim_list, trj_folder, shp_name, node='local')
#
# breakpoint()
# 1: Main analysis for a list of simulations;
sim_list = ['2017']  # Simulation identifier
trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
#trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
main_analysis(sim_list, trj_folder, shp_name, node='local')



split_analysis(sim_list, trj_folder, shp_name, node='local')
breakpoint()



sim_list = ['2017']  # Simulation identifier
trj_folder = 'A:/Cian_sinmod/meeso_sim/sim_'  # trajectory folder
#trj_folder = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
shp_name = 'ssmusPolygon.shp'  # Name of shape polygon if relevant
main_analysis(sim_list, trj_folder, shp_name, node='local')

breakpoint()

# 2: Comparisons between 2 simulations
sim_folder1 = 'DVM'  # Simulation identifier
sim_folder2 = '2017'  # Second simulation identifier
t_folder1 = 'D:/Cian_sinmod/sim_'  # trajectory folder- betzy
t_folder2 = 'A:/Cian_sinmod/sim_'  # trajectory folder- saga
compare_2_years(t_folder1, sim_folder1, t_folder2, sim_folder2, node='local')


