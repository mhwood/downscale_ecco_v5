
import os
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
from datetime import datetime, timedelta


########################################################################################################################

def create_beaufort_dataset(config_dir):

    print('Creating the beaufort dv reference dictionary')

    output_period = 3600
    timestep = 120
    iters_per_output = int(output_period/timestep)
    test_var_name = 'AREA'
    Nr = 1
    N = 207104

    dv_dir = os.path.join(config_dir,'N1_1080','run','dv')

    source_dict = {}

    for file_name in os.listdir(dv_dir):
        if 'beaufort_mask_'+test_var_name in file_name:
            file_path = os.path.join(config_dir, 'N1_1080', 'run', 'dv',file_name)
            grid = np.fromfile(file_path, '>f4')
            n_iters_in_this_file = int(np.size(grid) / (Nr * N))
            file_iter = int(file_name.split('.')[-2])

            first_iter = (file_iter-1)+int(iters_per_output/2)
            last_iter = (file_iter-1)+int(iters_per_output/2) + (n_iters_in_this_file-1)*iters_per_output

            first_date = datetime(1992,1,1) + timedelta(seconds=first_iter*timestep)
            last_date = datetime(1992, 1, 1) + timedelta(seconds=last_iter * timestep)

            source_dict[file_iter] = [first_iter,last_iter,first_date,last_date]

    output = 'iter_number\tfirst_iter\tlast_iter\tfirst_date\tlast_date'

    source_iters = sorted(list(source_dict.keys()))
    for iter_number in source_iters:
        output+='\n'+str(iter_number)+\
                '\t'+str(source_dict[iter_number][0])+\
                '\t'+str(source_dict[iter_number][1])+\
                '\t'+str(source_dict[iter_number][2])+\
                '\t'+str(source_dict[iter_number][3])

    output_file = os.path.join(config_dir,'N1_1080','run','dv','beaufort_mask_source.dict')
    f = open(output_file,'w')
    f.write(output)
    f.close()












if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="Path to the directory where config files are stored.", dest="config_dir",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    config_dir = args.config_dir

    create_beaufort_dataset(config_dir)