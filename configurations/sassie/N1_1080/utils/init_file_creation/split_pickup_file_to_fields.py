
import os
import numpy as np
from MITgcmutils import mds


def read_pickup_file_to_compact(input_init_dir,pickup_file='pickup.0000000001'):

    Nr = 90
    print('      Reading from '+input_init_dir + '/' + pickup_file)

    global_data, _, global_metadata = mds.rdmds(input_init_dir + '/' + pickup_file,
                                                returnmeta=True)
    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    all_var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        # var_grid_faces = llc_compact_to_faces(var_grid, less_output=True)

        all_var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,all_var_grids,global_metadata)



########################################################################################################################

def split_pickup_file():

    print('Splitting the pickup file for the N1_1080 model into components')

    print('    - Reading in variables from the N1_1080 pickup file')
    input_init_dir = os.path.join('..', 'input')
    var_names,row_bounds,all_var_grids,global_metadata = read_pickup_file_to_compact(input_init_dir)

    print('Writing Theta')
    theta_file = os.path.join(input_init_dir,'L1_1080_temp_init.bin')
    all_var_grids[var_names.index('Theta')].ravel(order='C').astype('>f4').tofile(theta_file)

    print('Writing Salt')
    theta_file = os.path.join(input_init_dir, 'L1_1080_salt_init.bin')
    all_var_grids[var_names.index('Salt')].ravel(order='C').astype('>f4').tofile(theta_file)

    print('Writing Uvel')
    theta_file = os.path.join(input_init_dir, 'L1_1080_uvel_init.bin')
    all_var_grids[var_names.index('Uvel')].ravel(order='C').astype('>f4').tofile(theta_file)

    print('Writing Theta')
    theta_file = os.path.join(input_init_dir, 'L1_1080_vvel_init.bin')
    all_var_grids[var_names.index('Vvel')].ravel(order='C').astype('>f4').tofile(theta_file)


if __name__ == '__main__':

    split_pickup_file()