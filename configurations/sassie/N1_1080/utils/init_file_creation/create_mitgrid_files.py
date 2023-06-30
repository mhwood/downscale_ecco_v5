
import os
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

########################################################################################################################

def create_mitgrid_files(ecco_path):

    print('Creating the mitgrid files for the N1_1080 configuration')

    config_dir = '../..'
    llc = 1080
    n = 680

    for face in range(1,6):

        print(' - Reading the LLC' + str(llc) + ' mitgrid for face '+str(face))
        llc_mitgrid_file = os.path.join(ecco_path, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles', 'tile00' + str(face)+'.mitgrid')
        llc_mitgrid_grid = np.fromfile(llc_mitgrid_file, '>f8')
        if face in [1,2]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16,3*llc+1, llc+1))
            llc_mitgrid_grid = llc_mitgrid_grid[:, -n-1:, :]
        if face in [3]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16,llc+1, llc+1))
        if face in [4,5]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16,llc+1, 3*llc+1))
            llc_mitgrid_grid = llc_mitgrid_grid[:, :, :n+1]

        # output the file
        output_file = os.path.join('..','input', 'tile00' + str(face)+'.mitgrid')
        llc_mitgrid_grid.ravel('C').astype('>f8').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_mitgrid_files(ecco_path)