
import os
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

########################################################################################################################

def create_bathymetry_file(ecco_path):

    print('Creating the bathymetry grid for the N1_1080 configuration')

    config_dir = '../..'

    llc = 2160
    n = 320

    print(' - Reading the LLC' + str(llc) + ' bathymetry onto the faces')
    llc_bathy_file = os.path.join(ecco_path, 'LLC' + str(llc) + '_Files', 'input_init', 'bathy_llc' + str(llc))
    llc_bathy_grid = np.fromfile(llc_bathy_file, '>f4')
    llc_bathy_grid = np.reshape(llc_bathy_grid, (13 * llc, llc))
    llc_bathy_faces = llc_compact_to_faces(llc_bathy_grid, less_output=True)

    total_rows = n
    zerod_rows = 0

    llc_truncated_faces = dict()
    llc_truncated_faces[1] = llc_bathy_faces[1][-total_rows:, :] + 0 # note + 0 makes a copy
    llc_truncated_faces[2] = llc_bathy_faces[2][-total_rows:, :] + 0
    llc_truncated_faces[3] = llc_bathy_faces[3] + 0
    llc_truncated_faces[4] = llc_bathy_faces[4][:, :total_rows] + 0
    llc_truncated_faces[5] = llc_bathy_faces[5][:, :total_rows] + 0

    if zerod_rows>0:
        llc_truncated_faces[1][:zerod_rows, :] = 0
        llc_truncated_faces[2][:zerod_rows, :] = 0
        llc_truncated_faces[4][:, -zerod_rows:] = 0
        llc_truncated_faces[5][:, -zerod_rows:] = 0

    # # fill in the med and the caspian, etc
    # llc_truncated_faces[1][:90-zero_row_shift, 460:] = 0
    # llc_truncated_faces[1][:2*90-zero_row_shift, 500:] = 0
    # llc_truncated_faces[1][:270-zero_row_shift, 750:] = 0
    # llc_truncated_faces[2][:270-zero_row_shift, :500] = 0
    #
    # # fill in the baltic
    # llc_truncated_faces[1][300-zero_row_shift:480-zero_row_shift, 600:840] = 0
    # llc_truncated_faces[1][480-zero_row_shift:600-zero_row_shift, 660:800] = 0

    # convert the faces to compact form
    llc_truncated_compact = sf.sassie_n1_faces_to_compact(llc_truncated_faces,llc,n)

    # output the file
    output_file = os.path.join('..','input', 'bathy_N2_no_walls')
    llc_truncated_compact.ravel('C').astype('>f4').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_bathymetry_file(ecco_path)