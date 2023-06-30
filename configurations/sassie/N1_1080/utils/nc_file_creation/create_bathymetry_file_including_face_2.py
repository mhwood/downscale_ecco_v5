
import os
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

########################################################################################################################

def create_bathymetry_file(ecco_path,coarse_boundary):

    print('Creating the bathymetry grid for the N1_1080 configuration')

    config_dir = '../..'
    LLC = 270
    llc = 1080

    # read the llc and LLC bathy files to faces
    if coarse_boundary:
        print(' - Reading the LLC'+str(LLC)+' bathymetry onto the faces')
        LLC_bathy_file = os.path.join(ecco_path,'LLC'+str(LLC)+'_Files','input_init','bathy_llc'+str(LLC))
        LLC_bathy_grid = np.fromfile(LLC_bathy_file,'>f4')
        LLC_bathy_grid = np.reshape(LLC_bathy_grid,(13*LLC,LLC))
        LLC_bathy_faces = llc_compact_to_faces(LLC_bathy_grid, less_output=True)

    print(' - Reading the LLC' + str(llc) + ' bathymetry onto the faces')
    llc_bathy_file = os.path.join(ecco_path, 'LLC' + str(llc) + '_Files', 'input_init', 'bathy_llc' + str(llc))
    llc_bathy_grid = np.fromfile(llc_bathy_file, '>f4')
    llc_bathy_grid = np.reshape(llc_bathy_grid, (13 * llc, llc))
    llc_bathy_faces = llc_compact_to_faces(llc_bathy_grid, less_output=True)

    llc_truncated_faces = dict()
    llc_truncated_faces[1] = llc_bathy_faces[1][-680:, :] + 0 # note + 0 makes a copy
    llc_truncated_faces[2] = llc_bathy_faces[2][-680:, :] + 0
    llc_truncated_faces[3] = llc_bathy_faces[3] + 0
    llc_truncated_faces[4] = llc_bathy_faces[4][:, :680] + 0
    llc_truncated_faces[5] = llc_bathy_faces[5][:, :680] + 0

    # add truncated bathymetry on the boundary (if requested)
    if coarse_boundary:
        raise ValueError('Didnt finish the coarse boundary part of the code')

    # convert the faces to compact form
    llc_truncated_compact = sf.sassie_n1_faces_to_compact(llc_truncated_faces)

    # output the file
    output_file = os.path.join('..','input', 'bathymetry.bin')
    llc_truncated_compact.ravel('C').astype('>f4').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    parser.add_argument("-c", "--coarse_boundary", action="store",
                        help="Choose whether to use 'blocky' bathymetry on the boundaries.", dest="coarse_boundary",
                        type=bool, required=False, default=False)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    coarse_boundary = args.coarse_boundary

    create_bathymetry_file(ecco_path,coarse_boundary)