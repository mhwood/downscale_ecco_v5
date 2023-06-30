

import os
import simplegrid as sg
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import xarray as xr
import xesmf as xe
import argparse
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf

########################################################################################################################


def create_regridding_files():
    config_dir = os.path.join('..','..')

    N0_XC_faces, N0_YC_faces = sf.read_sassie_n0_grid(ecco_path)

    N0_XC_faces_extended = {}
    N0_YC_faces_extended = {}
    for face in range(1,6):
        N0_XC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_XC_faces, face)
        N0_YC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_YC_faces, face)

    N1_XC_faces, N1_YC_faces = sf.read_sassie_n1_grid(input_dir=os.path.join('..', 'input'))

    # N1_XC_faces_extended = {}
    # N1_YC_faces_extended = {}
    # for face in range(1, 6):
    #     N1_XC_faces_extended[face] = sf.get_extended_var_grid_on_face(N1_XC_faces, face)
    #     N1_YC_faces_extended[face] = sf.get_extended_var_grid_on_face(N1_YC_faces, face)

    for face in range(1,6):
        print('    Creating the regridder for face '+str(face))
        print('        Face '+str(face)+' shape:')
        print('           L0: ',np.shape(N0_XC_faces_extended[face]))
        print('           L1: ', np.shape(N1_XC_faces[face]))

        filename = os.path.join(config_dir,'N1_1080','input','regridding_weights_face_'+str(face)+'.nc')

        ds_in = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N0_XC_faces_extended[face]).T),
                lat=(["x", "y"], np.asfortranarray(N0_YC_faces_extended[face]).T),
            ))

        ds_out = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N1_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N1_YC_faces[face]).T),
            ))
        regridder = xe.Regridder(ds_in, ds_out, method='bilinear', filename=filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_regridding_files()