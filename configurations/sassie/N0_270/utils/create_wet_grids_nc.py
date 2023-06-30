

import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
from pyproj import Transformer
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf


########################################################################################################################

def get_extended_bathy_grid(bathy_faces, face):

    bathy_face = bathy_faces[face]

    if face == 1 or face == 2:
        extended_bathy_face = np.zeros((np.shape(bathy_face)[0]+1,np.shape(bathy_face)[1]+2))
        extended_bathy_face[:-1, 1:-1] = bathy_face
        if face == 1:
            top = np.rot90(bathy_faces[3],k=3)
            left = np.rot90(bathy_faces[5])
            right = bathy_faces[2]
        if face == 2:
            top = bathy_faces[3]
            left = bathy_faces[1]
            right = np.rot90(bathy_faces[4])

        extended_bathy_face[:-1,0] = left[:,-1]
        extended_bathy_face[:-1, -1] = right[:, -1]
        extended_bathy_face[-1,1:-1] = top[0,:]

        # plt.subplot(2,3,2)
        # plt.imshow(top,origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(2, 3, 4)
        # plt.imshow(left, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(2, 3, 5)
        # plt.imshow(bathy_face, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(2, 3, 6)
        # plt.imshow(right, origin='lower',vmin=-10,vmax = 0)
        # plt.show()

    if face == 3:
        extended_bathy_face = np.zeros((np.shape(bathy_face)[0]+2,np.shape(bathy_face)[1]+2))
        extended_bathy_face[1:-1, 1:-1] = bathy_face

        top = np.rot90(bathy_faces[5],k=3)
        right = bathy_faces[4]
        left = np.rot90(bathy_faces[1])
        bottom = bathy_faces[2]

        extended_bathy_face[1:-1,0] = left[:,-1]
        extended_bathy_face[0, 1:-1] = bottom[-1, :]
        extended_bathy_face[-1,1:-1] = top[0,:]
        extended_bathy_face[1:-1, -1] = right[:, 0]

        # plt.subplot(3,3,2)
        # plt.imshow(top,origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 4)
        # plt.imshow(left, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 5)
        # plt.imshow(bathy_face, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 8)
        # plt.imshow(bottom, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 6)
        # plt.imshow(right, origin='lower', vmin=-10, vmax=0)
        # plt.show()

    if face == 4 or face == 5:
        extended_bathy_face = np.zeros((np.shape(bathy_face)[0]+2,np.shape(bathy_face)[1]+1))
        extended_bathy_face[1:-1, 1:] = bathy_face
        if face == 4:
            top = bathy_faces[5]
            left = bathy_faces[3]
            bottom = np.rot90(bathy_faces[2],k=3)
        if face == 5:
            top = np.rot90(bathy_faces[1],k=3)
            left = np.rot90(bathy_faces[3])
            bottom = bathy_faces[4]

        extended_bathy_face[1:-1,0] = left[:,-1]
        extended_bathy_face[0, 1:] = bottom[-1, :]
        extended_bathy_face[-1,1:] = top[0,:]

        # plt.subplot(3,2,2)
        # plt.imshow(top,origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 2, 3)
        # plt.imshow(left, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 2, 4)
        # plt.imshow(bathy_face, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 2, 6)
        # plt.imshow(bottom, origin='lower',vmin=-10,vmax = 0)
        # plt.show()
    return(extended_bathy_face)


def create_wetgrids_file(ecco_path,rows):

    print('Creating an nc file with the wet grids for the N0_270 domain')

    LLC = 270
    print(' - Reading the LLC' + str(LLC) + ' bathymetry onto the faces')
    LLC_bathy_file = os.path.join(ecco_path, 'LLC' + str(LLC) + '_Files', 'input_init', 'bathy_LLC' + str(LLC))
    LLC_bathy_grid = np.fromfile(LLC_bathy_file, '>f4')
    LLC_bathy_grid = np.reshape(LLC_bathy_grid, (13 * LLC, LLC))
    LLC_bathy_faces = llc_compact_to_faces(LLC_bathy_grid, less_output=True)

    delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                     10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                     31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                     93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                     139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                     341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    LLC_truncated_faces = dict()
    LLC_truncated_faces[1] = LLC_bathy_faces[1][-rows:, :] + 0  # note + 0 makes a copy
    LLC_truncated_faces[2] = LLC_bathy_faces[2][-rows:, :] + 0
    LLC_truncated_faces[3] = LLC_bathy_faces[3] + 0
    LLC_truncated_faces[4] = LLC_bathy_faces[4][:, :rows] + 0
    LLC_truncated_faces[5] = LLC_bathy_faces[5][:, :rows] + 0

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]),Z_bottom[:-1]])
    Z = (Z_bottom + Z_top)/2

    wet_grid_C_faces = {}
    wet_grid_S_faces = {}
    wet_grid_W_faces = {}

    for face in range(1,6):
        print('   - Calculating fields for face '+str(face))

        extended_bathy_face = get_extended_bathy_grid(LLC_truncated_faces, face)

        print('   - Calculating the hFacC wetgrid')
        wet_grid_C = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='C')
        wet_grid_C_faces[face] = wet_grid_C

        print('   - Calculating the hFacS wetgrid')
        wet_grid_S = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='S')
        wet_grid_S_faces[face] = wet_grid_S

        print('   - Calculating the hFacW wetgrid')
        wet_grid_W = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='W')
        wet_grid_W_faces[face] = wet_grid_W

    print('Get the script from the existing function that chops back the extended grid')

    # if bathy_file[:-4]+'.nc' in os.listdir(os.path.join('..', 'input')):
    #     os.remove(os.path.join('..',  'input', bathy_file[:-4]+'.nc'))
    #
    # output_file = os.path.join('..',  'input', bathy_file[:-4]+'.nc')
    # ds = nc4.Dataset(output_file,'w')
    #
    # ds.createDimension('x',np.shape(XC)[1])
    # ds.createDimension('y', np.shape(XC)[0])
    # ds.createDimension('z', len(delR))
    #
    # ds.source=bathy_file
    #
    # var = ds.createVariable('XC','f4',('y','x'))
    # var[:] = XC
    # var = ds.createVariable('YC', 'f4', ('y','x'))
    # var[:] = YC
    # var = ds.createVariable('z', 'f4', ('z',))
    # var[:] = Z
    #
    # var = ds.createVariable('bathymetry','f4',('y','x'))
    # var[:,:] = bathy_grid
    # var = ds.createVariable('wet_grid_C', 'f4', ('z','y', 'x'))
    # var[:,:,:] = wet_grid_C
    # var = ds.createVariable('wet_grid_W', 'f4', ('z','y', 'x'))
    # var[:,:,:] = wet_grid_W
    # var = ds.createVariable('wet_grid_S', 'f4', ('z','y', 'x'))
    # var[:,:,:] = wet_grid_S
    #
    # ds.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    parser.add_argument("-r", "--rows", action="store",
                        help="The number of rows that extend down from face 3.", dest="rows",
                        type=int, required=True, default=180)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    rows = args.rows

    create_wetgrids_file(ecco_path,rows)