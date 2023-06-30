

import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from pyproj import Transformer
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf


########################################################################################################################

def create_wetgrid_file(bathy_file):

    print('Outputting an nc file with the bathymetry and the wet grids')

    print(' - Reading from file '+bathy_file)

    n = 720
    llc = 1080

    input_dir = os.path.join('..', 'input')
    bathy_file = 'bathymetry.bin'
    bathy_compact = np.fromfile(os.path.join(input_dir, bathy_file), '>f4')
    bathy_compact = np.reshape(bathy_compact, (4 * n + llc, llc))
    bathy_faces = sf.sassie_n1_compact_to_faces(bathy_compact,n,llc)

    # delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    #                  10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    #                  31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    #                  93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    #                  139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    #                  341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    delR = np.array([1.00,    1.14,    1.30,    1.49,   1.70,
                      1.93,    2.20,    2.50,    2.84,   3.21,
                      3.63,    4.10,    4.61,    5.18,   5.79,
                      6.47,    7.20,    7.98,    8.83,   9.73,
                     10.69,   11.70,   12.76,   13.87,  15.03,
                     16.22,   17.45,   18.70,   19.97,  21.27,
                     22.56,   23.87,   25.17,   26.46,  27.74,
                     29.00,   30.24,   31.45,   32.65,  33.82,
                     34.97,   36.09,   37.20,   38.29,  39.37,
                     40.45,   41.53,   42.62,   43.73,  44.87,
                     46.05,   47.28,   48.56,   49.93,  51.38,
                     52.93,   54.61,   56.42,   58.38,  60.53,
                     62.87,   65.43,   68.24,   71.33,  74.73,
                     78.47,   82.61,   87.17,   92.21,  97.79,
                    103.96,  110.79,  118.35,  126.73, 136.01,
                    146.30,  157.71,  170.35,  184.37, 199.89,
                    217.09,  236.13,  257.21,  280.50, 306.24,
                    334.64,  365.93,  400.38,  438.23, 479.74])

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]),Z_bottom[:-1]])
    Z = (Z_bottom + Z_top)/2

    # hFacC_faces = {}
    # hFacS_faces = {}
    # hFacW_faces = {}

    faces = [1,2,3,4,5]

    hFacMin = 0.3
    hFacMinDr = 1.0

    for face in faces:
        print(' - Calculating fields for face '+str(face))

        extended_bathy_face = sf.get_extended_var_grid_on_face(bathy_faces, face)

        print('    - Calculating the hFacC wetgrid')
        wet_grid_C = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='C', hFacMin=hFacMin, hFacMinDr=hFacMinDr)
        # hFacC_faces[face] = wet_grid_C

        print('    - Calculating the hFacS wetgrid')
        wet_grid_S = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='S', hFacMin=hFacMin, hFacMinDr=hFacMinDr)
        # hFacS_faces[face] = wet_grid_S

        print('    - Calculating the hFacW wetgrid')
        wet_grid_W = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='W', hFacMin=hFacMin, hFacMinDr=hFacMinDr)
        # hFacW_faces[face] = wet_grid_W


        output_file = os.path.join('..',  'input', 'N1_extended_wet_grids_face_'+str(face)+'.nc')

        ds = nc4.Dataset(output_file,'w')

        ds.createDimension('levels', len(delR))
        ds.createDimension('rows',np.shape(wet_grid_C)[1])
        ds.createDimension('cols', np.shape(wet_grid_C)[2])

        var = ds.createVariable('wet_grid_C', 'f4', ('levels','rows', 'cols'))
        var[:,:,:] = wet_grid_C
        var = ds.createVariable('wet_grid_W', 'f4', ('levels','rows', 'cols'))
        var[:,:,:] = wet_grid_W
        var = ds.createVariable('wet_grid_S', 'f4', ('levels','rows', 'cols'))
        var[:,:,:] = wet_grid_S

        ds.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bathy_file", action="store",
                        help="The bathymetry file to read in.", dest="bathy_file",
                        type=str, required=False, default='bathymetry.bin')

    args = parser.parse_args()
    bathy_file = args.bathy_file

    create_wetgrid_file(bathy_file)