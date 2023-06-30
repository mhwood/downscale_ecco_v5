

import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import xarray as xr
import xesmf as xe
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf


########################################################################################################################




def get_N0_bathy_faces_on_N1_domain(N0_XC_faces, N0_YC_faces,N0_bathy_faces,N1_XC_faces, N1_YC_faces):

    N0_bathy_faces_on_N1 = {}

    for face in range(1,6):
        ds_in = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N0_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N0_YC_faces[face]).T),
            ))

        ds_out = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N1_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N1_YC_faces[face]).T),
            ))
        regridder = xe.Regridder(ds_in, ds_out, method='nearest_s2d')
        grid = regridder(np.asfortranarray(N0_bathy_faces[face]).T).T

        # plt.imshow(grid,origin='lower')
        # plt.title(str(face))
        # plt.show()

        N0_bathy_faces_on_N1[face]=grid

    return(N0_bathy_faces_on_N1)



def create_wetgrid_file(bathy_file,ecco_path):

    print('Outputting an nc file with the bathymetry and the wet grids')

    print(' - Reading from file '+bathy_file)

    n = 720
    llc = 1080

    input_dir = os.path.join('..', 'input')
    #     # bathy_file = 'bathymetry.bin'
    #     # bathy_compact = np.fromfile(os.path.join(input_dir, bathy_file), '>f4')
    #     # bathy_compact = np.reshape(bathy_compact, (4 * n + llc, llc))
    #     # bathy_faces = sf.sassie_n1_compact_to_faces(bathy_compact,n,llc)

    print(' - Interpolating the N0_270 bathymetry onto the new domain using nearest neighbors ')

    N0_bathy_faces = sf.read_sassie_n0_bathymetry(ecco_path)

    N0_XC_faces, N0_YC_faces = sf.read_sassie_n0_grid(ecco_path)

    N1_XC_faces, N1_YC_faces = sf.read_sassie_n1_grid(input_dir)

    N0_bathy_faces_on_N1 = get_N0_bathy_faces_on_N1_domain(N0_XC_faces, N0_YC_faces,N0_bathy_faces,N1_XC_faces, N1_YC_faces)

    delR = sf.read_sassie_delR(domain_level=1)

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

        extended_bathy_face = sf.get_extended_var_grid_on_face(N0_bathy_faces_on_N1, face)

        print('    - Calculating the hFacC wetgrid')
        wet_grid_C = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='C', hFacMin=hFacMin, hFacMinDr=hFacMinDr)
        # hFacC_faces[face] = wet_grid_C

        print('    - Calculating the hFacS wetgrid')
        wet_grid_S = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='S', hFacMin=hFacMin, hFacMinDr=hFacMinDr)
        # hFacS_faces[face] = wet_grid_S

        print('    - Calculating the hFacW wetgrid')
        wet_grid_W = df.create_3D_wet_grid(extended_bathy_face, delR, hFac='W', hFacMin=hFacMin, hFacMinDr=hFacMinDr)
        # hFacW_faces[face] = wet_grid_W

        output_file = os.path.join('..',  'input', 'N0_extended_wet_grids_on_N1_face_'+str(face)+'.nc')

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

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    parser.add_argument("-b", "--bathy_file", action="store",
                        help="The bathymetry file to read in.", dest="bathy_file",
                        type=str, required=False, default='bathymetry.bin')

    args = parser.parse_args()
    bathy_file = args.bathy_file
    ecco_path = args.ecco_path

    create_wetgrid_file(bathy_file,ecco_path)