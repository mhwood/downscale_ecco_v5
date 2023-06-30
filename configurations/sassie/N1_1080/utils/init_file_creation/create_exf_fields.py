
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import xarray as xr
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
# sys.path.insert(1,'/Users/michwood/Documents/Research/Projects/Ocean_Modeling/Downscale_Baffin_Bay/MITgcm/configurations/baffin_bay_ECCOv4r4_forcing/utils')
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf

def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[mask_name]
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()
    return(faces,rows,cols)

def read_surface_variables_onto_N0_grid(N0_run_dir, var_name, nTimesteps, nTimestepsOut, faces, rows, cols, LLC, N):

    print('       - Reading N0 exf field ' + var_name +' output onto N0 domain ')

    # make a blank grid of zeros
    grid_faces = {}
    for i in [1,2,3,4,5]:
        if i < 3:
            grid_face = np.zeros((nTimestepsOut,LLC*3,LLC))
        if i == 3:
            grid_face = np.zeros((nTimestepsOut,LLC,LLC))
        if i > 3:
            grid_face = np.zeros((nTimestepsOut,LLC,LLC*3))
        grid_faces[i] = grid_face

    n_points = len(faces)

    var_file = os.path.join(N0_run_dir, 'mask_arctic_surface_' + var_name + '.0000000002.bin')
    var_grid = np.fromfile(var_file, dtype='>f4')
    var_grid = np.reshape(var_grid,(nTimesteps,n_points))
    var_grid = var_grid[:nTimestepsOut,:]

    for n in range(n_points-1):

        # if faces[n] in [1,2]:
        #     assign_row = rows[n]-2*LLC-(LLC-N)
        # else:
        #     assign_row = rows[n]

        grid_faces[faces[n]][:,rows[n],cols[n]] = var_grid[:,n]

    return(grid_faces)


def read_mask_faces_from_nc(nc_file_prefix, hFac='C',limit_extension=True):
    mask_faces = {}
    for face in range(1,6):
        ds = nc4.Dataset(nc_file_prefix+'_face_'+str(face)+'.nc')
        mask = ds.variables['wet_grid_'+hFac][:,:,:]
        if limit_extension:
            if face in [1,2]:
                mask = mask[:,:-1,1:-1]
            if face in [3]:
                mask = mask[:,1:-1,1:-1]
            if face in [4,5]:
                mask = mask[:,1:-1,1:]
        mask = mask[0,:,:]

        # plt.imshow(mask,origin='lower')
        # plt.title(nc_file_prefix+', '+str(face))
        # plt.show()

        mask_faces[face] = mask
        ds.close()
    return(mask_faces)


def downscale_N0_exf_field_to_N1(nTimestepsOut, N0_surface_var_faces_extended, regridding_weights,
                                 N0_wet_grid_2D_on_N1_faces, N1_wet_grid_2D_faces):

    # make a grid of zeros to fill in for each face
    N1_surface_var_faces = {}
    for face in range(1,6):
        N1_surface_var_faces[face] = np.zeros((nTimestepsOut,
                                               np.shape(N1_wet_grid_2D_faces[face])[0],
                                               np.shape(N1_wet_grid_2D_faces[face])[1]))

    for timestep in range(nTimestepsOut):

        print('          Working on timestep ' + str(timestep + 1) + ' of ' + str(nTimestepsOut))

        var_grid_faces = {}
        for face in range(1,6):
            var_grid_faces[face] = N0_surface_var_faces_extended[face][timestep,:,:]

        N1_var_faces = df.downscale_2D_field_faces(var_grid_faces, regridding_weights,
                                                          N0_wet_grid_2D_on_N1_faces, N1_wet_grid_2D_faces,
                                                          spread_horizontally=False)

        for face in range(1,6):
            N1_surface_var_faces[face][timestep,:,:] = N1_var_faces[face]

    return(N1_surface_var_faces)

def output_exf_variable_as_compact(output_dir, var_grid_faces, var_name,n,llc):

    var_grid_compact = sf.sassie_n1_faces_to_compact(var_grid_faces,llc,n,levels = np.shape(var_grid_faces[1])[0])

    print('   Creating exf field '+var_name+' of shape '+str(np.shape(var_grid_compact)))
    output_file = os.path.join(output_dir,'N1_exf_' + var_name + '.bin')
    var_grid_compact.ravel('C').astype('>f4').tofile(output_file)


########################################################################################################################

def create_exf_fields_via_interpolation(nTimesteps,nTimestepsOut,ecco_path):

    LLC = 270
    llc = 1080

    N = 180
    n = 720

    # this is the dir where the exf output will be stored
    input_dir = os.path.join('..', 'input')
    if 'exf' not in os.listdir(input_dir):
        os.mkdir(os.path.join(input_dir, 'exf'))
    output_dir = os.path.join(input_dir, 'exf')

    print('Creating the exf fields for the N1_'+str(llc)+' model from the output of the N0_'+str(LLC)+' model')

    print('    - Reading the mask to reference the variable to the LLC grid')
    nc_dict_file = os.path.join('..','..','N0_'+str(LLC),'input','dv_mask_reference_dict.nc')
    faces, rows, cols = read_mask_reference_from_nc_dict(nc_dict_file, 'mask_arctic_surface')

    N0_run_dir = os.path.join('..','..', 'N0_'+str(LLC), 'run', 'dv')
    surface_var_name_set = ['APRESS', 'AQH', 'ATEMP', 'LWDOWN', 'PRECIP', 'SWDOWN',
                            'USTRESS', 'VSTRESS', 'WSPEED', 'RUNOFF']
    # surface_var_name_set = ['WSPEED']

    print('    - Reading in the geometry of the N0_'+str(LLC)+' domain')
    N0_XC_faces, N0_YC_faces = sf.read_sassie_n0_grid(ecco_path)

    N0_XC_faces_extended = {}
    N0_YC_faces_extended = {}
    for face in range(1, 6):
        N0_XC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_XC_faces, face)
        N0_YC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_YC_faces, face)

    print('    - Reading in the geometry of the N1_' + str(LLC) + ' domain')
    N1_XC_faces, N1_YC_faces = sf.read_sassie_n1_grid(input_dir)

    print('    - Reading in the regridder functions')
    regridding_weights = df.read_regridder_functions_from_nc(input_dir, N0_XC_faces_extended, N0_YC_faces_extended, N1_XC_faces, N1_YC_faces)

    for ss in range(len(surface_var_name_set)):

        N0_surface_var_faces = read_surface_variables_onto_N0_grid(N0_run_dir, surface_var_name_set[ss], nTimesteps, nTimestepsOut,
                                                             faces,rows, cols, LLC, N)

        for face in range(1, 6):
            if face < 3:
                N0_surface_var_faces[face] = N0_surface_var_faces[face][:,-N:, :]
            if face > 3:
                N0_surface_var_faces[face] = N0_surface_var_faces[face][:,:, :N]

        N0_stitched_grid = sf.stitch_faces_to_single_grid(N0_surface_var_faces,270,180)
        # plt.imshow(stitched_grid[plot_timestep,:,:])
        # plt.show()

        N0_surface_var_faces_extended = {}
        for face in range(1, 6):
            N0_surface_var_faces_extended[face] = sf.get_extended_var_grid_on_face_3D(N0_surface_var_faces, face, np.shape(N0_surface_var_faces[face])[0])

        # stitched_grid = sf.stitch_faces_to_single_grid(N0_surface_var_faces_extended, 272, 181)
        # plt.imshow(stitched_grid[-1,:,:])
        # plt.show()

        if surface_var_name_set[ss] == 'VSTRESS':
            N1_wet_grid_2D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),hFac='S')
            N0_wet_grid_2D_on_N1_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='S')
        elif surface_var_name_set[ss] == 'USTRESS':
            N1_wet_grid_2D_faces = read_mask_faces_from_nc(os.path.join('..','input', 'N1_extended_wet_grids'), hFac='W')
            N0_wet_grid_2D_on_N1_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='W')
        else:
            N1_wet_grid_2D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'), hFac='C')
            N0_wet_grid_2D_on_N1_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='C')

        surface_var_name = surface_var_name_set[ss]
        print('        - Creating the '+surface_var_name+' conditions')

        N1_surface_var_faces = downscale_N0_exf_field_to_N1(nTimestepsOut, N0_surface_var_faces_extended, regridding_weights,
                                                            N0_wet_grid_2D_on_N1_faces, N1_wet_grid_2D_faces)

        # ###################################################
        # # Make a plot if desired
        # plot_timestep = 3
        #
        # N0_stitched_grid = sf.stitch_faces_to_single_grid(N0_surface_var_faces, LLC, N)
        # N1_stitched_grid = sf.stitch_faces_to_single_grid(N1_surface_var_faces, llc, n)
        #
        # plt.subplot(1,2,1)
        # C = plt.imshow(N0_stitched_grid[plot_timestep,:,:], origin='lower')
        # plt.colorbar(C)
        # plt.subplot(1, 2, 2)
        # C = plt.imshow(N1_stitched_grid[plot_timestep, :, :], origin='lower')
        # plt.colorbar(C)
        # plt.title(surface_var_name)
        # plt.show()

        output_exf_variable_as_compact(output_dir, N1_surface_var_faces, surface_var_name,n,llc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--nTimesteps", action="store",
                        help="The number of timesteps in the diagnostics_vec output from the N0_270 domain.",
                        dest="nTimesteps",
                        type=int, required=True)

    parser.add_argument("-o", "--nTimestepsOut", action="store",
                        help="The number of timesteps to output for the N1_1080 domain (counting by steps in the N0_270 domain).", dest="nTimestepsOut",
                        type=int, required=False, default=-1)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    nTimesteps = args.nTimesteps
    nTimestepsOut = args.nTimestepsOut
    ecco_path = args.ecco_path

    if nTimestepsOut==-1:
        nTimestepsOut = nTimesteps

    create_exf_fields_via_interpolation(nTimesteps,nTimestepsOut,ecco_path)