

import os
import simplegrid as sg
import numpy as np
import netCDF4 as nc4
import xarray as xr
import xesmf as xe
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df


########################################################################################################################

def read_mask_reference_from_nc_dict(nc_dict_file, mask_name):
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[mask_name]
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()
    return (faces, rows, cols)

def read_L0_faces_geometry(ecco_dir, llc):
    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in [1, 2, 3, 4, 5]:
        if i < 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),
                                                   llc,
                                                   3 * llc)
        if i == 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),
                                                   llc,
                                                   llc)
        if i > 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),
                                                   3 * llc, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    # bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init','bathy_llc'+str(llc))
    # bathy_compact = np.fromfile(bathy_file,'>f4')
    # bathy_compact = np.reshape(bathy_compact,(13*llc,llc))
    # bathy_faces = llc_compact_to_faces(bathy_compact,less_output=True)

    return (XC_faces, YC_faces)

def read_surface_variables_onto_L0_grid(L0_run_dir, var_name, file_name, nTimesteps, faces, rows, cols, llc):
    print('       - Reading L0 exf field ' + var_name + ' output onto N0 domain ')

    # make a blank grid of zeros
    grid_faces = {}
    for i in [1, 2, 3, 4, 5]:
        if i < 3:
            grid_face = np.zeros((llc * 3, llc))
        if i == 3:
            grid_face = np.zeros((llc, llc))
        if i > 3:
            grid_face = np.zeros((llc, llc * 3))
        grid_faces[i] = grid_face

    var_file = os.path.join(L0_run_dir, 'mask_arctic_surface', var_name, file_name)
    var_grid = np.fromfile(var_file, dtype='>f4')
    N = int(np.size(var_grid) / nTimesteps)
    var_grid = np.reshape(var_grid, (nTimesteps, N))
    var_grid = var_grid[0, :]

    for n in range(N):
        grid_faces[faces[n]][rows[n], cols[n]] = var_grid[n]

    return (grid_faces)

def read_L0_variables(config_dir):
    LLC = 270

    nc_dict_file = os.path.join(config_dir, 'N0_' + str(LLC), 'input', 'dv_mask_reference_dict.nc')
    faces, rows, cols = read_mask_reference_from_nc_dict(nc_dict_file, 'mask_arctic_surface')

    ecco_path = '/Users/michwood/Documents/Research/Projects/Ocean_Modeling/ECCO'
    L0_XC_faces, L0_YC_faces = read_L0_faces_geometry(ecco_path, LLC)

    var_name = 'WSPEED'
    L0_run_dir = config_dir + '/N0_270/run/dv'
    file_name = 'mask_arctic_surface_' + var_name + '.0000000002.bin'
    nTimesteps = 72
    llc = 270
    L0_var_faces = read_surface_variables_onto_L0_grid(L0_run_dir, var_name, file_name, nTimesteps, faces, rows,
                                                       cols, llc)

    N = 680
    for face in range(1, 6):
        if face < 3:
            L0_XC_faces[face] = L0_XC_faces[face][-int(N / 4):, :]
            L0_YC_faces[face] = L0_YC_faces[face][-int(N / 4):, :]
            L0_var_faces[face] = L0_var_faces[face][-int(N / 4):, :]
        if face > 3:
            L0_XC_faces[face] = L0_XC_faces[face][:, :int(N / 4):]
            L0_YC_faces[face] = L0_YC_faces[face][:, :int(N / 4):]
            L0_var_faces[face] = L0_var_faces[face][:, :int(N / 4):]

    return (L0_XC_faces, L0_YC_faces, L0_var_faces)

def read_L1_faces_geometry(config_dir, llc, rows, cols):
        N = int((rows - llc) / 4)

        grid_file_dir = os.path.join(config_dir, 'N1_1080', 'input')
        XC_faces = {}
        YC_faces = {}
        for i in [1, 2, 3, 4, 5]:
            if i < 3:
                grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),
                                                       llc,
                                                       N)
            if i == 3:
                grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),
                                                       llc,
                                                       llc)
            if i > 3:
                grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'),
                                                       N,
                                                       llc)
            XC_face = grid_dict['XC'].T
            YC_face = grid_dict['YC'].T
            XC_faces[i] = XC_face
            YC_faces[i] = YC_face

        return (XC_faces, YC_faces, N)

def read_L1_variables(config_dir):


    llc = 1080
    n_rows_L1 = 3800
    n_cols_L1 = 1080

    print('    - Reading in the geometry of the L1_' + str(llc) + ' domain')
    L1_XC_faces, L1_YC_faces, N = read_L1_faces_geometry(config_dir, llc, n_rows_L1, n_cols_L1)

    return (L1_XC_faces, L1_YC_faces)

def create_regridding_files():
    config_dir = os.path.join('..','..')

    # read in example grids (face 1 and face 4)
    L0_XC_faces, L0_YC_faces, L0_var_faces = read_L0_variables(config_dir)

    L1_XC_faces, L1_YC_faces = read_L1_variables(config_dir)

    for face in range(1,6):
        print('    Creating the regridder for face '+str(face))
        # print('        Face '+str(face)+' shape:')
        # print('           L0: ',np.shape(L0_XC_faces[face]),np.shape(L0_YC_faces[face]))
        # print('           L1: ', np.shape(L1_XC_faces[face]),np.shape(L0_YC_faces[face]))

        filename = os.path.join(config_dir,'N1_1080','input','face_'+str(face)+'_downscale_weights.nc')

        ds_in = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(L0_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(L0_YC_faces[face]).T),
            ))

        ds_out = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(L1_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(L1_YC_faces[face]).T),
            ))
        regridder = xe.Regridder(ds_in, ds_out, method='bilinear', filename=filename)

if __name__ == '__main__':

    create_regridding_files()