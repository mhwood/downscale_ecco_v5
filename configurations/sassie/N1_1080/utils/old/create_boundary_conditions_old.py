
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata, interp1d
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
# sys.path.insert(1,'/Users/michwood/Documents/Research/Projects/Ocean_Modeling/Downscale_Baffin_Bay/MITgcm/configurations/baffin_bay_ECCOv4r4_forcing/utils')
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf

def read_N0_faces_geometry(ecco_dir,llc):
    grid_file_dir = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in [1,2,3,4,5]:
        if i<3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, 3*llc)
        if i==3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, llc)
        if i>3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3*llc, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    # bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init','bathy_llc'+str(llc))
    # bathy_compact = np.fromfile(bathy_file,'>f4')
    # bathy_compact = np.reshape(bathy_compact,(13*llc,llc))
    # bathy_faces = llc_compact_to_faces(bathy_compact,less_output=True)

    return(XC_faces,YC_faces)

def read_N1_faces_geometry(llc,rows,cols):
    N = int((rows - llc)/4)

    grid_file_dir = os.path.join('..','input')
    XC_faces = {}
    YC_faces = {}
    for i in [1,2,3,4,5]:
        if i<3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, N)
        if i==3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, llc)
        if i>3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), N, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    return(XC_faces,YC_faces,N)

def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[mask_name]
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()
    return(faces,rows,cols)

def read_boundary_variable_onto_N0_grid(N0_run_dir, mask_name, var_name,year,
                                        nTimesteps, nTimestepsOut, faces, rows, cols, LLC, N, Nr,
                                        N0_boundary_var_faces):

    n_points = len(faces)

    if year==1992:
        file_time = 2
    else:
        file_time = int(1 + (year-1992)*24*60*60*365/1200)

    var_file = os.path.join(N0_run_dir,mask_name,var_name, mask_name + '_' + var_name + '.'+'{:010d}'.format(file_time)+'.bin')
    var_grid = np.fromfile(var_file, dtype='>f4')
    var_grid = np.reshape(var_grid, (nTimesteps, Nr, n_points))
    var_grid = var_grid[:nTimestepsOut,:,:]

    for n in range(n_points):

        if faces[n] in [1,2]:
            assign_row = rows[n]-2*LLC-(LLC-N)
        else:
            assign_row = rows[n]

        # if n%100==0:
        #     print(n,N)
        if faces[n]!=0:
            N0_boundary_var_faces[faces[n]][:,:,assign_row,cols[n]] = var_grid[:,:,n]

    return(N0_boundary_var_faces)

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

        # plt.imshow(mask,origin='lower')
        # plt.title(nc_file_prefix+', '+str(face))
        # plt.show()

        mask_faces[face] = mask
        ds.close()
    return(mask_faces)


def downscale_N0_bc_field_to_N1_old(nTimestepsOut, var_name,
                              N0_XC_faces,N0_YC_faces,N0_boundary_var_faces,
                              N0_wet_grid_3D, N0_wet_grid_on_N1_3D,
                              N1_XC_faces, N1_YC_faces, N1_wet_grid_3D):

    N1_boundary_var = np.zeros((nTimestepsOut,
                               np.shape(N1_wet_grid_3D[1])[0],
                               np.shape(N1_XC_faces[1])[1]+np.shape(N1_XC_faces[2])[1]+
                               np.shape(N1_XC_faces[4])[0]+np.shape(N1_XC_faces[5])[0]))

    #############################################################################################################
    # first, run the downscaling on faces 1 + 2

    # the data is only in the bottom 4 rows for the N0
    N0_XC_subset = np.concatenate([N0_XC_faces[1][:4, :],N0_XC_faces[2][:4, :]],axis=1)
    N0_YC_subset = np.concatenate([N0_YC_faces[1][:4, :],N0_YC_faces[2][:4, :]],axis=1)
    N0_boundary_var_subset = np.concatenate([N0_boundary_var_faces[1][:, :, :4, :],
                                             N0_boundary_var_faces[2][:, :, :4, :]],axis=3)

    # the output will be on the bottom row of the N1
    N1_XC_subset = np.concatenate([N1_XC_faces[1][2, :],N1_XC_faces[2][2, :]],axis=0)
    N1_YC_subset = np.concatenate([N1_YC_faces[1][2, :],N1_YC_faces[2][2, :]],axis=0)

    N1_XC_subset = np.reshape(N1_XC_subset, (1,np.size(N1_XC_subset)))
    N1_YC_subset = np.reshape(N1_YC_subset, (1,np.size(N1_YC_subset)))

    N1_wet_grid_subset = np.concatenate([N1_wet_grid_3D[1][:, 2, :],
                                         N1_wet_grid_3D[2][:, 2, :]],axis=1)
    N1_wet_grid_subset = np.reshape(N1_wet_grid_subset,
                                    (np.shape(N1_wet_grid_subset)[0],1,np.shape(N1_wet_grid_subset)[1]))

    N0_wet_grid_subset = []
    N0_wet_grid_on_N1_subset = []

    # print(np.shape(N0_XC_subset),np.shape(N0_YC_subset))
    # print(np.shape(N0_boundary_var_subset))
    # print(np.shape(N1_XC_subset),np.shape(N1_YC_subset))
    # print(np.shape(N1_wet_grid_subset))

    # plt.plot(N0_XC_subset,N0_YC_subset,'k.')
    # plt.plot(N1_XC_subset, N1_YC_subset, 'b.')
    # plt.show()

    # plt.imshow(N0_boundary_var_subset[0,:,1,:])
    # plt.show()

    # plt.imshow(N1_wet_grid_subset[:, 0, :])
    # plt.show()

    for timestep in range(nTimestepsOut):
        #if timestep%50==0:
        print('          Working on faces 1+2 in timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name)

        downscaled_field = df.downscale_3D_field(N0_XC_subset, N0_YC_subset,
                                                 N0_boundary_var_subset[timestep, :, :, :],
                                                 N0_wet_grid_subset, N0_wet_grid_on_N1_subset,
                                                 N1_XC_subset, N1_YC_subset, N1_wet_grid_subset)

        N1_boundary_var[timestep, :, :np.shape(N1_XC_faces[1])[1]+np.shape(N1_XC_faces[2])[1]] = downscaled_field[:,0,:]

        # C = plt.imshow(downscaled_field[:,0,:],aspect='equal')
        # plt.colorbar(C)
        # plt.title(timestep)
        # plt.show()

    #############################################################################################################
    # next, run the downscaling on faces 4 + 5

    # the data is only in the bottom 4 rows for the N0
    N0_XC_subset = np.concatenate([N0_XC_faces[4][:, -4:],N0_XC_faces[5][:, -4:]],axis=0)
    N0_YC_subset = np.concatenate([N0_YC_faces[4][:, -4:],N0_YC_faces[5][:, -4:]],axis=0)
    N0_boundary_var_subset = np.concatenate([N0_boundary_var_faces[4][:, :, :, -4:],
                                             N0_boundary_var_faces[5][:, :, :, -4:]],axis=2)

    N0_XC_subset[N0_XC_subset<0]+=360

    # the output will be on the bottom row of the N1
    N1_XC_subset = np.concatenate([N1_XC_faces[4][:, -3],N1_XC_faces[5][:, -3]],axis=0)
    N1_YC_subset = np.concatenate([N1_YC_faces[4][:, -3],N1_YC_faces[5][:, -3]],axis=0)

    N1_XC_subset = np.reshape(N1_XC_subset, (np.size(N1_XC_subset),1))
    N1_YC_subset = np.reshape(N1_YC_subset, (np.size(N1_YC_subset),1))

    N1_XC_subset[N1_XC_subset < 0] += 360

    N1_wet_grid_subset = np.concatenate([N1_wet_grid_3D[4][:, :, -3],
                                         N1_wet_grid_3D[5][:, :, -3]],axis=1)
    N1_wet_grid_subset = np.reshape(N1_wet_grid_subset,
                                    (np.shape(N1_wet_grid_subset)[0],np.shape(N1_wet_grid_subset)[1],1))

    N0_wet_grid_subset = []
    N0_wet_grid_on_N1_subset = []

    # print(np.shape(N0_XC_subset),np.shape(N0_YC_subset))
    # print(np.shape(N0_boundary_var_subset))
    # print(np.shape(N1_XC_subset),np.shape(N1_YC_subset))
    # print(np.shape(N1_wet_grid_subset))

    # plt.plot(N0_XC_subset,N0_YC_subset,'k.')
    # plt.plot(N1_XC_subset, N1_YC_subset, 'b.')
    # plt.show()
    #
    # plt.imshow(N0_boundary_var_subset[0,:,:,0])
    # plt.show()
    #
    # plt.imshow(N1_wet_grid_subset[:, :, 0])
    # plt.show()

    for timestep in range(nTimestepsOut):
        #if timestep%50==0:
        print('          Working on faces 4+5 in timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name)

        downscaled_field = df.downscale_3D_field(N0_XC_subset, N0_YC_subset,
                                                 N0_boundary_var_subset[timestep, :, :, :],
                                                 N0_wet_grid_subset, N0_wet_grid_on_N1_subset,
                                                 N1_XC_subset, N1_YC_subset, N1_wet_grid_subset)

        N1_boundary_var[timestep, :, np.shape(N1_XC_faces[1])[1] + np.shape(N1_XC_faces[2])[1]:] = downscaled_field[:,:,0]

        # C = plt.imshow(downscaled_field[:,:,0],aspect='equal')
        # plt.colorbar(C)
        # plt.title(timestep)
        # plt.show()

    return(N1_boundary_var)

def downscale_N0_bc_field_to_N1(nTimestepsOut, var_name,
                                var_grid_faces_extended, regridding_weights,
                                N0_wet_grid_3D_on_N1_faces, N1_wet_grid_3D_faces,
                                spread_horizontally=False):

    N1_boundary_var = np.zeros((nTimestepsOut,
                               np.shape(N1_wet_grid_3D[1])[0],
                               np.shape(N1_XC_faces[1])[1]+np.shape(N1_XC_faces[2])[1]+
                               np.shape(N1_XC_faces[4])[0]+np.shape(N1_XC_faces[5])[0]))

    #############################################################################################################
    # first, run the downscaling on faces 1 + 2

    # the data is only in the bottom 4 rows for the N0
    N0_XC_subset = np.concatenate([N0_XC_faces[1][:4, :],N0_XC_faces[2][:4, :]],axis=1)
    N0_YC_subset = np.concatenate([N0_YC_faces[1][:4, :],N0_YC_faces[2][:4, :]],axis=1)
    N0_boundary_var_subset = np.concatenate([N0_boundary_var_faces[1][:, :, :4, :],
                                             N0_boundary_var_faces[2][:, :, :4, :]],axis=3)

    # the output will be on the bottom row of the N1
    N1_XC_subset = np.concatenate([N1_XC_faces[1][2, :],N1_XC_faces[2][2, :]],axis=0)
    N1_YC_subset = np.concatenate([N1_YC_faces[1][2, :],N1_YC_faces[2][2, :]],axis=0)

    N1_XC_subset = np.reshape(N1_XC_subset, (1,np.size(N1_XC_subset)))
    N1_YC_subset = np.reshape(N1_YC_subset, (1,np.size(N1_YC_subset)))

    N1_wet_grid_subset = np.concatenate([N1_wet_grid_3D[1][:, 2, :],
                                         N1_wet_grid_3D[2][:, 2, :]],axis=1)
    N1_wet_grid_subset = np.reshape(N1_wet_grid_subset,
                                    (np.shape(N1_wet_grid_subset)[0],1,np.shape(N1_wet_grid_subset)[1]))

    N0_wet_grid_subset = []
    N0_wet_grid_on_N1_subset = []

    # print(np.shape(N0_XC_subset),np.shape(N0_YC_subset))
    # print(np.shape(N0_boundary_var_subset))
    # print(np.shape(N1_XC_subset),np.shape(N1_YC_subset))
    # print(np.shape(N1_wet_grid_subset))

    # plt.plot(N0_XC_subset,N0_YC_subset,'k.')
    # plt.plot(N1_XC_subset, N1_YC_subset, 'b.')
    # plt.show()

    # plt.imshow(N0_boundary_var_subset[0,:,1,:])
    # plt.show()

    # plt.imshow(N1_wet_grid_subset[:, 0, :])
    # plt.show()

    for timestep in range(nTimestepsOut):
        #if timestep%50==0:
        print('          Working on faces 1+2 in timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name)

        downscaled_field = df.downscale_3D_field(N0_XC_subset, N0_YC_subset,
                                                 N0_boundary_var_subset[timestep, :, :, :],
                                                 N0_wet_grid_subset, N0_wet_grid_on_N1_subset,
                                                 N1_XC_subset, N1_YC_subset, N1_wet_grid_subset)

        N1_boundary_var[timestep, :, :np.shape(N1_XC_faces[1])[1]+np.shape(N1_XC_faces[2])[1]] = downscaled_field[:,0,:]

        # C = plt.imshow(downscaled_field[:,0,:],aspect='equal')
        # plt.colorbar(C)
        # plt.title(timestep)
        # plt.show()

    #############################################################################################################
    # next, run the downscaling on faces 4 + 5

    # the data is only in the bottom 4 rows for the N0
    N0_XC_subset = np.concatenate([N0_XC_faces[4][:, -4:],N0_XC_faces[5][:, -4:]],axis=0)
    N0_YC_subset = np.concatenate([N0_YC_faces[4][:, -4:],N0_YC_faces[5][:, -4:]],axis=0)
    N0_boundary_var_subset = np.concatenate([N0_boundary_var_faces[4][:, :, :, -4:],
                                             N0_boundary_var_faces[5][:, :, :, -4:]],axis=2)

    N0_XC_subset[N0_XC_subset<0]+=360

    # the output will be on the bottom row of the N1
    N1_XC_subset = np.concatenate([N1_XC_faces[4][:, -3],N1_XC_faces[5][:, -3]],axis=0)
    N1_YC_subset = np.concatenate([N1_YC_faces[4][:, -3],N1_YC_faces[5][:, -3]],axis=0)

    N1_XC_subset = np.reshape(N1_XC_subset, (np.size(N1_XC_subset),1))
    N1_YC_subset = np.reshape(N1_YC_subset, (np.size(N1_YC_subset),1))

    N1_XC_subset[N1_XC_subset < 0] += 360

    N1_wet_grid_subset = np.concatenate([N1_wet_grid_3D[4][:, :, -3],
                                         N1_wet_grid_3D[5][:, :, -3]],axis=1)
    N1_wet_grid_subset = np.reshape(N1_wet_grid_subset,
                                    (np.shape(N1_wet_grid_subset)[0],np.shape(N1_wet_grid_subset)[1],1))

    N0_wet_grid_subset = []
    N0_wet_grid_on_N1_subset = []

    # print(np.shape(N0_XC_subset),np.shape(N0_YC_subset))
    # print(np.shape(N0_boundary_var_subset))
    # print(np.shape(N1_XC_subset),np.shape(N1_YC_subset))
    # print(np.shape(N1_wet_grid_subset))

    # plt.plot(N0_XC_subset,N0_YC_subset,'k.')
    # plt.plot(N1_XC_subset, N1_YC_subset, 'b.')
    # plt.show()
    #
    # plt.imshow(N0_boundary_var_subset[0,:,:,0])
    # plt.show()
    #
    # plt.imshow(N1_wet_grid_subset[:, :, 0])
    # plt.show()

    for timestep in range(nTimestepsOut):
        #if timestep%50==0:
        print('          Working on faces 4+5 in timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name)

        downscaled_field = df.downscale_3D_field(N0_XC_subset, N0_YC_subset,
                                                 N0_boundary_var_subset[timestep, :, :, :],
                                                 N0_wet_grid_subset, N0_wet_grid_on_N1_subset,
                                                 N1_XC_subset, N1_YC_subset, N1_wet_grid_subset)

        N1_boundary_var[timestep, :, np.shape(N1_XC_faces[1])[1] + np.shape(N1_XC_faces[2])[1]:] = downscaled_field[:,:,0]

        # C = plt.imshow(downscaled_field[:,:,0],aspect='equal')
        # plt.colorbar(C)
        # plt.title(timestep)
        # plt.show()

    return(N1_boundary_var)


def intepolate_boundary_condition_to_new_timestep(Lf_var_subset,var_name,downscale_factor):
    if downscale_factor!=2:
        raise ValueError('This interpolation only works for downscale_factor=2')

    if var_name == 'ETAN':
        interpolated_Lf_var_subset = np.zeros((downscale_factor * (np.shape(Lf_var_subset)[0] - 1) + 1,
                                               np.shape(Lf_var_subset)[1],
                                               np.shape(Lf_var_subset)[2]))
        for i in range(np.shape(Lf_var_subset)[0]-1):
            interpolated_Lf_var_subset[2*i,:,:] = Lf_var_subset[i,:,:]
            interpolated_grid = (Lf_var_subset[i,:,:]+Lf_var_subset[i+1,:,:])/2
            interpolated_Lf_var_subset[2*i+1,:,:] = interpolated_grid
        interpolated_Lf_var_subset[-1,:,:] = Lf_var_subset[-1,:,:]
    else:
        interpolated_Lf_var_subset = np.zeros((downscale_factor * (np.shape(Lf_var_subset)[0] - 1) + 1,
                                               np.shape(Lf_var_subset)[1],
                                               np.shape(Lf_var_subset)[2],
                                               np.shape(Lf_var_subset)[3]))
        for i in range(np.shape(Lf_var_subset)[0]-1):
            interpolated_Lf_var_subset[2*i,:,:,:] = Lf_var_subset[i,:,:,:]
            interpolated_grid = (Lf_var_subset[i,:,:,:]+Lf_var_subset[i+1,:,:,:])/2
            interpolated_Lf_var_subset[2*i+1,:,:,:] = interpolated_grid
        interpolated_Lf_var_subset[-1,:,:,:] = Lf_var_subset[-1,:,:,:]

    return(interpolated_Lf_var_subset)


def output_pbc_field(mask_name,var_name,N1_boundary_var,prescribe_layers):

    print(np.shape(N1_boundary_var))

    if var_name=='ETAN':
        if mask_name=='north':
            prescribe_grid = N1_boundary_var[:, -prescribe_layers:, :]
        if mask_name=='north_etan':
            prescribe_grid = N1_boundary_var[:, -prescribe_layers-1:, :]
        if mask_name=='south':
            prescribe_grid = N1_boundary_var[:, :prescribe_layers, :]
        if mask_name=='south_etan':
            prescribe_grid = N1_boundary_var[:, :prescribe_layers+1, :]
        if mask_name=='west':
            prescribe_grid = N1_boundary_var[:, :, :prescribe_layers]
        if mask_name == 'west_etan':
            prescribe_grid = N1_boundary_var[:, :, :prescribe_layers+1]
    else:
        if mask_name=='north':
            prescribe_grid = N1_boundary_var[:, :, -prescribe_layers:, :]
        if mask_name=='south':
            prescribe_grid = N1_boundary_var[:, :, :prescribe_layers, :]
        if mask_name=='west':
            prescribe_grid = N1_boundary_var[:, :, :, :prescribe_layers]

    if mask_name == 'north' or mask_name == 'north_etan':
        prescribe_grid = np.flip(prescribe_grid, axis=-2)

    if mask_name in ['west', 'east','west_etan', 'east_etan']:
        if var_name == 'ETAN':
            prescribe_grid = np.rot90(prescribe_grid, axes=(1, 2))
        else:
            prescribe_grid = np.rot90(prescribe_grid, axes=(2, 3))

    # if mask_name == 'east':
    #     prescribe_grid = np.flip(prescribe_grid, axis=-2)

    print('        - N1_BC_' + mask_name + '_' + var_name + '.bin output shape:' + str(np.shape(prescribe_grid)))

    print('        - N1_BC_'+mask_name+'_'+var_name+'.bin output shape:'+str(np.shape(prescribe_grid)))

    output_file = os.path.join('..', 'input', 'pbcs', 'N1_BC_' + mask_name + '_' + var_name + '.bin')
    prescribe_grid.ravel(order='C').astype('>f4').tofile(output_file)


########################################################################################################################

def create_bc_fields_via_interpolation(config_dir,ecco_path):

    LLC = 270
    N = 180

    llc = 1080
    n = 720

    print('Creating the bc fields for the N1_1080 model from the output of the N0_270 model')

    N0_Nr = 50
    N1_Nr = 90

    # this is the dir where the boundary conditions will be stored
    if 'obcs' not in os.listdir(os.path.join(config_dir,'N1_1080','input')):
        os.mkdir(os.path.join(config_dir,'N1_1080','input','obcs'))
    output_dir = os.path.join(config_dir,'N1_1080','input','obcs')

    print('    - Reading in the geometry of the N0_' + str(LLC) + ' domain')
    N0_XC_faces, N0_YC_faces = sf.read_sassie_n0_grid(ecco_path)

    print('    - Reading in the geometry of the N1_' + str(LLC) + ' domain')
    input_dir = os.path.join('..', 'input')
    N1_XC_faces, N1_YC_faces = sf.read_sassie_n1_grid(input_dir)

    print(np.shape(N0_XC_faces[1]))
    print(np.shape(N1_XC_faces[1]))

    N0_delR = sf.read_sassie_delR(domain_level=0)
    N1_delR = sf.read_sassie_delR(domain_level=1)

    var_names = ['THETA']#['UVEL','VVEL','SALT','THETA']

    masks = ['mask_41N','mask_41N_i1','mask_41N_i2','mask_41N_i3']

    years = np.arange(1993, 1994)

    for var_name in var_names:

        print(' + Working on variable '+var_name)

        print('   - Reading in the wet grids to use for in the interpolation')

    #     if var_name == 'VVEL':
    #         N1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),hFac='S')
    #         N0_wet_grid_3D_on_N1_faces = read_mask_faces_from_nc(
    #             os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='S')
    #     elif var_name == 'UVEL':
    #         N1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),hFac='W')
    #         N0_wet_grid_3D_on_N1_faces = read_mask_faces_from_nc(
    #             os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='W')
    #     else:
    #         N1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),hFac='C')
    #         N0_wet_grid_3D_on_N1_faces = read_mask_faces_from_nc(
    #             os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='C')
    #
    #     ##########################################################################################################
    #
    #     for year in years:
    #
    #         print('   + Working on data for year ' + str(year))
    #
    #         N0_boundary_var_faces = {}
    #         for face in range(1,6):
    #             N0_boundary_var_faces[face] = np.zeros((nTimestepsOut, N0_Nr,
    #                                                     np.shape(N0_XC_faces[face])[0],np.shape(N0_XC_faces[face])[1]))
    #
    #         print('     - Reading the mask variables onto the N0 geometry')
    #         for mask_name in masks:
    #
    #             print('        - Reading the '+mask_name+' mask reference faces, rows, and columns')
    #             nc_dict_file = os.path.join(config_dir,'N0_270','input','dv_mask_reference_dict.nc')
    #             mask_reference_name = mask_name
    #             faces, rows, cols = read_mask_reference_from_nc_dict(nc_dict_file, mask_reference_name)
    #
    #             print('        - Reading the ' + mask_name + ' mask variable onto the grid')
    #             N0_run_dir = os.path.join(config_dir, 'N0_270', 'run', 'dv')
    #             N0_boundary_var_faces = read_boundary_variable_onto_N0_grid(N0_run_dir, mask_name, var_name,year,
    #                                                                         nTimesteps, nTimestepsOut,
    #                                                                         faces, rows, cols, LLC, N,
    #                                                                         N0_Nr, N0_boundary_var_faces)
    #
    #         print('      - Extending the faces for an interpolation buffer')
    #         N0_boundary_var_faces_extended = {}
    #         for face in range(1, 6):
    #             N0_boundary_var_faces_extended[face] = sf.get_extended_var_grid_on_face_4D(N0_boundary_var_faces, face, N0_Nr,nTimestepsOut)
    #
    #         print('      - Reading the boundary to new vertical levels')
    #         N0_boundary_var_faces_extended = df.interpolate_var_grid_faces_to_new_depth_levels(N0_boundary_var_faces_extended,
    #                                                                                 N0_wet_grid_3D_on_N1_faces,
    #                                                                                 N0_delR, N1_delR)
    #
    #
    #         # for face in [1,2,4,5]:
    #         #     rows,cols = np.where(N0_boundary_var_faces_extended[face][0,0,:,:]!=0)
    #         #     print(face,np.min(rows),np.max(rows),np.min(cols),np.max(cols))
    #         #
    #         # print('    - Limiting geometry of the N0_' + str(LLC) + ' domain')
    #         # for face in range(1, 6):
    #         #     if face < 3:
    #         #         N0_XC_faces[face] = N0_XC_faces[face][-170:, :]
    #         #         N0_YC_faces[face] = N0_YC_faces[face][-170:, :]
    #         #         N0_boundary_var_faces[face] = N0_boundary_var_faces[face][:,:,-170:, :]
    #         #     if face > 3:
    #         #         N0_XC_faces[face] = N0_XC_faces[face][:, :170]
    #         #         N0_YC_faces[face] = N0_YC_faces[face][:, :170]
    #         #         N0_boundary_var_faces[face] = N0_boundary_var_faces[face][:,:,:, :170]
    #         #
    #         # N0_wet_grid_3D = []
    #         # N0_wet_grid_on_N1_3D = []
    #         #
    #         # print('  - Downscaling the output to the new boundary')
    #         # N1_boundary_var = downscale_N0_bc_field_to_N1(nTimestepsOut, var_name,
    #         #                                               N0_XC_faces,N0_YC_faces,N0_boundary_var_faces,
    #         #                                               N0_wet_grid_3D, N0_wet_grid_on_N1_3D,
    #         #                                               N1_XC_faces, N1_YC_faces, N1_wet_grid_3D)
    #         #
    #         # N1_boundary_var = downscale_N0_bc_field_to_N1(nTimestepsOut, var_name,
    #         #                                               N0_boundary_var_faces_extended, regridding_weights,
    #         #                                               N0_wet_grid_3D_on_N1_faces, N1_wet_grid_3D_faces,
    #         #                                               spread_horizontally=False)
    #         #
    #         # output_file = os.path.join(output_dir, 'N1_BC_' + var_name[0]+ '.bin')
    #         # N1_boundary_var.ravel('C').astype('>f4').tofile(output_file)
    #
    # # #         output_pbc_field(mask_name, var_name, N1_boundary_var, prescribe_layers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_path = args.ecco_path

    create_bc_fields_via_interpolation(config_dir, ecco_path)