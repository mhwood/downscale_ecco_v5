
import os
import numpy as np
import netCDF4 as nc4
import ast
from datetime import datetime
import simplegrid as sg
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata, interp1d
import sys
import argparse


def read_grid_tile_geometry(sf,config_dir,var_name,llc):

    input_dir = os.path.join(config_dir,'N1_1080','input')

    n = 680
    Nr = 90

    grid_file_dir = os.path.join('..', 'input')
    XC_faces = {}
    YC_faces = {}
    for i in [1, 2, 3, 4, 5]:
        if i < 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(input_dir, 'tile00' + str(i) + '.mitgrid'), llc, n)
        if i == 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(input_dir, 'tile00' + str(i) + '.mitgrid'), llc,llc)
        if i > 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(input_dir, 'tile00' + str(i) + '.mitgrid'), n, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    # for face in range(1,6):
    #     plt.subplot(1,2,1)
    #     C = plt.imshow(XC_faces[face], origin='lower')
    #     plt.colorbar(C)
    #     plt.subplot(1, 2, 2)
    #     C = plt.imshow(YC_faces[face], origin='lower')
    #     plt.colorbar(C)
    #     plt.title('face '+str(face))
    #     plt.show()

    file_names = ['AngleCS.data','AngleSN.data']
    if var_name == 'UVEL':
        file_names.append('hFacW.data')
    elif var_name == 'VVEL':
        file_names.append('hFacS.data')
    else:
        file_names.append('hFacC.data')

    for fn in range(len(file_names)):
        file_name = file_names[fn]

        grid_compact = np.fromfile(os.path.join(input_dir, file_name), '>f4')
        if 'hFac' in file_name:
            grid_compact = np.reshape(grid_compact, (Nr, 4 * n + llc, llc))
        else:
            grid_compact = np.reshape(grid_compact, (4 * n + llc, llc))

        # grid_compact, _, global_metadata = mds.rdmds(os.path.join(input_dir, file_name), returnmeta=True)
        grid_faces = sf.sassie_n1_compact_to_faces(grid_compact, llc, n)

        # if fn==2:
        #     for face in range(1, 6):
        #         C = plt.imshow(grid_faces[face][0, :, :], origin='lower')
        #         plt.colorbar(C)
        #         plt.title(file_name+', face ' + str(face))
        #         plt.show()

        if fn==0:
            AngleCS_faces = grid_faces
        if fn==1:
            AngleSN_faces = grid_faces
        if fn==2:
            hfac_faces = grid_faces

    delR = sf.read_sassie_delR(domain_level=1)

    return(XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, hfac_faces, delR)


def create_src_dest_list(config_dir,boundary,var_name):

    suffixes = []
    for file_name in os.listdir(os.path.join(config_dir,'N0_270','run','dv','mask_41N','SALT')):
        if file_name[-4:] == '.bin':
            suffix = '.'.join(file_name.split('.')[-2:])
            suffixes.append(suffix)

    print(suffixes)

    if '0000709561.bin' in suffixes:
        suffixes.pop(suffixes.index('0000709561.bin'))

    suffixes = sorted(suffixes)
    dest_files = []
    for i in range(len(suffixes)):
        dest_files.append('N1_1080_BC_'+boundary+'_' + var_name + '.'+str(1992+i)+'.bin')

    return(dest_files,suffixes)


def read_L0_boundary_variable_points(config_dir, dest_file, boundary, var_name, suffix,
                                     n_timesteps, llc, Nr, ecco_XC_faces, ecco_YC_faces,
                                     ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces):


    mask_names = ['mask_41N','mask_41N_i1','mask_41N_i2','mask_41N_i3']

    # print('        - Counting the number of points to read')
    n_points = 0
    nc_dict_file = os.path.join(config_dir, 'N0_270', 'input', 'dv_mask_reference_dict.nc')
    ds = nc4.Dataset(nc_dict_file)
    for mask_name in mask_names:
        grp = ds.groups[mask_name]
        faces = grp.variables['source_faces'][:]
        n_points += np.size(faces)
    ds.close()

    points = np.zeros((n_points, 2))
    values = np.zeros((n_timesteps, Nr, n_points))
    hfac_points = np.zeros((Nr, n_points))
    points_counted = 0

    for mask_name in mask_names:

        print('            - Reading the '+mask_name+' mask to reference the variable to the llc grid')
        nc_dict_file = os.path.join(config_dir, 'N0_270', 'input', 'dv_mask_reference_dict.nc')
        ds = nc4.Dataset(nc_dict_file)
        grp = ds.groups[mask_name]
        faces = grp.variables['source_faces'][:]
        rows = grp.variables['source_rows'][:]
        cols = grp.variables['source_cols'][:]
        ds.close()

        # print('          - Reading the BC data grid')
        N = np.size(faces)

        if var_name in ['UVEL','VVEL']:
            u_var_file = os.path.join(config_dir, 'N0_270', 'run', 'dv', mask_name, 'UVEL', mask_name + '_UVEL' + suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            u_var_grid = np.reshape(u_var_grid,(n_timesteps,Nr,N))
            v_var_file = os.path.join(config_dir, 'N0_270', 'run', 'dv', mask_name, 'VVEL', mask_name + '_VVEL' + suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            v_var_grid = np.reshape(v_var_grid, (n_timesteps, Nr, N))
        else:
            file_path = os.path.join(config_dir, 'N0_270', 'run', 'dv', mask_name, var_name,
                                     mask_name + '_' + var_name + suffix)
            var_grid = np.fromfile(file_path, '>f4')
            var_grid = np.reshape(var_grid, (n_timesteps, Nr, N))

        for n in range(N):
            if faces[n] != 0:
                points[points_counted+n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                points[points_counted+n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                hfac_points[:, points_counted+n] = ecco_hFacC_faces[faces[n]][:, rows[n], cols[n]]
                if var_name in ['UVEL', 'VVEL']:
                    angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                    angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                    if var_name=='UVEL':
                        zonal_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                        for k in range(np.shape(zonal_velocity)[1]):
                            zonal_velocity[:,k] = angle_cos * u_var_grid[:, k, n] - angle_sin * v_var_grid[:, k, n]
                        values[:, :, points_counted+n] = zonal_velocity
                    if var_name=='VVEL':
                        meridional_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                        for k in range(np.shape(meridional_velocity)[1]):
                            meridional_velocity[:,k] = angle_sin * u_var_grid[:, k, n] + angle_cos * v_var_grid[:, k, n]
                        values[:, :, points_counted+n] = meridional_velocity
                else:
                    values[:, :, points_counted+n] = var_grid[:, :, n]

        points_counted += N

    return(points,values,hfac_points)

def subset_tile_geometry_to_boundary(boundary, var_name, face_number,
                                     XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, hfac_faces):

    # plt.subplot(1,2,1)
    # C = plt.imshow(XC_faces[face_number],origin='lower')
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(YC_faces[face_number], origin='lower')
    # plt.colorbar(C)
    # plt.show()

    # subset to the boundary

    if boundary=='south':
        XC_subset = XC_faces[face_number][:1, :]
        YC_subset = YC_faces[face_number][:1, :]
        AngleCS_subset = AngleCS_faces[face_number][:1, :]
        AngleSN_subset = AngleSN_faces[face_number][:1, :]
        if var_name=='VVEL':
            hFac_subset = hfac_faces[face_number][:,1:2, :]
        else:
            hFac_subset = hfac_faces[face_number][:, :1, :]

    if boundary=='east':
        XC_subset = XC_faces[face_number][:, -1:]
        YC_subset = YC_faces[face_number][:, -1:]
        AngleCS_subset = AngleCS_faces[face_number][:, -1:]
        AngleSN_subset = AngleSN_faces[face_number][:, -1:]
        hFac_subset = hfac_faces[face_number][:,:, -1:]

    return(XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset)

def interpolate_L0_boundary_variables_to_new_depth_levels(var_grid,wet_grid,delR_in,delR_out):

    Z_bottom_in = np.cumsum(delR_in)
    Z_top_in = np.concatenate([np.array([0]), Z_bottom_in[:-1]])
    Z_in = (Z_bottom_in + Z_top_in) / 2

    Z_bottom_out = np.cumsum(delR_out)
    Z_top_out = np.concatenate([np.array([0]), Z_bottom_out[:-1]])
    Z_out = (Z_bottom_out + Z_top_out) / 2

    new_var_grid = np.zeros((np.shape(var_grid)[0], np.size(delR_out), np.shape(var_grid)[2]))
    for timestep in range(np.shape(new_var_grid)[0]):
        for j in range(np.shape(var_grid)[2]):
            test_profile = var_grid[timestep, :, j]
            if np.sum(test_profile != 0) > 1:
                set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                          bounds_error=False, fill_value=np.nan)
                new_profile = set_int_linear(Z_out)

                new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                    first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                    bottom_value = new_profile[~np.isnan(new_profile)][-1]
                    new_profile[np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                new_var_grid[timestep, :, j] = new_profile

            if np.sum(test_profile == 0) == 1:
                new_var_grid[timestep, 0, j] = var_grid[timestep, 0, j]

    if np.shape(wet_grid)[0]!=len(delR_out):
        new_wet_grid = np.zeros((np.size(delR_out), np.shape(wet_grid)[1]))
        for i in range(np.shape(wet_grid)[1]):
            test_profile = wet_grid[:, i]
            if np.sum(test_profile != 0) > 1:
                set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                          bounds_error=False, fill_value=np.nan)
                new_profile = set_int_linear(Z_out)

                new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                    first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                    bottom_value = new_profile[~np.isnan(new_profile)][-1]
                    new_profile[np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                new_wet_grid[:, i] = new_profile

            if np.sum(test_profile == 0) == 1:
                new_wet_grid[0, i] = wet_grid[0, i]
        new_wet_grid[np.isnan(new_wet_grid)] = 0
        new_wet_grid = np.round(new_wet_grid).astype(int)

        new_var_grid[np.isnan(new_var_grid)] = 0
    else:
        new_wet_grid = wet_grid


    return(new_var_grid,new_wet_grid)


def create_bc_fields_via_interpolation(config_dir,ecco_dir):

    sys.path.insert(1, os.path.join(config_dir, 'utils'))
    import downscale_functions as df
    import ecco_functions as ef
    import sassie_functions as sf

    LLC = 270

    llc = 1080
    n = 680
    n_timesteps = 73

    Nr_in = 50
    Nr_out = 90

    for var_name in ['THETA','SALT','UVEL','VVEL']: #'UVEL','VVEL','THETA','SALT'
        print('    - Creating the ' + var_name + ' BC files from ECCOv5 data')

        # step 0: get the model domain
        print('    - Reading in the model geometry')
        XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, hfac_faces, delR = \
            read_grid_tile_geometry(sf,config_dir,var_name,llc)

        # step 1: get the ecco faces geometry
        print('    - Reading in the ECCO geometry')
        ecco_dir = '/Users/michwood/Documents/Research/Projects/Ocean_Modeling/ECCO'
        # ecco_XC, ecco_YC, ecco_AngleCS, ecco_AngleSN, ecco_hfacC, ecco_delR = \
        #     ef.read_ecco_grid_geometry(ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations)
        ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = ef.read_ecco_geometry_to_faces(ecco_dir, LLC)
        ecco_delR = sf.read_sassie_delR(domain_level=0)

        N0_run_dir = os.path.join(config_dir, 'N0_270', 'run')

        ####################################################################################################################
        # Loop through the boundaries

        for boundary in ['south','east']:

            if boundary == 'south':
                faces_to_compute = [1,2]
            if boundary == 'east':
                faces_to_compute = [4,5]

            print('    - Running the downscale routine for the '+boundary+' boundary')

            if boundary not in os.listdir(os.path.join(config_dir,'N1_1080','input','obcs')):
                os.mkdir(os.path.join(config_dir,'N1_1080','input','obcs',boundary))
            if var_name not in os.listdir(os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary)):
                os.mkdir(os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary, var_name))

            dest_files, source_file_suffixes = create_src_dest_list(config_dir, boundary, var_name)

            # print(dest_files)
            # print(source_file_suffixes)

            for d in range(len(dest_files)-3,len(dest_files)):

                dest_file = dest_files[d]
                print('        - Working on '+dest_file)
                suffix = '.'+source_file_suffixes[d]
                print('        - Reading from suffix '+suffix)

                print('            - Reading in the L0_540 diagnostics_vec output')
                L0_boundary_points, L0_boundary_values, L0_boundary_point_hFacC = \
                    read_L0_boundary_variable_points(config_dir, dest_file, boundary, var_name, suffix,
                                                     n_timesteps, llc, Nr_in, ecco_XC_faces, ecco_YC_faces,
                                                     ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces)

                # print('                - L0_540 extent: x:'+str(np.min(L0_boundary_points[:,0]))+' '+str(np.max(L0_boundary_points[:,0])))
                # print('                - L0_540 extent: y:' + str(np.min(L0_boundary_points[:, 1])) + ' ' + str(np.max(L0_boundary_points[:, 1])))

                print('        - Downscaling the timesteps to be stored in file ' + str(dest_file))

                output_points_added = 0
                output_grid = np.zeros((n_timesteps, Nr_out, 2*n + 3*llc))

                L0_boundary_values = L0_boundary_values[:, :, L0_boundary_points[:, 0] != 0]
                L0_boundary_point_hFacC = L0_boundary_point_hFacC[:, L0_boundary_points[:, 0] != 0]
                L0_boundary_points = L0_boundary_points[L0_boundary_points[:, 0] != 0, :]

                L0_boundary_point_mask = np.copy(L0_boundary_point_hFacC)
                L0_boundary_point_mask[L0_boundary_point_mask > 0] = 1

                L0_boundary_values,L0_boundary_point_mask =\
                    interpolate_L0_boundary_variables_to_new_depth_levels(L0_boundary_values, L0_boundary_point_mask,
                                                                          ecco_delR, delR)

                # plt.plot(L0_boundary_points[:,0],L0_boundary_points[:,1],'g.')
                # plt.show()

                for face_number in range(1,6):
                    # face_number = [1,2,4,5][fn]

                    if face_number in faces_to_compute:
                        print('            - Downscaling data for face number '+str(face_number))

                        XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset = \
                            subset_tile_geometry_to_boundary(boundary, var_name, face_number,
                                                             XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, hfac_faces)

                        # print('                - XC_subset extent: x:' + str(np.min(XC_subset)) + ' ' + str(np.max(XC_subset)))
                        # print('                - YC_subset extent: y:' + str(np.min(YC_subset)) + ' ' + str(np.max(YC_subset)))

                        mask_subset = np.copy(hFac_subset)
                        mask_subset[mask_subset>0]=1

                        # plt.imshow(mask_subset[:,0,:])
                        # plt.show()

                        # plt.plot(L0_boundary_points[:, 0], L0_boundary_points[:, 1], 'g.')
                        # plt.plot(XC_subset.ravel(),YC_subset.ravel(),'k.')
                        # plt.show()

                        for timestep in range(np.shape(L0_boundary_values)[0]):
                            if timestep%5 == 0:
                                print('                - Downscaling timestep '+str(timestep)+' of '+str(np.shape(L0_boundary_values)[0]))

                            # if timestep == 0:
                            #     print('            - L0_boundary_points shape: '+str(np.shape(L0_boundary_points)))
                            #     print('            - L0_boundary_values shape: ' + str(np.shape(L0_boundary_values)))
                            #     print('            - L0_boundary_point_mask shape: ' + str(np.shape(L0_boundary_point_mask)))
                            #     print('            - XC_subset shape: ' + str(np.shape(XC_subset)))
                            #     print('            - YC_subset shape: ' + str(np.shape(YC_subset)))
                            #     print('            - mask_subset shape: ' + str(np.shape(mask_subset)))

                                # plt.subplot(1,2,1)
                                # plt.imshow(L0_boundary_values[0,:,:],aspect=15)
                                # plt.subplot(1, 2, 2)
                                # plt.imshow(L0_boundary_point_mask,aspect=15)
                                # plt.show()

                            interp_field = df.downscale_3D_boundary_points(L0_boundary_points, L0_boundary_values[timestep,:,:], L0_boundary_point_mask,
                                                         XC_subset, YC_subset, mask_subset,
                                                         mean_vertical_difference=0, fill_downward=True,
                                                         printing=False, remove_zeros=False)

                            if boundary in ['north', 'south']:
                                output_grid[timestep, :, output_points_added:output_points_added + np.shape(interp_field)[2]] = interp_field[:,0,:]
                                # if timestep
                                # output_points_added += np.shape(interp_field)[2]
                            else:
                                output_grid[timestep, :, output_points_added:output_points_added + np.shape(interp_field)[1]] = interp_field[:,:,0]
                                # output_points_added += np.shape(interp_field)[1]

                            # if timestep==0:
                            #     if boundary in ['north', 'south']:
                            #         plt.subplot(1,2,1)
                            #         plt.imshow(mask_subset[:,0,:],aspect = 5)
                            #         plt.subplot(1, 2, 2)
                            #         plt.imshow(interp_field[:,0,:],aspect = 5)
                            #     else:
                            #         plt.subplot(1, 2, 1)
                            #         plt.imshow(mask_subset[:, :, 0],aspect = 5)
                            #         plt.subplot(1, 2, 2)
                            #         plt.imshow(interp_field[:, :, 0],aspect = 5)
                            #     plt.show()

                    if face_number in [1,2] and boundary=='east':
                        output_points_added += n
                    if face_number in [4,5] and boundary=='south':
                        output_points_added += n
                    if face_number in [1,2] and boundary=='south':
                        output_points_added += llc
                    if face_number in [4,5] and boundary=='east':
                        output_points_added += llc
                    if face_number in [3]:
                        output_points_added += llc

                # if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                #     plt.plot(output_grid[0, :])
                #     plt.show()
                # else:
                #     plt.imshow(output_grid[0, :, :],aspect = 7,vmin=-0.1,vmax=0.1,cmap='seismic')
                #     plt.title(boundary+', '+var_name)
                #     plt.show()

                if var_name in ['UVEL','VVEL']:
                    output_file = os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary, var_name, dest_file[:-4]+'_rotated.bin')
                else:
                    output_file = os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary, var_name, dest_file)
                output_grid.ravel(order='C').astype('>f4').tofile(output_file)








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