
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys


def get_dest_file_list(boundary, var_name, timesteps_per_file, Nr, llc, n, start_year, final_year):
    dest_files = []
    dest_file_shapes = {}
    total_timesteps = 0
    for year in range(start_year,final_year+1):
        dest_file = 'N1_1080_BC_' + boundary + '_' + var_name + '.' + str(year) + '.bin'
        dest_files.append(dest_file)
        total_timesteps += timesteps_per_file
        dest_file_shapes[dest_file] = (timesteps_per_file, Nr, 3*llc+2*n)

    return(dest_files,dest_file_shapes,total_timesteps)

def read_grid_angle_faces(sf,config_dir,llc,n):

    input_dir = os.path.join(config_dir,'N1_1080','input')

    file_names = ['AngleCS.data','AngleSN.data']

    for fn in range(len(file_names)):
        file_name = file_names[fn]

        grid_compact = np.fromfile(os.path.join(input_dir, file_name), '>f4')
        grid_compact = np.reshape(grid_compact, (4 * n + llc, llc))
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

    return(AngleCS_faces, AngleSN_faces)

def subset_grid_angle_faces_to_boundary(boundary, llc, n,
                                        AngleCS_faces, AngleSN_faces):

    AngleCS_subset = np.zeros((llc * 3 + n * 2,))
    AngleSN_subset = np.zeros((llc * 3 + n * 2,))

    if boundary=='south':
        AngleCS_subset[:1080] = AngleCS_faces[1][:1,:]
        AngleCS_subset[1080:2*1080] = AngleCS_faces[2][:1, :]
        AngleCS_subset[2 * 1080:3 * 1080] = AngleCS_faces[3][:1, :]
        AngleCS_subset[3*1080:3*1080+680] = AngleCS_faces[4][:1, :]
        AngleCS_subset[3*1080+680:] = AngleCS_faces[5][:1, :]

        AngleSN_subset[:1080] = AngleSN_faces[1][:1, :]
        AngleSN_subset[1080:2 * 1080] = AngleSN_faces[2][:1, :]
        AngleSN_subset[2*1080:3 * 1080] = AngleSN_faces[3][:1, :]
        AngleSN_subset[3 * 1080:3 * 1080 + 680] = AngleSN_faces[4][:1, :]
        AngleSN_subset[3 * 1080 + 680:] = AngleSN_faces[5][:1, :]

    if boundary=='east':
        AngleCS_subset[:680] = AngleCS_faces[1][:,-1]
        AngleCS_subset[680:2*680] = AngleCS_faces[2][:, -1]
        AngleCS_subset[2 * 680:2 * 680 + 1080] = AngleCS_faces[3][:, -1]
        AngleCS_subset[2*680+1080:2*680+2*1080] = AngleCS_faces[4][:, -1]
        AngleCS_subset[2*680+2*1080:] = AngleCS_faces[5][:, -1]

        AngleSN_subset[:680] = AngleSN_faces[1][:, -1]
        AngleSN_subset[680:2 * 680] = AngleSN_faces[2][:, -1]
        AngleSN_subset[2 * 680:2 * 680 + 1080] = AngleSN_faces[3][:, -1]
        AngleSN_subset[2 * 680+1080:2 * 680 + 2*1080] = AngleSN_faces[4][:, -1]
        AngleSN_subset[2 * 680 + 2*1080:] = AngleSN_faces[5][:, -1]

    return(AngleCS_subset, AngleSN_subset)

def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
    uvel = np.zeros_like(zonal_vel)
    vvel = np.zeros_like(meridional_vel)
    for k in range(np.shape(uvel)[0]):
        uvel[k, :] = angle_cos * zonal_vel[k, :] + angle_sin * meridional_vel[k, :]
        vvel[k, :] = -1 * angle_sin * zonal_vel[k, :] + angle_cos * meridional_vel[k, :]
    return (uvel, vvel)

def rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset):
    rotated_u = np.zeros_like(u_var_grid)
    rotated_v = np.zeros_like(v_var_grid)
    for timestep in range(np.shape(u_var_grid)[0]):
        u, v = rotate_velocity_vectors_to_domain(AngleCS_subset, AngleSN_subset,
                                                 u_var_grid[timestep, :, :], v_var_grid[timestep, :, :])
        rotated_u[timestep, :, :] = u
        rotated_v[timestep, :, :] = v
    if var_name[0]=='U':
        return(rotated_u)
    if var_name[0]=='V':
        return(rotated_v)

def stack_files_to_one(config_dir, boundary, var_name,
                       dest_files, dest_file_shapes, total_timesteps,
                       AngleCS_subset, AngleSN_subset):
    depth_levels = dest_file_shapes[dest_files[0]][1]
    n_points = dest_file_shapes[dest_files[0]][2]

    # the 2 is added because we will duplicate the first and last field
    output_grid = np.zeros((total_timesteps + 2, depth_levels, n_points))
    timesteps_added = 1
    for dest_file in dest_files:
        print('      - Adding timesteps from file ' + dest_file)
        if var_name in ['UVEL','VVEL']:
            if var_name=='UVEL':
                u_dest_file = dest_file[:-4]+'_rotated.bin'
                v_dest_file = u_dest_file.replace('UVEL','VVEL')
            if var_name=='VVEL':
                v_dest_file = dest_file[:-4]+'_rotated.bin'
                u_dest_file = v_dest_file.replace('VVEL','UVEL')
            u_var_grid = np.fromfile(os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary, 'U'+var_name[1:], u_dest_file),'>f4')
            v_var_grid = np.fromfile(os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary, 'V'+var_name[1:], v_dest_file),'>f4')
            u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
            v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])
            print('          - Rotating the grid to the domain orientation')
            var_grid = rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset)
        else:
            var_grid = np.fromfile(os.path.join(config_dir, 'N1_1080', 'input', 'obcs', boundary, var_name, dest_file), '>f4')
            var_grid = np.reshape(var_grid, dest_file_shapes[dest_file])

        print('          - Timesteps added to levels ' + str(timesteps_added) + ' through ' + str(timesteps_added + np.shape(var_grid)[0]))
        print('          - The mean value on this day is ' + str(np.mean(var_grid)))
        print('          - There are ' + str(np.sum(np.isnan(var_grid[0, :, :]))) + ' nan values')
        # rows = np.where(np.isnan(var_grid[0,:,:]))
        # print(rows)
        # print(cols)
        output_grid[timesteps_added:timesteps_added + np.shape(var_grid)[0], :, :] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    # here we duplicate the first and last field
    # this is done so that all timesteps can be interpolated by the model
    output_grid[0, :, :] = output_grid[1, :, :]
    output_grid[-1, :, :] = output_grid[-2, :, :]
    return (output_grid)


def combine_BC_files(config_dir):
    sys.path.insert(1, os.path.join(config_dir, 'utils'))
    import sassie_functions as sf

    start_year = 2014
    final_year = 2021

    llc = 1080
    n = 680
    timesteps_per_file = 73
    Nr = 90

    for var_name in ['THETA','SALT','UVEL','VVEL']:#

        for boundary in ['east','south']:


            print('  Combining the daily files for ' + var_name+' on the '+boundary+' boundary')

            print('      - Determining the size of the input/output files')
            dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(boundary, var_name, timesteps_per_file,
                                                                               Nr, llc, n, start_year, final_year)

            if var_name in ['UVEL','VVEL',]:
                print('      - Reading the grid orientation')
                AngleCS_faces, AngleSN_faces = read_grid_angle_faces(sf,config_dir,llc,n)
                AngleCS_subset, AngleSN_subset = subset_grid_angle_faces_to_boundary(boundary, llc, n,
                                                                                     AngleCS_faces, AngleSN_faces)
            else:
                AngleCS_subset = AngleSN_subset = []

            print('      - Stacking all of the files into a big global file')
            output_grid = stack_files_to_one(config_dir, boundary, var_name,
                                              dest_files, dest_file_shapes, total_timesteps,
                                              AngleCS_subset, AngleSN_subset)


            output_file = os.path.join(config_dir, 'N1_1080', 'input', 'obcs', 'N1_1080_BC_'+boundary+'_' + var_name + '.bin')
            print('      - Outputting to ' + output_file)
            output_grid.ravel('C').astype('>f4').tofile(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L2, and L3 configurations are stored.",
                        dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    combine_BC_files(config_dir)


