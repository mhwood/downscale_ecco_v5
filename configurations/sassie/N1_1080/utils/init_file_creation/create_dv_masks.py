
import os
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import netCDF4 as nc4
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

########################################################################################################################


def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list):

    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names_list)):
        # print(mask_names_list[m],np.shape(all_mask_dicts[m]))

        grp = ds.createGroup(mask_names_list[m])
        grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
        var = grp.createVariable('source_faces', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 0].astype(int)
        var = grp.createVariable('source_rows', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 1].astype(int)
        var = grp.createVariable('source_cols', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 2].astype(int)

    ds.close()

def create_mask_files():

    print('Creating the dv masks N1_1080 configuration')

    llc = 1080
    n = 680

    if 'plots' not in os.listdir('..'):
        os.mkdir(os.path.join('..', 'plots'))
    if 'init_files' not in os.listdir(os.path.join('..', 'plots')):
        os.mkdir(os.path.join('..', 'plots', 'init_files'))

    f = open(os.path.join('..', '..', 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L1_size = size_dict['N1_1080']
    n_rows = L1_size[0]
    n_cols = L1_size[1]
    n_rows = 1080 + 4 * 680

    bathy_file = os.path.join('..', 'input', 'bathy_N1_walls')
    bathy_grid = np.fromfile(bathy_file, '>f4')
    bathy_grid_compact = np.reshape(bathy_grid, (n_rows, n_cols))
    bathy_grid_faces = sf.sassie_n1_compact_to_faces(bathy_grid_compact, llc, n)



    ####################################################################################################################
    # make the beafort mask

    beaufort_points = []

    mask_faces = {}
    for face in range(1, 6):
        mask_faces[face] = np.zeros_like(bathy_grid_faces[face])

    counter = 1
    for col in range(80,-1,-1):
        for row in range(480,1000):
            if bathy_grid_faces[4][row,col]<0:
                beaufort_points.append([4,row,col])
                mask_faces[4][row,col] = counter
                counter +=1

    for col in range(1079,680,-1):
        for row in range(480,1000):
            if bathy_grid_faces[3][row,col]<0:
                beaufort_points.append([3,row,col])
                mask_faces[3][row,col] = counter
                counter +=1

    beaufort_points = np.array(beaufort_points)

    all_mask_dicts = [beaufort_points]
    mask_names_list = ['beaufort']

    # convert the faces to compact form
    compact_mask = sf.sassie_n1_faces_to_compact(mask_faces, llc, n)

    # output the file
    output_file = os.path.join('..', 'input', 'dv', 'beaufort_mask.bin')
    compact_mask.ravel('C').astype('>f4').tofile(output_file)

    ####################################################################################################################
    # make the N2 mask

    N2_boundary_points = []
    N2_llc = 2160
    N2_n = 320

    mask_faces = {}
    for face in range(1, 6):
        mask_faces[face] = np.zeros_like(bathy_grid_faces[face])

    # # decided to do this manually so didnt finish this section
    # for face in [1]:
    #     n2_face_file = os.path.join('..','..','N2_2160','input','tile'+'{:03d}'.format(face)+'.mitgrid')
    #     n2_face_grid = np.fromfile(n2_face_file,'>f8')
    #     if face in [1,2]:
    #         n2_face_grid = np.reshape()

    counter = 1
    for col in range(1080):
        for row in range(-162,-159):
            if bathy_grid_faces[1][row, col] < 0:
                N2_boundary_points.append([1, row, col])
                mask_faces[1][row, col] = counter
                counter += 1

    for col in range(1080):
        for row in range(-162,-159):
            if bathy_grid_faces[2][row, col] < 0:
                N2_boundary_points.append([2, row, col])
                mask_faces[2][row, col] = counter
                counter += 1

    for row in range(1080):
        for col in range(159,162):
            if bathy_grid_faces[4][row, col] < 0:
                N2_boundary_points.append([4, row, col])
                mask_faces[4][row, col] = counter
                counter += 1

    for row in range(1080):
        for col in range(159,162):
            if bathy_grid_faces[5][row, col] < 0:
                N2_boundary_points.append([5, row, col])
                mask_faces[5][row, col] = counter
                counter += 1

    N2_boundary_points = np.array(N2_boundary_points)

    all_mask_dicts.append(N2_boundary_points)
    mask_names_list.append('N2_boundary')

    # convert the faces to compact form
    compact_mask = sf.sassie_n1_faces_to_compact(mask_faces, llc, n)

    # output the file
    output_file = os.path.join('..', 'input', 'dv', 'N2_boundary_mask.bin')
    compact_mask.ravel('C').astype('>f4').tofile(output_file)

    ##############################################################################################

    output_dir = os.path.join('..', 'input')
    output_file_name = 'N1_1080_dv_mask_ref.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list)




if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    #
    # parser.add_argument("-e", "--ecco_directory", action="store",
    #                     help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
    #                     type=str, required=False, default='../../ECCO')
    #
    # args = parser.parse_args()
    # ecco_path = args.ecco_path

    create_mask_files()