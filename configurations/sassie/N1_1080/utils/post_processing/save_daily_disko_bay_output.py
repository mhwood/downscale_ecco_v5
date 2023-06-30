
import os
import time
from datetime import datetime, timedelta
import simplegrid as sg
import matplotlib.path as mpath
import netCDF4 as nc4

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean.cm as cm
import argparse
from matplotlib.patches import Rectangle


def field_number_to_set_and_name(field_number):

    # this is a dictionary referencing whats in the diagnostics file

    # 3D files
    field_number_dict = {1: ['ocean_state_3D_day_mean', 'SALT', 0, False],
                            2: ['ocean_state_3D_day_mean', 'THETA', 1, False],
                            3: ['ocean_vel_day_mean', 'UVEL', 0, True],
                            4: ['ocean_vel_day_mean', 'VVEL', 1, True]}

    field_set = field_number_dict[field_number][0]
    field_name = field_number_dict[field_number][1]
    field_set_index = field_number_dict[field_number][2]
    is_vector = field_number_dict[field_number][3]

    return(field_set,field_name,field_set_index,is_vector)

def get_filename_date_pairs(config_dir, field_set):

    def file_iter_to_date(iter_number):
        timestep = 120
        seconds_elapsed = iter_number*timestep
        date = datetime(1992,1,1) + timedelta(seconds=seconds_elapsed)
        return(date.year, date.month, date.day)

    diag_dir = os.path.join(config_dir,'N1_1080','run','diags',field_set)
    all_file_names = os.listdir(diag_dir)

    file_prefixes = []
    dates = []

    for file_name in all_file_names:
        if file_name[0]!='.' and file_name[-5:]=='.data':
            file_prefix = file_name[:-5]
            file_prefixes.append(file_prefix)
            iter_number = int(file_prefix.split('.')[-1])
            year, month, day = file_iter_to_date(iter_number)
            dates.append(str(year)+'{:02d}'.format(month)+'{:02d}'.format(day))

    return(file_prefixes, dates)

def read_faces_geometry(config_dir, is_vector):
    llc = 1080
    N = 680

    grid_file_dir = os.path.join(config_dir,'N1_1080','input')
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

    AngleCS_faces = {}
    AngleSN_faces = {}
    if is_vector:
        AngleCS_grid = np.fromfile(os.path.join(grid_file_dir, 'AngleCS.data'), '>f4')
        N = int(np.size(AngleCS_grid) / (102600 * 40))
        AngleCS_grid = np.reshape(AngleCS_grid, (N, 102600, 40))

        face_1 = AngleCS_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1 = face_1[0, :, :]
        AngleCS_faces[1] = face_1

        face_2 = AngleCS_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2 = face_2[0, :, :]
        AngleCS_faces[2] = face_2

        face_3 = AngleCS_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3 = face_3[0, :, :]
        AngleCS_faces[3] = face_3

        face_4 = AngleCS_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4 = face_4[0, :, :]
        AngleCS_faces[4] = face_4

        face_5 = AngleCS_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5 = face_5[0, :, :]
        AngleCS_faces[5] = face_5

        AngleSN_grid = np.fromfile(os.path.join(grid_file_dir, 'AngleSN.data'), '>f4')
        N = int(np.size(AngleSN_grid) / (102600 * 40))
        AngleSN_grid = np.reshape(AngleSN_grid, (N, 102600, 40))

        face_1 = AngleSN_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1 = face_1[0, :, :]
        AngleSN_faces[1] = face_1

        face_2 = AngleSN_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2 = face_2[0, :, :]
        AngleSN_faces[2] = face_2

        face_3 = AngleSN_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3 = face_3[0, :, :]
        AngleSN_faces[3] = face_3

        face_4 = AngleSN_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4 = face_4[0, :, :]
        AngleSN_faces[4] = face_4

        face_5 = AngleSN_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5 = face_5[0, :, :]
        AngleSN_faces[5] = face_5

    delR = np.array([1.00, 1.14, 1.30, 1.49, 1.70,
                         1.93, 2.20, 2.50, 2.84, 3.21,
                         3.63, 4.10, 4.61, 5.18, 5.79,
                         6.47, 7.20, 7.98, 8.83, 9.73,
                         10.69, 11.70, 12.76, 13.87, 15.03,
                         16.22, 17.45, 18.70, 19.97, 21.27,
                         22.56, 23.87, 25.17, 26.46, 27.74,
                         29.00, 30.24, 31.45, 32.65, 33.82,
                         34.97, 36.09, 37.20, 38.29, 39.37,
                         40.45, 41.53, 42.62, 43.73, 44.87,
                         46.05, 47.28, 48.56, 49.93, 51.38,
                         52.93, 54.61, 56.42, 58.38, 60.53,
                         62.87, 65.43, 68.24, 71.33, 74.73,
                         78.47, 82.61, 87.17, 92.21, 97.79,
                         103.96, 110.79, 118.35, 126.73, 136.01,
                         146.30, 157.71, 170.35, 184.37, 199.89,
                         217.09, 236.13, 257.21, 280.50, 306.24,
                         334.64, 365.93, 400.38, 438.23, 479.74, ])

    Z_bottom_out = np.cumsum(delR)
    Z_top_out = np.concatenate([np.array([0]), Z_bottom_out[:-1]])
    Z_out = (Z_bottom_out + Z_top_out) / 2


    return(XC_faces,YC_faces,AngleCS_faces,AngleSN_faces,Z_out)

def subset_geometry_to_Disko(XC_faces, YC_faces):

    face = 5

    # these lines calculate the bounds, only ran this once

    # min_lon = -54.5
    # min_lat = 68.3
    # max_lon = -50
    # max_lat = 70.1
    #
    # ll_dist = (XC_faces[face]-min_lon)**2 + (YC_faces[face]-min_lat)**2
    # ll_row, ll_col = np.where(ll_dist==np.min(ll_dist))
    # ll_row = ll_row[0]
    # ll_col = ll_col[0]
    # print(ll_row,ll_col)
    # print(XC_faces[face][ll_row,ll_col],YC_faces[face][ll_row,ll_col])
    #
    # lr_dist = (XC_faces[face] - max_lon) ** 2 + (YC_faces[face] - min_lat) ** 2
    # lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    # lr_row = lr_row[0]
    # lr_col = lr_col[0]
    # print(lr_row, lr_col)
    # print(XC_faces[face][lr_row, lr_col], YC_faces[face][lr_row, lr_col])
    #
    # ur_dist = (XC_faces[face] - max_lon) ** 2 + (YC_faces[face] - max_lat) ** 2
    # ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))
    # ur_row = ur_row[0]
    # ur_col = ur_col[0]
    # print(ur_row, ur_col)
    # print(XC_faces[face][ur_row, ur_col], YC_faces[face][ur_row, ur_col])
    #
    # ul_dist = (XC_faces[face] - min_lon) ** 2 + (YC_faces[face] - max_lat) ** 2
    # ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))
    # ul_row = ul_row[0]
    # ul_col = ul_col[0]
    # print(ul_row, ul_col)
    # print(XC_faces[face][ul_row, ul_col], YC_faces[face][ul_row, ul_col])
    #
    # min_col = np.min([ll_col,lr_col,ur_col,ul_col])
    # max_col = np.max([ll_col, lr_col, ur_col, ul_col])
    # min_row = np.min([ll_row, lr_row, ur_row, ul_row])
    # max_row = np.max([ll_row, lr_row, ur_row, ul_row])

    # these were calculated from the above lines
    min_row = 829
    max_row = 896
    min_col = 0
    max_col = 62

    return(face, min_row,max_row,min_col,max_col)

def read_field_to_disko_subset(config_dir, field_set, file_name, field_set_index,
                               face, min_row,max_row,min_col,max_col):

    file_path = os.path.join(config_dir,'N1_1080','run','diags',field_set,file_name+'.data')

    dim = 3

    if dim==2:
        diag_grid = np.fromfile(file_path, '>f4')
        N = int(np.size(diag_grid) / (102600 * 40))
        diag_grid = np.reshape(diag_grid, (N, 102600, 40))

        face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1 = face_1[field_set_index, :, :]

        face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2 = face_2[field_set_index, :, :]

        face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3 = face_3[field_set_index, :, :]

        face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4 = face_4[field_set_index, :, :]

        face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5 = face_5[field_set_index, :, :]

    if dim==3:
        depth_level = 10
        Nr = 90

        diag_grid = np.fromfile(file_path, '>f4')
        N = int(np.size(diag_grid) / (102600 * 40 * Nr))
        diag_grid = np.reshape(diag_grid, (N * Nr, 102600, 40))

        if face==1:
            face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
            face_1 = face_1.reshape((N* Nr, 680, 1080))
            face_1 = face_1[field_set_index*Nr:(field_set_index+1)*Nr, :, :]

        if face == 2:
            face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
            face_2 = face_2.reshape((N * Nr, 680, 1080))
            face_2 = face_2[field_set_index*Nr:(field_set_index+1)*Nr, :, :]

        if face == 3:
            face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
            face_3 = face_3.reshape((N * Nr, 1080, 1080))
            face_3 = face_3[field_set_index*Nr:(field_set_index+1)*Nr, :, :]

        if face == 3:
            face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
            face_4 = face_4.reshape((N * Nr, 1080, 680))
            face_4 = face_4[field_set_index*Nr:(field_set_index+1)*Nr, :, :]

        if face == 5:
            face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
            face_5 = face_5.reshape((N * Nr, 1080, 680))
            face_5 = face_5[field_set_index*Nr:(field_set_index+1)*Nr, :, :]
            field_subset = face_5[:,min_row:max_row,min_col:max_col]
            field_subset = np.rot90(field_subset,axes=(1,2))

    # field_subset = np.ma.masked_where(field_subset == 0, field_subset)
    # C = plt.imshow(field_subset[0,:,:])
    # plt.colorbar(C)
    # plt.show()


    return (field_subset)

def save_field_to_nc(output_dir, output_file, XC, YC, depth, Angle_CS, Angle_SN, field_name, field_subset, is_vector):

    depth = depth[:45]
    field_subset = field_subset[:45,:,:]

    output_path = os.path.join(output_dir, output_file)
    ds = nc4.Dataset(output_path,'w')

    ds.createDimension('depth', np.shape(field_subset)[0])
    ds.createDimension('rows', np.shape(field_subset)[1])
    ds.createDimension('cols', np.shape(field_subset)[2])

    dvar = ds.createVariable('depth','f4',('depth'))
    dvar[:] = depth

    xvar = ds.createVariable('longitude','f4',('rows','cols'))
    xvar[:,:] = XC

    yvar = ds.createVariable('latitude', 'f4', ('rows', 'cols'))
    yvar[:, :] = YC

    if is_vector:
        xvar = ds.createVariable('angle_cs', 'f4', ('rows', 'cols'))
        xvar[:, :] = Angle_CS

        yvar = ds.createVariable('angle_sn', 'f4', ('rows', 'cols'))
        yvar[:, :] = Angle_SN

    var = ds.createVariable(field_name,'f4',('depth','rows','cols'))
    var[:,:,:] = field_subset

    ds.close()

def save_disko_fields(config_dir, field_number, remove_old, overwrite):

    field_set, field_name, field_set_index, is_vector = field_number_to_set_and_name(field_number)

    # time.sleep(field_number)
    if 'disko_bay' not in os.listdir(os.path.join(config_dir, 'N1_1080','results')):
        os.mkdir(os.path.join(config_dir, 'N1_1080','results','disko_bay'))
    if field_name not in os.listdir(os.path.join(config_dir, 'N1_1080','results','disko_bay')):
        os.mkdir(os.path.join(config_dir, 'N1_1080','results','disko_bay', field_name))
    output_dir = os.path.join(config_dir, 'N1_1080','results','disko_bay', field_name)

    file_names, dates = get_filename_date_pairs(config_dir, field_set)

    XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, depth = read_faces_geometry(config_dir, is_vector)

    # find the disko indices
    face, min_row,max_row,min_col,max_col = subset_geometry_to_Disko(XC_faces, YC_faces)

    # subset the geometry
    XC = XC_faces[face][min_row:max_row,min_col:max_col]
    YC = YC_faces[face][min_row:max_row, min_col:max_col]
    XC = np.rot90(XC)
    YC = np.rot90(YC)
    if is_vector:
        Angle_CS = AngleCS_faces[face][min_row:max_row, min_col:max_col]
        Angle_SN = AngleSN_faces[face][min_row:max_row, min_col:max_col]
        Angle_CS = np.rot90(Angle_CS)
        Angle_SN = np.rot90(Angle_SN)
    else:
        Angle_CS = {}
        Angle_SN = {}

    for f in range(len(file_names)):
        file_name = file_names[f]
        date = dates[f]
        print('        - Working on nc file for '+field_name+' on '+date+' ('+file_name+')')

        output_file = 'N1_1080_' + field_name + '_' + date + '.nc'

        if overwrite:
            continue_to_processing = True
        else:
            if output_file in os.listdir(output_dir):
                continue_to_processing = False
            else:
                continue_to_processing=True

        if continue_to_processing:

            field_subset = read_field_to_disko_subset(config_dir, field_set, file_name, field_set_index,
                                               face, min_row,max_row,min_col,max_col)

            save_field_to_nc(output_dir, output_file, XC, YC, depth, Angle_CS, Angle_SN, field_name, field_subset, is_vector)



        #     print('    - Creating plot for ' + output_file)
        #     create_plot(output_dir, output_file, XC_faces, YC_faces, field_set, field_name, date, field_faces, dim)
        #
        #     del field_faces




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config_dir", action="store",
                        help="The directory where the N1, N2, and N3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_number", action="store",
                        help="The number of the field to plot (1-4, used for parallel saving on pleiades).", dest="field_number",
                        type=int, required=True)

    parser.add_argument("-r", "--remove_old", action="store",
                        help="Choose whether to remove all of the old files (1 is true, 0 is false).", dest="remove_old",
                        type=int, required=False, default = 0)

    parser.add_argument("-o", "--overwrite", action="store",
                        help="Choose whether to overwrite old files (1 is true, 0 is false).", dest="overwrite",
                        type=int, required=False, default=0)

    args = parser.parse_args()
    config_dir = args.config_dir
    field_number = args.field_number
    remove_old = args.remove_old
    overwrite = args.overwrite

    if remove_old==0:
        remove_old = False
    else:
        remove_old = True

    if overwrite==0:
        overwrite = False
    else:
        overwrite = True

    save_disko_fields(config_dir, field_number, remove_old, overwrite)
   

