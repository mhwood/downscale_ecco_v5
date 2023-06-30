
import os
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import time
import argparse
import ast
import sys
from datetime import datetime, timedelta


def var_name_to_metadata(var_name):
    meta_dict = {'THETA':['ocean_state_3D_day_mean',2,1,'Potential Temperature','degrees C',True],
                 'SALT':['ocean_state_3D_day_mean',2,0,'Practical Salinity','g/kg',True]}
    return(meta_dict[var_name])


def read_subset_geometry(config_dir):

    min_col = 600
    max_col = 900
    max_row = 850
    min_row = 625

    # get the geometry
    face_3_grid = np.fromfile(os.path.join(config_dir,'N1_1080','input','tile003.mitgrid'),'>f8')
    face_3_grid = np.reshape(face_3_grid,(16,1081,1081))
    face_3_XC = face_3_grid[0, :-1 ,:-1]
    face_3_YC = face_3_grid[1, :-1 ,:-1]

    face_5_grid = np.fromfile(os.path.join(config_dir, 'N1_1080', 'input', 'tile005.mitgrid'), '>f8')
    face_5_grid = np.reshape(face_5_grid, (16, 1081, 681))
    face_5_XC = face_5_grid[0,  :-1 ,:-1]
    face_5_YC = face_5_grid[1,  :-1 ,:-1]

    XC = np.concatenate([np.rot90(face_5_XC,k=1),
                         np.rot90(face_3_XC,k=2)], axis=0)
    YC = np.concatenate([np.rot90(face_5_YC,k=1),
                         np.rot90(face_3_YC,k=2)], axis=0)

    XC = XC[min_row:max_row, min_col:max_col]
    YC = YC[min_row:max_row, min_col:max_col]

    # plt.subplot(1,2,1)
    # C = plt.imshow(XC,origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(YC, origin='lower')
    # plt.colorbar(C)

    plt.show()

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

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z_out = (Z_bottom + Z_top) / 2

    return(XC,YC,Z_out)


def generate_month_file_iters(year, month):

    # calculate the number of days in the month
    if month in [1,3,5,7,8,10,12]:
        n_days = 31
    elif month in [4,6,9,11]:
        n_days = 30
    else:
        if year%4==0:
            n_days = 29
        else:
            n_days = 28
    days = np.arange(1,n_days+1)

    seconds_per_iter = 120

    iter_numbers = []
    for day in days:
        date = datetime(year, month, day)
        total_seconds = (date - datetime(1992, 1, 1)).total_seconds()
        iter_number = total_seconds / seconds_per_iter
        iter_numbers.append(int(iter_number))

    return(days,iter_numbers)


def check_all_iters_present(config_dir,diag_subset,iter_numbers):
    diag_dir = os.path.join(config_dir, 'N1_1080', 'run', 'diags', diag_subset)

    all_iters_found = True
    total_count = 0
    for iter_number in iter_numbers:
        file_name = diag_subset + '.' + '{:010d}'.format(iter_number) + '.data'
        if file_name not in os.listdir(diag_dir):
            all_iters_found = False
            print('         - Missing '+file_name)
        else:
            total_count +=1

    print('        - Found '+str(total_count)+' of '+str(len(iter_numbers))+' files')

    return(all_iters_found)


def read_diag_to_W_Greenland_grid(sf, file_path, n_field_levels, field_level, Nr):

    llc = 1080
    n = 680
    total_points = (n*4*llc + llc*llc)*Nr


    sassie_n1_compact = np.fromfile(file_path, '>f4')
    sassie_n1_compact = sassie_n1_compact[field_level*total_points:(field_level+1)*total_points]
    sassie_n1_compact = np.reshape(sassie_n1_compact,(Nr,4*n+llc,llc))
    var_faces = sf.sassie_n1_compact_to_faces(sassie_n1_compact, llc, n)

    W_Greenland_subset = np.concatenate([np.rot90(var_faces[5],axes=(1,2),k=1),
                                         np.rot90(var_faces[3],axes=(1,2),k=2)], axis=1)

    min_x_index = 600
    max_x_index = 900
    max_y_index = 850
    min_y_index = 625

    W_Greenland_subset = W_Greenland_subset[:,min_y_index:max_y_index,min_x_index:max_x_index]

    # plt.subplot(1,3,1)
    # plt.pcolormesh(var_faces[3][0,:,:])
    #
    # plt.subplot(1, 3, 2)
    # plt.pcolormesh(var_faces[5][0,:,:])
    #
    # plt.subplot(1, 3, 3)
    # plt.pcolormesh(W_Greenland_subset[0,:,:])
    #
    # plt.show()

    return(W_Greenland_subset)


def output_file_to_nc(config_dir, var_name, description, units, has_Nr, year, month, days, iter_numbers, XC, YC, Z, var_grid):

    if var_name not in os.listdir(os.path.join(config_dir, 'N1_1080', 'results', 'W_Greenland')):
        os.mkdir(os.path.join(config_dir, 'N1_1080', 'results', 'W_Greenland',var_name))

    output_file = os.path.join(config_dir, 'N1_1080', 'results', 'W_Greenland', var_name,
                               'W_Greenland_'+var_name+'_'+str(year)+'{:02d}'.format(month)+'.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time', np.shape(var_grid)[0])
    ds.createDimension('cols',np.shape(XC)[1])
    ds.createDimension('rows',np.shape(YC)[0])

    xvar = ds.createVariable('longitude','f4',('rows','cols'))
    xvar[:,:] = XC
    xvar.description = 'longitude of the cell center'

    yvar = ds.createVariable('latitude', 'f4', ('rows', 'cols'))
    yvar[:, :] = YC
    yvar.description = 'latitude of the cell center'

    ivar = ds.createVariable('iter_numbers', 'i8', ('time',))
    ivar[:] = iter_numbers
    ivar.description = 'iteration number of the model run referenced to 1992/1/1 at timesteps of 120 seconds'

    yrvar = ds.createVariable('years', 'i4', ('time',))
    yrvar[:] = year*np.ones_like(days)
    yrvar.description = 'nominal year of the output field'

    mvar = ds.createVariable('months', 'i4', ('time',))
    mvar[:] = month*np.ones_like(days)
    mvar.description = 'nominal month of the output field'

    dvar = ds.createVariable('days', 'i4', ('time',))
    dvar[:] = days
    dvar.description = 'nominal day of year of the output field'

    if has_Nr:
        Z_out = Z[:72]
        var_grid = var_grid[:,:72,:,:]

        ds.createDimension('depth', np.shape(Z_out)[0])
        zvar = ds.createVariable('depth', 'f4', ('depth',))
        zvar[:] = Z_out
        var = ds.createVariable(var_name,'f4',('time','depth','rows','cols'))
        var[:,:,:,:] = var_grid

    else:
        var_grid = var_grid[:,0,:,:]
        var = ds.createVariable(var_name, 'f4', ('time', 'rows', 'cols'))
        var[:, :, :] = var_grid

    var.description = description
    var.units = units

    ds.close()




########################################################################################################################

def create_W_Greenland_dataset(config_dir, field_number, remove_old, overwrite):

    sys.path.insert(1, os.path.join('..', '..', 'utils'))
    import sassie_functions as sf

    var_names_2D = ['AREA', 'EMPMR', 'ETAN', 'HEFF', 'HSNOW',
                    'QNET', 'SICELOAD',
                    'UICE', 'USTRESS', 'UWIND', 'VICE', 'VSTRESS', 'VWIND',
                    'AQH', 'ATEMP', 'EVAP', 'HL', 'HS', 'LWDOWN', 'LWFLUX', 'SWFLUX', 'SWDOWN', 'PRECIP', 'RUNOFF',
                    'KPPFRAC', 'KPPHBL']
    var_names_3D = ['SALT', 'THETA', 'UVEL', 'VVEL', 'WVEL', 'KPPDIFFS', 'KPPDIFFT', 'KPPVISCA']
    var_names = var_names_2D + var_names_3D
    var_name = var_names[field_number - 1]

    print('Creating the West Greenland dataset for '+var_name+' (field number '+str(field_number)+')')

    print('    - Reading in the variable metadata')
    metadata = var_name_to_metadata(var_name)
    diag_subset = metadata[0]
    n_field_levels = metadata[1]
    field_level = metadata[2]
    description = metadata[3]
    units = metadata[4]
    has_Nr = metadata[5]

    print('    - Reading in the W Greenland geometry')
    XC, YC, ZC = read_subset_geometry(config_dir)
    Nr = len(ZC)
    # Nr = 90

    for year in range(2015,2017):
        for month in range(1,13):

            print('    - Reading in the W Greenland subset in '+str(year)+'/'+str(month))
            days, iter_numbers = generate_month_file_iters(year, month)
            print('         - min iter: ' + str(min(iter_numbers)))
            print('         - max iter: ' + str(max(iter_numbers)))

            # print('    - Checking that all iter numbers are present')
            # all_iters_found = check_all_iters_present(config_dir,diag_subset,iter_numbers)
            #
            #
            # all_iters_found = False
            #
            # # days = days[:2]
            # # iter_numbers = iter_numbers[:2]
            # if all_iters_found:
            #
            #     for i in range(len(iter_numbers)):
            #         iter_number = iter_numbers[i]
            #         print('        - Working on iter number '+str(iter_number))
            #
            #         file_path = os.path.join(config_dir, 'N1_1080', 'run', 'diags', diag_subset,
            #                                  diag_subset+'.'+'{:010d}'.format(iter_number)+'.data')
            #         W_Greenland_subset = read_diag_to_W_Greenland_grid(sf, file_path, n_field_levels, field_level, Nr)
            #
            #         if i==0:
            #             output_subset = np.reshape(W_Greenland_subset,(1,np.shape(W_Greenland_subset)[0],
            #                                                            np.shape(W_Greenland_subset)[1],
            #                                                            np.shape(W_Greenland_subset)[2]))
            #         else:
            #             W_Greenland_subset = np.reshape(W_Greenland_subset, (1, np.shape(W_Greenland_subset)[0],
            #                                                                  np.shape(W_Greenland_subset)[1],
            #                                                                  np.shape(W_Greenland_subset)[2]))
            #             output_subset = np.concatenate([output_subset,W_Greenland_subset],axis=0)
            #
            #     output_file_to_nc(config_dir, var_name, description, units, has_Nr,
            #                       year, month, days, iter_numbers, XC, YC, ZC, output_subset)
            #
            # else:
            #     print('    - Missing iters - '+str(year)+'/'+str(month)+' skipped')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config_dir", action="store",
                        help="The directory where the N1, N2, and N3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_number", action="store",
                        help="The number of the field to plot (1-24, used for parallel plotting on pleiades).", dest="field_number",
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

    create_W_Greenland_dataset(config_dir, field_number, remove_old, overwrite)
