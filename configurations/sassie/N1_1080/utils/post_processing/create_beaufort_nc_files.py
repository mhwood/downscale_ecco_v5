
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

def create_source_file_dictionary(config_dir):

    print('Creating the beaufort dv reference dictionary')

    output_period = 3600
    timestep = 120
    iters_per_output = int(output_period/timestep)
    test_var_name = 'AREA'
    Nr = 1
    N = 207104

    dv_dir = os.path.join(config_dir,'N1_1080','run','dv')

    source_dict = {}

    for file_name in os.listdir(dv_dir):
        if 'beaufort_mask_'+test_var_name in file_name:
            file_path = os.path.join(config_dir, 'N1_1080', 'run', 'dv',file_name)
            grid = np.fromfile(file_path, '>f4')
            n_iters_in_this_file = int(np.size(grid) / (Nr * N))
            file_iter = int(file_name.split('.')[-2])

            first_iter = (file_iter-1)+int(iters_per_output/2)
            last_iter = (file_iter-1)+int(iters_per_output/2) + (n_iters_in_this_file-1)*iters_per_output

            first_date = datetime(1992,1,1) + timedelta(seconds=first_iter*timestep)
            last_date = datetime(1992, 1, 1) + timedelta(seconds=last_iter * timestep)

            source_dict[file_iter] = [first_iter,last_iter,first_date,last_date]

    output = 'iter_number\tfirst_iter\tlast_iter\tfirst_date\tlast_date'

    source_iters = sorted(list(source_dict.keys()))
    for iter_number in source_iters:
        output+='\n'+str(iter_number)+\
                '\t'+str(source_dict[iter_number][0])+\
                '\t'+str(source_dict[iter_number][1])+\
                '\t'+str(source_dict[iter_number][2])+\
                '\t'+str(source_dict[iter_number][3])

    output_file = os.path.join(config_dir,'N1_1080','run','dv','beaufort_mask_source.dict')
    f = open(output_file,'w')
    f.write(output)
    f.close()

def read_location_dictionary(config_dir):

    file_path = os.path.join(config_dir,'N1_1080','input','N1_1080_dv_mask_ref.nc')

    ds = nc4.Dataset(file_path)
    faces = ds.groups['beaufort'].variables['source_faces'][:]
    rows = ds.groups['beaufort'].variables['source_rows'][:]
    cols = ds.groups['beaufort'].variables['source_cols'][:]
    ds.close()

    return(faces, rows, cols)

def read_subset_geometry(config_dir, faces, rows, cols):

    # determine the extents within the faces
    # face 3
    face_3_min_row = np.min(rows[faces == 3])
    face_3_max_row = np.max(rows[faces == 3])
    face_3_min_col = np.min(cols[faces == 3])
    face_3_max_col = np.max(cols[faces == 3])
    # print('geometry',face_3_min_row,face_3_max_row,face_3_min_col,face_3_max_col)

    # face 4
    face_4_min_row = np.min(rows[faces == 4])
    face_4_max_row = np.max(rows[faces == 4])
    face_4_min_col = np.min(cols[faces == 4])
    face_4_max_col = np.max(cols[faces == 4])
    # print('geometry',face_4_min_row, face_4_max_row, face_4_min_col, face_4_max_col)

    min_row = np.min([face_3_min_row, face_4_min_row])
    max_row = np.max([face_3_max_row, face_4_max_row])

    # get the geometry
    face_3_grid = np.fromfile(os.path.join(config_dir,'N1_1080','input','tile003.mitgrid'),'>f8')
    face_3_grid = np.reshape(face_3_grid,(16,1081,1081))
    face_3_XC = face_3_grid[0, min_row:max_row + 1, face_3_min_col:-1]
    face_3_YC = face_3_grid[1, min_row:max_row + 1, face_3_min_col:-1]
    # print(np.shape(face_3_XC))

    face_4_grid = np.fromfile(os.path.join(config_dir, 'N1_1080', 'input', 'tile004.mitgrid'), '>f8')
    face_4_grid = np.reshape(face_4_grid, (16, 1081, 681))
    face_4_XC = face_4_grid[0, min_row:max_row + 1, :face_4_max_col+1]
    face_4_YC = face_4_grid[1, min_row:max_row + 1, :face_4_max_col+1]
    # print(np.shape(face_4_XC))

    XC = np.concatenate([face_3_XC, face_4_XC], axis=1)
    XC[XC<0] += 360
    YC = np.concatenate([face_3_YC, face_4_YC], axis=1)

    XC = np.rot90(XC)
    YC = np.rot90(YC)

    # plt.subplot(1,2,1)
    # C = plt.imshow(XC,origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(YC, origin='lower')
    # plt.colorbar(C)
    #
    # plt.show()

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

def read_source_iter_dict(config_dir):

    file_name = os.path.join(config_dir,'N1_1080','run','dv','beaufort_mask_source.dict')
    f = open(file_name)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)
    source_dict = {}
    for line in lines:
        line = line.split('\t')
        if len(line)>2:
            source_dict[int(line[0])] = [int(line[1]),int(line[2])]

    return(source_dict)

def get_source_iter_numbers(year,month,day, source_iter_dict):

    output_period = 3600
    timestep = 120
    iters_per_output = int(output_period / timestep)

    iter_centerpoint_datetimes = []
    for hour in range(24):
        iter_centerpoint_datetimes.append(datetime(year,month,day,hour,minute=30))

    iter_numbers = []
    for date in iter_centerpoint_datetimes:
        seconds_elapsed = (date - datetime(1992, 1, 1)).total_seconds()
        iter_number = seconds_elapsed / timestep
        iter_numbers.append(iter_number)
    # iter_numbers = np.array(iter_numbers)

    iters_found = np.zeros((24,))
    iter_sources = np.zeros((24, 2))
    unique_sources = []

    for file_iter in list(source_iter_dict.keys()):
        if np.sum(iters_found) != 24:
            # print('        - Looking for iters in file with iter ' + str(file_iter))
            file_iterations = np.arange(source_iter_dict[file_iter][0], source_iter_dict[file_iter][1] + 1,
                                        iters_per_output)
            for iter_number in file_iterations:
                if iter_number in iter_numbers:
                    ind = iter_numbers.index(iter_number)
                    iters_found[ind] = 1
                    iter_sources[ind,0] = file_iter
                    iter_sources[ind,1] = np.argmin(np.abs(file_iterations-iter_number))
                    if file_iter not in unique_sources:
                        unique_sources.append(file_iter)
    iters_found = np.sum(iters_found)

    return(unique_sources, iter_sources, iters_found)

def dv_field_name_to_diag_field_name(var_name):

    has_Nr = False

    if var_name=='AQH':
        diag_var_name = 'EXFaqh'
        description = 'surface (2-m) specific humidity'
        units = 'kg/kg'
    elif var_name=='AREA':
        diag_var_name = 'SIarea'
        description = 'SEAICE fractional ice-covered area [0 to 1]'
        units = 'm^2/m^2'
    elif var_name=='ATEMP':
        diag_var_name = 'EXFatemp'
        description = 'surface (2-m) air temperature'
        units = 'degK'
    elif var_name=='EMPMR':
        diag_var_name = 'EXFempmr'
        description = 'net upward freshwater flux, > 0 increases salinity'
        units = 'm/s'
    elif var_name=='ETAN':
        diag_var_name = 'ETAN'
        description = 'Surface Height Anomaly'
        units = 'm'
    elif var_name=='EVAP':
        diag_var_name = 'EXFevap'
        description = 'evaporation, > 0 increases salinity'
        units = 'm/s'
    elif var_name=='HEFF':
        diag_var_name = 'SIheff'
        description = 'SEAICE effective ice thickness'
        units = 'm'
    elif var_name=='HL':
        diag_var_name = 'EXFhl'
        description = 'Latent heat flux into ocean, >0 increases theta'
        units = 'W/m^2'
    elif var_name=='HS':
        diag_var_name = 'EXFhs'
        description = 'Sensible heat flux into ocean, >0 increases theta'
        units = 'W/m^2'
    elif var_name=='HSNOW':
        diag_var_name = 'SIhsnow'
        description = 'SEAICE effective snow thickness'
        units = 'm'
    elif var_name=='KPPDIFFS':
        diag_var_name = 'KPPdiffS'
        description = 'Vertical diffusion coefficient for salt & tracers'
        units = 'm^2/s'
        has_Nr = True
    elif var_name=='KPPDIFFT':
        diag_var_name = 'KPPdiffT'
        description = 'Vertical diffusion coefficient for heat'
        units = 'm^2/s'
        has_Nr = True
    elif var_name=='KPPFRAC':
        diag_var_name = 'KPPfrac'
        description = 'Short-wave flux fraction penetrating mixing layer'
        units = ''
    elif var_name=='KPPHBL':
        diag_var_name = 'KPPhbl'
        description = 'KPP boundary layer depth, bulk Ri criterion'
        units = 'm'
    elif var_name=='KPPVISCA':
        diag_var_name = 'KPPviscA'
        description = 'KPP vertical eddy viscosity coefficient'
        units = 'm^2/s'
        has_Nr = True
    elif var_name=='LWDOWN':
        diag_var_name = 'EXFlwdn'
        description = 'Downward longwave radiation, >0 increases theta'
        units = 'W/m^2'
    elif var_name=='LWFLUX':
        diag_var_name = 'EXFlwnet'
        description = 'Net upward longwave radiation, >0 decreases theta'
        units = 'W/m^2'
    elif var_name=='PRECIP':
        diag_var_name = 'EXFpreci'
        description = 'recipitation, > 0 decreases salinity'
        units = 'm/s'
    elif var_name=='QNET':
        diag_var_name = 'EXFqnet'
        description = 'Net upward heat flux (turb+rad), >0 decreases theta'
        units = 'W/m^2'
    elif var_name=='RUNOFF':
        diag_var_name = 'EXFroff'
        description = 'river runoff, > 0 decreases salinity'
        units = 'm/s'
    elif var_name=='SALT':
        diag_var_name = 'SALT'
        description = 'Salinity'
        units = 'psu'
        has_Nr = True
    elif var_name=='SICELOAD':
        diag_var_name = 'sIceLoad'
        description = 'sea-ice loading (in Mass of ice+snow / area unit)'
        units = 'kg/m^2'
    elif var_name=='SWDOWN':
        diag_var_name = 'EXFswdn'
        description = 'Downward shortwave radiation, >0 increases theta'
        units = 'W/m^2'
    elif var_name=='SWFLUX':
        diag_var_name = 'EXFswnet'
        description = 'Net upward shortwave radiation, >0 decreases theta'
        units = 'W/m^2'
    elif var_name=='THETA':
        diag_var_name = 'THETA'
        description = 'Potential Temperature'
        units = 'degC'
        has_Nr = True
    elif var_name=='UICE':
        diag_var_name = 'SIuice'
        description = 'SEAICE zonal ice velocity, >0 from West to East'
        units = 'm/s'
    elif var_name=='USTRESS':
        diag_var_name = 'EXFtaux'
        description = 'zonal surface wind stress, >0 increases uVel'
        units = 'N/m^2'
    elif var_name=='UVEL':
        diag_var_name = 'UVEL'
        description = 'Zonal Component of Velocity (m/s)'
        units = 'm/s'
        has_Nr = True
    elif var_name=='UWIND':
        diag_var_name = 'EXFuwind'
        description = 'zonal 10-m wind speed, >0 increases uVel'
        units = 'm/s'
    elif var_name=='VICE':
        diag_var_name = 'SIvice'
        description = 'SEAICE merid. ice velocity, >0 from South to North'
        units = 'm/s'
    elif var_name=='VSTRESS':
        diag_var_name = 'EXFtauy'
        description = 'meridional surface wind stress, >0 increases vVel'
        units = 'N/m^2'
    elif var_name=='VVEL':
        diag_var_name = 'VVEL'
        description = 'Zonal Component of Velocity (m/s)'
        units = 'm/s'
        has_Nr = True
    elif var_name=='VWIND':
        diag_var_name = 'EXFvwind'
        description = 'meridional 10-m wind speed, >0 increases vVel'
        units = 'm/s'
    elif var_name=='WVEL':
        diag_var_name = 'WVEL'
        description = 'Vertical Component of Velocity (m/s)'
        units = 'm/s'
        has_Nr = True
    else:
        raise ValueError('Havent defined this field yet')

    return(diag_var_name,description,units,has_Nr)

def read_field_from_dv_files(config_dir, var_name, has_Nr, faces, rows, cols, unique_sources, iter_sources):

    if has_Nr:
        Nr = 90
    else:
        Nr = 1
    grid_faces = {}
    grid_faces[3] = np.zeros((24, Nr, 1080, 1080))
    grid_faces[4] = np.zeros((24, Nr, 1080, 680))

    for file_iter in unique_sources:

        print('        - Reading in file '+'beaufort_mask_'+var_name+'.'+'{:010d}'.format(file_iter)+'.bin')
        file_path = os.path.join(config_dir,'N1_1080','run','dv',
                                 'beaufort_mask_'+var_name+'.'+'{:010d}'.format(file_iter)+'.bin')
        grid = np.fromfile(file_path,'>f4')
        ntimesteps = int(np.size(grid) / (Nr * len(faces)))
        grid = np.reshape(grid, (ntimesteps, Nr, len(faces)))

        print(np.sum(grid))

        min_src_ind = int(np.min(iter_sources[iter_sources[:, 0] == file_iter, 1]))
        max_src_ind = int(np.max(iter_sources[iter_sources[:, 0] == file_iter, 1]))
        # print('src',min_src_ind,max_src_ind)
        print('        - Adding indices '+str(min_src_ind)+' to '+str(max_src_ind)+' from '+'beaufort_mask_'+var_name+'.'+'{:010d}'.format(file_iter)+'.bin'+' to output grid')

        min_dst_ind = np.argmin(np.abs(iter_sources[:, 1] - min_src_ind))
        max_dst_ind = np.argmin(np.abs(iter_sources[:, 1] - max_src_ind))
        # print('dst', min_dst_ind, max_dst_ind)

        for j in range(len(faces)):
            if has_Nr:
                if j%10000==0:
                    print('            - Completed points '+str(j)+' out of '+str(len(faces)))
            grid_faces[faces[j]][min_dst_ind:max_dst_ind+1,:,rows[j],cols[j]] = grid[min_src_ind:max_src_ind+1,:,j]

    if has_Nr:
        print('         - Concatenating the grids along axis 3')
    var_grid = np.concatenate([grid_faces[3],grid_faces[4]],axis=3)


    #######################################################
    # this section will output the grid even if its all 0's

    min_row = 480
    max_row = 999
    min_col = 681
    max_col = 1160

    var_grid = var_grid[:, :, min_row:max_row + 1, min_col:max_col + 1]
    var_grid = np.rot90(var_grid, axes=(2, 3))

    do_output = True

    # #######################################################
    # # this section won't output the grid if its all 0's
    #
    # if np.any(var_grid[0,0,:,:]!=0):
    #
    #     # rs, cs = np.where(var_grid[0,0,:,:]!=0)
    #     # min_row = np.min(rs)
    #     # max_row = np.max(rs)
    #     # min_col = np.min(cs)
    #     # max_col = np.max(cs)
    #
    #     # C = plt.imshow(var_grid[0,0,:,:],origin='lower')
    #     # plt.colorbar(C)
    #     # plt.show()
    #
    #     var_grid = var_grid[:, :, min_row:max_row+1,min_col:max_col+1]
    #     var_grid = np.rot90(var_grid,axes=(2,3))
    #
    #     # plt.imshow(var_grid[0,0,:,:],origin='lower')
    #     # plt.show()
    #
    #     do_output = True
    # else:
    #     do_output = False

    return(var_grid, do_output)

def output_file_to_nc(config_dir, diag_var_name, description, units, has_Nr, year, month, day, XC, YC, Z, var_grid):

    output_file = os.path.join(config_dir, 'N1_1080', 'results', 'beaufort', diag_var_name,
                               'beaufort_'+diag_var_name+'_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time', np.shape(var_grid)[0])
    ds.createDimension('cols',np.shape(XC)[1])
    ds.createDimension('rows',np.shape(YC)[0])

    xvar = ds.createVariable('Longitude','f4',('rows','cols'))
    xvar[:,:] = XC
    yvar = ds.createVariable('Latitude', 'f4', ('rows', 'cols'))
    yvar[:, :] = YC

    if has_Nr:
        Z_out = Z[:82]
        var_grid = var_grid[:,:82,:,:]

        ds.createDimension('depth', np.shape(Z_out)[0])
        zvar = ds.createVariable('Depth', 'f4', ('depth',))
        zvar[:] = Z_out
        var = ds.createVariable(diag_var_name,'f4',('time','depth','rows','cols'))
        var[:,:,:,:] = var_grid
    else:
        var_grid = var_grid[:,0,:,:]
        var = ds.createVariable(diag_var_name, 'f4', ('time', 'rows', 'cols'))
        var[:, :, :] = var_grid

    var.description = description
    var.units = units

    ds.close()




########################################################################################################################

def create_beaufort_dataset(config_dir, field_number, remove_old, overwrite):

    print('Creating the beaufort dataset for field number '+str(field_number))

    if field_number==1:
        print('    - Creating the source iteration dict')
        if 'beaufort_mask_source.dict' not in os.listdir(os.path.join(config_dir,'N1_1080','run','dv')):
            create_source_file_dictionary(config_dir)
    else:
        # if proc 1 is still working on the dictionary, wait for it to finish
        for attempt in range(100):
            if 'beaufort_mask_source.dict' not in os.listdir(os.path.join(config_dir, 'N1_1080', 'run', 'dv')):
                time.sleep(10)

    print('    - Reading in the reference dictionary')
    faces, rows, cols = read_location_dictionary(config_dir)

    print('    - Reading in the geometry')
    XC, YC, Z = read_subset_geometry(config_dir, faces, rows, cols)
    print('        - Domain shape: '+str(np.shape(XC)))

    print('    - Reading in the source iteration dict')
    source_iter_dict = read_source_iter_dict(config_dir)

    var_names_2D = ['AREA','EMPMR','ETAN','HEFF','HSNOW',
                    'QNET','SICELOAD',
                    'UICE','USTRESS','UWIND','VICE','VSTRESS','VWIND',
                    'AQH', 'ATEMP', 'EVAP', 'HL', 'HS', 'LWDOWN','LWFLUX','SWFLUX', 'SWDOWN', 'PRECIP','RUNOFF',
                    'KPPFRAC', 'KPPHBL']

    # something is wrong with this one?
    # 'SWFLUX',

    var_names_3D = ['SALT','THETA','UVEL','VVEL','WVEL','KPPDIFFS', 'KPPDIFFT', 'KPPVISCA']

    var_names = var_names_2D+var_names_3D

    var_name = var_names[field_number-1]

    diag_var_name, description, units, has_Nr = dv_field_name_to_diag_field_name(var_name)

    if 'results' not in os.listdir(os.path.join(config_dir, 'N1_1080')):
        os.mkdir(os.path.join(config_dir, 'N1_1080', 'results'))
    if 'beaufort' not in os.listdir(os.path.join(config_dir, 'N1_1080', 'results')):
        os.mkdir(os.path.join(config_dir, 'N1_1080', 'results', 'beaufort'))
    if diag_var_name not in os.listdir(os.path.join(config_dir, 'N1_1080', 'results', 'beaufort')):
        os.mkdir(os.path.join(config_dir, 'N1_1080', 'results', 'beaufort', diag_var_name))

    print('    - Creating the daily nc files for variable '+var_name)
    print('            - Output name: '+diag_var_name)
    print('            - Description: '+description)
    print('            - Units: '+units)
    print('            - Has Depth: '+str(has_Nr))

    dates = [datetime(2014,6,7), datetime(2014,6,8)]

    for date in dates:
        year = date.year
        month = date.month
        day = date.day

        output_file = 'beaufort_'+diag_var_name+'_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.nc'
        if output_file not in os.listdir(os.path.join(config_dir, 'N1_1080', 'results', 'beaufort', diag_var_name)):

            print('        - Working on date '+str(year)+'/'+str(month)+'/'+str(day))

            # unique sources is a list of file iters
            # iter_sources is a 2d list where each row is the file iter and the index within the file
            unique_sources, iter_sources, iters_found = get_source_iter_numbers(year, month, day, source_iter_dict)

            var_grid, do_output = read_field_from_dv_files(config_dir, var_name, has_Nr, faces, rows, cols, unique_sources, iter_sources)

            if do_output:
                print('            - The output grid has shape '+str(np.shape(var_grid)))
                output_file_to_nc(config_dir, diag_var_name, description, units, has_Nr, year, month, day, XC, YC, Z, var_grid)
        else:
            print('        - Skipping date '+str(year)+'/'+str(month)+'/'+str(day)+' - already created')



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

    create_beaufort_dataset(config_dir, field_number, remove_old, overwrite)
