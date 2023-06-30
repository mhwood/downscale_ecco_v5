import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
from MITgcmutils import mds
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import xesmf as xe
import xarray as xr
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf

def read_pickup_file_to_faces(input_init_dir,pickup_file='pickup.0000000001'):

    Nr = 50
    print('      Reading from '+input_init_dir + '/' + pickup_file)

    global_data, _, global_metadata = mds.rdmds(input_init_dir + '/' + pickup_file,
                                                returnmeta=True)
    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    all_var_grid_faces = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grid_faces = llc_compact_to_faces(var_grid, less_output=True)

        all_var_grid_faces.append(var_grid_faces)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,all_var_grid_faces,global_metadata)


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

def read_regridder_functions_from_nc(L0_XC_faces, L0_YC_faces, L1_XC_faces, L1_YC_faces):

    regridding_faces = {}

    for face in range(1,6):
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

        filename = os.path.join('..', 'input', 'face_'+str(face)+'_downscale_weights.nc')
        regridder = xe.Regridder(ds_in, ds_out, method='bilinear', reuse_weights=True, filename = filename)
        regridding_faces[face] = regridder

    return(regridding_faces)

def create_mean_vertical_difference_profile(var_XC, var_YC, var_grid, XC_subset, YC_subset):

    # make a rough subset of the variable to sample points in the local vicinity
    min_x_index = np.argmin(np.abs(var_XC[int(np.shape(var_XC)[0]/2),:] - np.min(XC_subset)))
    max_x_index = np.argmin(np.abs(var_XC[int(np.shape(var_XC)[0]/2),:] - np.max(XC_subset)))
    min_y_index = np.argmin(np.abs(var_YC[:,int(np.shape(var_XC)[1]/2)] - np.min(YC_subset)))
    max_y_index = np.argmin(np.abs(var_YC[:,int(np.shape(var_XC)[1]/2)] - np.max(YC_subset)))
    if max_y_index<min_y_index:
        tmp = np.copy(min_y_index)
        min_y_index = np.copy(max_y_index)
        max_y_index = tmp

    var_grid_subset = var_grid[:,min_y_index:max_y_index,:]
    var_grid_subset = var_grid_subset[:,:,min_x_index:max_x_index]

    # plt.imshow(var_grid_subset[0,:,:])
    # plt.show()

    average_profile = np.zeros((np.shape(var_grid)[0],))
    count_profile = np.zeros((np.shape(var_grid)[0],))
    average_diff_profile = np.zeros((np.shape(var_grid)[0]-1,))
    count_diff_profile = np.zeros((np.shape(var_grid)[0] - 1,))

    for row in range(np.shape(var_grid_subset)[1]):
        for col in range(np.shape(var_grid_subset)[2]):
            profile = var_grid_subset[:,row,col]
            difference_profile = np.zeros_like(profile)[:-1]
            for p in range(len(profile)-1):
                if profile[p+1]!=0 and profile[p]!=0:
                    difference_profile[p] = profile[p+1]-profile[p]

            average_profile[profile != 0] += profile[profile != 0]
            count_profile[profile != 0] += 1
            average_diff_profile[difference_profile != 0] += difference_profile[difference_profile != 0]
            count_diff_profile[difference_profile != 0] += 1

    average_profile[count_profile!=0] = average_profile[count_profile!=0]/count_profile[count_profile!=0]
    average_diff_profile[count_diff_profile != 0] = average_diff_profile[count_diff_profile != 0] / count_diff_profile[count_diff_profile != 0]

    average_diff_profile = np.zeros_like(average_profile)[:-1]
    for p in range(len(average_profile) - 1):
        if average_profile[p + 1] != 0 and average_profile[p] != 0:
            average_diff_profile[p] = average_profile[p + 1] - average_profile[p]

    mean_vertical_difference = np.mean(average_diff_profile[average_diff_profile!=0])

    # plt.subplot(1,2,1)
    # plt.plot(average_profile)
    # plt.subplot(1, 2, 2)
    # plt.plot(average_diff_profile)
    # plt.show()

    return(mean_vertical_difference)

def stack_grids_to_pickup(var_names, interp_grids):
    counter = 0

    print('  - Stacking all of the downscaled grids together...')
    for g in range(len(interp_grids)-1,-1,-1):
        grid = interp_grids[g]
        print('      + Adding in '+str(var_names[g]))

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([grid, pickup_grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    subset_metadata['nrecords'][0] = np.shape(pickup_grid)[0]

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [   "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    if subset_metadata['timeinterval'][0][0]<1.0:
        time_interval_base = 0
        time_interval_exponent = 0
    else:
        time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
        time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [   "+str(subset_metadata['nflds'][0])+" ];\n"
    output += " fldList = {\n "
    for var_name in subset_metadata['fldlist']:
        output += "'"+var_name
        for i in range(8-len(var_name)):
            output+= " "
        output+="' "
    output += "\n };"

    f = open(output_file+'.meta','w')
    f.write(output)
    f.close()


########################################################################################################################

def create_pickup_file(ecco_path):

    LLC = 270
    N = 180

    llc = 1080
    n = 720

    Nr_in = 50
    Nr_out = 90

    print('Creating the pickup file for the N1_'+str(llc)+' model from the output of the N0_'+str(LLC)+' model')

    print('    - Reading in variables from the N0_'+str(LLC)+' pickup file')
    input_init_dir = os.path.join('..','..', 'N0_'+str(LLC), 'input')
    pickup_file = 'pickup.after_adjustments'
    var_names,row_bounds,all_var_grid_faces,global_metadata = read_pickup_file_to_faces(input_init_dir,pickup_file)

    print('    - Reading in the geometry of the N0_'+str(LLC)+' domain')
    N0_XC_faces, N0_YC_faces = sf.read_sassie_n0_grid(ecco_path)

    N0_XC_faces_extended = {}
    N0_YC_faces_extended = {}
    for face in range(1, 6):
        N0_XC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_XC_faces, face)
        N0_YC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_YC_faces, face)

    print('    - Reading in the geometry of the N1_' + str(LLC) + ' domain')
    input_dir = os.path.join('..','input')
    N1_XC_faces, N1_YC_faces = sf.read_sassie_n1_grid(input_dir)

    print('    - Reading in the regridder functions')
    regridding_weights = df.read_regridder_functions_from_nc(input_dir, N0_XC_faces_extended, N0_YC_faces_extended,
                                                             N1_XC_faces, N1_YC_faces)

    N0_delR = sf.read_sassie_delR(domain_level=0)
    N1_delR = sf.read_sassie_delR(domain_level=1)

    print('    - Downscaling the pickup grids')
    interp_grids = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        print('      - Downscaling ' + var_name)
        if var_name not in []:#['Theta','Salt','Vvel','GvNm1','GvNm2','Uvel','GuNm1','GuNm2','EtaN','dEtaHdt','EtaH']:

            if var_name in ['Vvel','GvNm1','GvNm2']:
                N1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'), hFac='S')
                N0_wet_grid_3D_on_N1_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='S')
            elif var_name in ['Uvel','GuNm1','GuNm2']:
                N1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'), hFac='W')
                N0_wet_grid_3D_on_N1_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='W')
            else:
                N1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'), hFac='C')
                N0_wet_grid_3D_on_N1_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='C')

            var_grid_faces = all_var_grid_faces[vn]

            for face in range(1, 6):
                if face < 3:
                    var_grid_faces[face] = var_grid_faces[face][:,-N:, :]
                if face > 3:
                    var_grid_faces[face] = var_grid_faces[face][:,:, :N]

            var_grid_faces_extended = {}
            for face in range(1, 6):
                if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                    var_grid_faces_extended[face] = sf.get_extended_var_grid_on_face_3D(var_grid_faces, face, Nr_in)
                else:
                    var_grid_faces_extended[face] = sf.get_extended_var_grid_on_face_3D(var_grid_faces, face, 1)

            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                var_grid_faces_extended = df.interpolate_var_grid_faces_to_new_depth_levels(var_grid_faces_extended,
                                                                                            N0_wet_grid_3D_on_N1_faces,
                                                                                            N0_delR, N1_delR)

            if var_name.lower() not in ['etan','detahdt','etah']:
                row_col_levels_per_face = {1:(Nr_out, n, llc),
                                           2:(Nr_out, n, llc),
                                           3:(Nr_out, llc, llc),
                                           4:(Nr_out, llc, n),
                                           5:(Nr_out, llc, n)}
            else:
                row_col_levels_per_face = {1: (1, n, llc),
                                           2: (1, n, llc),
                                           3: (1, llc, llc),
                                           4: (1, llc, n),
                                           5: (1, llc, n)}

            print('          - Running the downscaling scheme')
            N1_pickup_var_faces = df.downscale_3D_field_faces(var_grid_faces_extended, regridding_weights,
                                                              N0_wet_grid_3D_on_N1_faces, N1_wet_grid_3D_faces,
                                                              row_col_levels_per_face,
                                                              spread_horizontally=False)

            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                interp_grid_compact = sf.sassie_n1_faces_to_compact(N1_pickup_var_faces, llc, n)
            else:
                interp_grid_compact = sf.sassie_n1_faces_to_compact(N1_pickup_var_faces, llc, n, levels=1)

            # ###################################################
            # # Make a plot if desired
            # plot_faces = {}
            # for face in range(1,6):
            #     plot_faces[face] = L1_pickup_var_faces[face][0,:,:]
            # plot_grid = sf.stitch_faces_to_single_grid(plot_faces, llc, rows=680)
            #
            # C = plt.imshow(plot_grid, origin='lower')
            # plt.colorbar(C)
            # plt.title(var_name)
            # plt.show()

        else:
            if var_name.lower() not in ['etan','detahdt','etah']:
                interp_grid_compact = np.zeros((Nr_out,4*n+llc,llc))
            else:
                interp_grid_compact = np.zeros((1,4*n+llc,llc))

        interp_grids.append(interp_grid_compact)

        print('          - '+var_name+' compact shape: ',np.shape(interp_grid_compact))

    pickup_grid = stack_grids_to_pickup(var_names, interp_grids)
    print('    Total pickup grid shape: ', np.shape(pickup_grid))

    # pickup_grid = []

    output_dir = os.path.join('..', 'input')
    output_file = os.path.join(output_dir, 'pickup.0000000001')
    dtype = '>f8'
    write_pickup_file(output_file, dtype, pickup_grid, global_metadata)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()


    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_pickup_file(ecco_path)