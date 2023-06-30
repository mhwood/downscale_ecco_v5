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

def read_seaice_pickup_file_to_faces(input_init_dir,pickup_file='pickup_seaice.0000000001'):

    print('      Reading from '+input_init_dir + '/' + pickup_file)

    global_data, _, global_metadata = mds.rdmds(input_init_dir + '/' + pickup_file,
                                                returnmeta=True)

    var_names = []
    row_bounds = []
    all_var_grid_faces = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grid_faces = llc_compact_to_faces(var_grid[0,:,:], less_output=True)

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
        mask = mask[0,:,:]

        # plt.imshow(mask,origin='lower')
        # plt.title(nc_file_prefix+', '+str(face))
        # plt.show()

        mask_faces[face] = mask
        ds.close()
    return(mask_faces)

def read_regridder_functions_from_nc(N0_XC_faces, N0_YC_faces, N1_XC_faces, N1_YC_faces):

    regridding_faces = {}

    for face in range(1,6):
        ds_in = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N0_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N0_YC_faces[face]).T),
            ))

        ds_out = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N1_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N1_YC_faces[face]).T),
            ))

        filename = os.path.join('..', 'input', 'face_'+str(face)+'_downscale_weights.nc')
        regridder = xe.Regridder(ds_in, ds_out, method='bilinear', reuse_weights=True, filename = filename)
        regridding_faces[face] = regridder

    return(regridding_faces)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [   "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
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

def create_seaice_pickup_file(ecco_path):

    config_dir = '../..'
    LLC = 270
    N = 180

    llc = 1080
    n = 720

    # these are outlined stats of the different domains
    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    N1_size = size_dict['N1_' + str(llc)]
    n_rows_N1 = N1_size[0]
    n_cols_N1 = N1_size[1]

    print('Creating the pickup file for the N1_'+str(llc)+' model from the output of the N0_'+str(LLC)+' model')

    print('    - Reading in variables from the N0_'+str(LLC)+' pickup file')
    input_init_dir = os.path.join('..','..', 'N0_'+str(LLC), 'input')
    pickup_file = 'pickup_seaice.0000000001'
    var_names,row_bounds,all_var_grid_faces,global_metadata = read_seaice_pickup_file_to_faces(input_init_dir,pickup_file)

    print('    - Reading in the geometry of the N0_' + str(LLC) + ' domain')
    N0_XC_faces, N0_YC_faces = sf.read_sassie_n0_grid(ecco_path)

    N0_XC_faces_extended = {}
    N0_YC_faces_extended = {}
    for face in range(1, 6):
        N0_XC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_XC_faces, face)
        N0_YC_faces_extended[face] = sf.get_extended_var_grid_on_face(N0_YC_faces, face)

    print('    - Reading in the geometry of the N1_' + str(LLC) + ' domain')
    input_dir = os.path.join('..', 'input')
    N1_XC_faces, N1_YC_faces = sf.read_sassie_n1_grid(input_dir)

    print('    - Reading in the regridding functions')
    regridding_weights = df.read_regridder_functions_from_nc(input_dir, N0_XC_faces_extended, N0_YC_faces_extended,
                                                             N1_XC_faces, N1_YC_faces)

    print('    - Downscaling the pickup grids')
    interp_grids = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        print('      - Downscaling ' + var_name)
        if var_name not in []:

            if var_name in ['siVICE']:
                N1_wet_grid_2D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),
                                                               hFac='S')
                N0_wet_grid_2D_on_N1_faces = read_mask_faces_from_nc(
                    os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='S')
            elif var_name in ['siUICE']:
                N1_wet_grid_2D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),
                                                               hFac='W')
                N0_wet_grid_2D_on_N1_faces = read_mask_faces_from_nc(
                    os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='W')
            else:
                N1_wet_grid_2D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_extended_wet_grids'),
                                                               hFac='C')
                N0_wet_grid_2D_on_N1_faces = read_mask_faces_from_nc(
                    os.path.join('..', 'input', 'N0_extended_wet_grids_on_N1'), hFac='C')


            var_grid_faces = all_var_grid_faces[vn]

            for face in range(1, 6):
                if face < 3:
                    var_grid_faces[face] = var_grid_faces[face][-N:, :]
                if face > 3:
                    var_grid_faces[face] = var_grid_faces[face][:, :N]

            var_grid_faces_extended = {}
            for face in range(1, 6):
                var_grid_faces_extended[face] = sf.get_extended_var_grid_on_face(var_grid_faces, face)

            N1_pickup_var_faces = df.downscale_2D_field_faces(var_grid_faces_extended, regridding_weights,
                                                              N0_wet_grid_2D_on_N1_faces, N1_wet_grid_2D_faces,
                                                              spread_horizontally=False)


            interp_grid_compact = sf.sassie_n1_faces_to_compact(N1_pickup_var_faces,llc,n)

            interp_grid_compact = interp_grid_compact.reshape((1,
                                                               np.shape(interp_grid_compact)[0],
                                                               np.shape(interp_grid_compact)[1]))

            # ###################################################
            # # Make a plot if desired
            # plot_faces = {}
            # for face in range(1,6):
            #     plot_faces[face] = N1_pickup_var_faces[face]
            # plot_grid = sf.stitch_faces_to_single_grid(plot_faces, llc, rows=680)
            #
            # C = plt.imshow(plot_grid, origin='lower')
            # plt.colorbar(C)
            # plt.title(var_name)
            # plt.show()

        else:
            interp_grid_compact = np.zeros((1,n_rows_N1,n_cols_N1))

        interp_grids.append(interp_grid_compact)

        print('    '+var_name+' compact shape: ',np.shape(interp_grid_compact))

    pickup_grid = stack_grids_to_pickup(interp_grids)

    print('    Total pickup grid shape: ', np.shape(pickup_grid))

    output_dir = os.path.join('..', 'input')
    output_file = os.path.join(output_dir, 'pickup_seaice.0000000001')
    dtype = '>f8'
    write_pickup_file(output_file, dtype, pickup_grid, global_metadata)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()


    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_seaice_pickup_file(ecco_path)