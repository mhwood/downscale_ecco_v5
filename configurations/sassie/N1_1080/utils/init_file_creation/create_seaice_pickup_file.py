import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
from pyproj import Transformer
from MITgcmutils import mds
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
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
        var_grid_faces = llc_compact_to_faces(var_grid, less_output=True)

        all_var_grid_faces.append(var_grid_faces)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,all_var_grid_faces,global_metadata)

def read_L0_faces_geometry(ecco_dir,llc):
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

def read_L1_faces_geometry(llc,rows):
    N = 680

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

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 4326 and outputCRS == 3411:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

def reproject_faces_to_polar(XC_faces, YC_faces):

    polar_XC_faces = {}
    polar_YC_faces = {}
    for face in [1,2,3,4,5]:
        points = np.column_stack([XC_faces[face].ravel(),YC_faces[face].ravel()])
        polar_points = reproject_points(points, 4326, 3411)
        polar_XC_faces[face] = polar_points[:, 0].reshape(np.shape(XC_faces[face]))
        polar_YC_faces[face] = polar_points[:, 1].reshape(np.shape(YC_faces[face]))

        # plt.subplot(1,2,1)
        # C = plt.imshow(polar_XC_faces[face])
        # plt.colorbar(C)
        #
        # plt.subplot(1, 2, 2)
        # C = plt.imshow(polar_YC_faces[face])
        # plt.colorbar(C)
        #
        # plt.show()

    return(polar_XC_faces, polar_YC_faces)

def read_L0_wetgrid_faces(config_dir,hFac):

    input_dir = os.path.join(config_dir,'N0_270','input')

    N = 680

    file_name = 'hFac'+hFac#+'.data'

    grid_compact = mds.rdmds(input_dir + '/' + file_name, returnmeta=False)
    grid_faces = llc_compact_to_faces(grid_compact,less_output=True)

    for face in range(1,6):
        if face<3:
            grid_faces[face] = grid_faces[face][:,-int(N / 4):, :]
        if face>3:
            grid_faces[face] = grid_faces[face][:,:,:int(N / 4):]

    # if fn==2:
    #     for face in range(1, 6):
    #         C = plt.imshow(grid_faces[face][0, :, :], origin='lower')
    #         plt.colorbar(C)
    #         plt.title(file_name+', face ' + str(face))
    #         plt.show()


    return(grid_faces)

def read_L1_wetgrid_faces(sf,config_dir,hFac,llc):

    input_dir = os.path.join(config_dir,'N1_1080','input')

    n = 680
    Nr = 90

    file_name = 'hFac'+hFac+'.data'

    grid_compact = np.fromfile(os.path.join(input_dir, file_name), '>f4')
    if 'hFac' in file_name:
        grid_compact = np.reshape(grid_compact, (Nr, 4 * n + llc, llc))
    else:
        grid_compact = np.reshape(grid_compact, (4 * n + llc, llc))

    grid_compact[grid_compact>0] = 1

    grid_faces = sf.sassie_n1_compact_to_faces(grid_compact, llc, n)

    # if fn==2:
    #     for face in range(1, 6):
    #         C = plt.imshow(grid_faces[face][0, :, :], origin='lower')
    #         plt.colorbar(C)
    #         plt.title(file_name+', face ' + str(face))
    #         plt.show()

    delR = sf.read_sassie_delR(domain_level=1)

    return(grid_faces)

def read_mask_faces_from_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['wet_grid_'+hFac][:,:,:]
    mask_faces = sf.sassie_n1_compact_to_faces(mask,dim=3,Nr=90)
    ds.close()
    return(mask_faces)

def downscale_L0_pickup_field_to_L1(surface_var_name,
                                  L0_XC_faces, L0_YC_faces, L0_pickup_var_faces, L0_wet_grid_3D_faces,
                                  L1_XC_faces, L1_YC_faces, L1_wet_grid_3D_faces):

    # make a grid of zeros to fill in for each face
    L1_pickup_var_faces = {}

    # print('    Variable shapes:')
    # print('        L0_XC: '+str(np.shape(L0_XC)))
    # print('        L0_YC: ' + str(np.shape(L0_YC)))
    # print('        L0_surface_var: ' + str(np.shape(L0_surface_var)))
    # print('        L0_wet_grid: ' + str(np.shape(L0_wet_grid)))
    # print('        XC_subset: ' + str(np.shape(XC_subset)))
    # print('        YC_subset: ' + str(np.shape(YC_subset)))

    # plt.subplot(1,2,1)
    # C = plt.imshow(L0_surface_var[5,:,:])
    # plt.colorbar(C)
    # plt.title(surface_var_name)
    # plt.subplot(1,2,2)
    # plt.imshow(L0_wet_grid_subset)
    # plt.show()

    L0_wet_grid = []
    L0_wet_grid_on_L1 = []

    for face in range(1,6):
        print('            Working on face ' + str(face) + ' of 5')
        L0_XC = L0_XC_faces[face]
        L0_YC = L0_YC_faces[face]
        L0_pickup_var = L0_pickup_var_faces[face]
        L0_wet_grid = L0_wet_grid_3D_faces[face]

        L1_XC = L1_XC_faces[face]
        L1_YC = L1_YC_faces[face]
        L1_wet_grid = L1_wet_grid_3D_faces[face]
        L0_wet_grid_on_L1 = np.copy(L1_wet_grid)

        # print('    Variable shapes:')
        # print('        L0_XC: '+str(np.shape(L0_XC)))
        # print('        L0_YC: ' + str(np.shape(L0_YC)))
        # print('        L0_pickup_var: ' + str(np.shape(L0_pickup_var)))
        # print('        L0_wet_grid: ' + str(np.shape(L0_wet_grid)))
        # print('        L1_XC: ' + str(np.shape(L1_XC)))
        # print('        L1_YC: ' + str(np.shape(L1_YC)))
        # print('        L1_wet_grid: ' + str(np.shape(L1_wet_grid)))

        L0_wet_grid = L0_wet_grid[:1,:,:]

        downscaled_field = df.downscale_3D_field_with_zeros(L0_XC, L0_YC, L0_pickup_var, L0_wet_grid, L0_wet_grid_on_L1,
                                                            L1_XC, L1_YC, L1_wet_grid,
                                                            mean_vertical_difference=0, fill_downward=True,
                                                            printing=True)

        # plt.imshow(downscaled_field[0,:,:],origin='lower')
        # plt.title(surface_var_name)
        # plt.show()

        # downscaled_field = df.downscale_3D_field(L0_XC, L0_YC, L0_pickup_var,
        #                                          L0_wet_grid, L0_wet_grid_on_L1,
        #                                          L1_XC, L1_YC, L1_wet_grid)
        L1_pickup_var_faces[face] = downscaled_field

    return(L1_pickup_var_faces)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_seaice_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')
    print(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [    "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [    "+str(subset_metadata['nflds'][0])+" ];\n"
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

    config_dir = '../..'
    LLC = 270
    llc = 1080
    n = 680

    n_rows_L1 = llc+4*n
    n_cols_L1 = llc

    print('Creating the pickup file for the N1_'+str(llc)+' model from the output of the N0_'+str(LLC)+' model')

    print('    - Reading in variables from the N0_'+str(LLC)+' pickup file')
    input_init_dir = os.path.join('..','..', 'N0_'+str(LLC), 'run')
    pickup_file = 'pickup_seaice.0000578592'
    var_names,row_bounds,all_var_grid_faces,global_metadata = read_seaice_pickup_file_to_faces(input_init_dir,pickup_file)

    # for i in range(2,3):#len(pickup_var_names)):
    #     C = plt.imshow(var_grid_faces[i][3][0, :, :],origin='lower')
    #     plt.title(var_names[i])
    #     plt.colorbar(C)
    #     plt.show()

    print('    - Reading in the geometry of the L0_'+str(LLC)+' domain')
    L0_XC_faces, L0_YC_faces = read_L0_faces_geometry(ecco_path,LLC)
    L0_XC_faces, L0_YC_faces = reproject_faces_to_polar(L0_XC_faces, L0_YC_faces)

    print('    - Reading in the geometry of the L1_'+str(llc)+' domain')
    L1_XC_faces, L1_YC_faces, N = read_L1_faces_geometry(llc,n_rows_L1)
    L1_XC_faces, L1_YC_faces = reproject_faces_to_polar(L1_XC_faces, L1_YC_faces)

    print('    - Limiting geometry of the N0_' + str(LLC) + ' domain')
    for face in range(1,6):
        if face<3:
            L0_XC_faces[face] = L0_XC_faces[face][-int(N / 4):, :]
            L0_YC_faces[face] = L0_YC_faces[face][-int(N / 4):, :]
        if face>3:
            L0_XC_faces[face] = L0_XC_faces[face][:,:int(N / 4):]
            L0_YC_faces[face] = L0_YC_faces[face][:,:int(N / 4):]

    print('    - Downscaling the pickup grids')
    interp_grids = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        print('      - Downscaling ' + var_name)
        if var_name not in []:

            if var_name in ['siVICE']:
                L0_wet_grid_3D_faces = read_L0_wetgrid_faces(config_dir,hFac='S')
                L1_wet_grid_3D_faces = read_L1_wetgrid_faces(sf,config_dir,hFac='S',llc=llc)
            elif var_name in ['siUICE']:
                L0_wet_grid_3D_faces = read_L0_wetgrid_faces(config_dir,hFac='W')
                L1_wet_grid_3D_faces = read_L1_wetgrid_faces(sf,config_dir,hFac='W',llc=llc)
            else:
                L0_wet_grid_3D_faces = read_L0_wetgrid_faces(config_dir,hFac='C')
                L1_wet_grid_3D_faces = read_L1_wetgrid_faces(sf,config_dir,hFac='C',llc=llc)

            var_grid_faces = all_var_grid_faces[vn]

            for face in range(1, 6):
                if face < 3:
                    var_grid_faces[face] = var_grid_faces[face][:,-int(N / 4):, :]
                if face > 3:
                    var_grid_faces[face] = var_grid_faces[face][:,:, :int(N / 4):]

        #     print('        Variable shapes:')
        #     print('            llc_XC_pacific_subset: '+str(np.shape(llc_XC_pacific_subset)))
        #     print('            llc_YC_pacific_subset: ' + str(np.shape(llc_YC_pacific_subset)))
        #     print('            llc_var_pacific_subset: ' + str(np.shape(llc_var_pacific_subset)))
        #     print('            L0_wet_cells_3D_subset: ' + str(np.shape(L0_wet_cells_3D_subset)))
        #     print('            L0_wet_cells_on_L1_3D: ' + str(np.shape(L0_wet_cells_on_L1_3D)))
        #     print('            XC_subset: ' + str(np.shape(XC_subset)))
        #     print('            YC_subset: ' + str(np.shape(YC_subset)))
        #     print('            L1_wet_cells_3D ' + str(np.shape(L1_wet_cells_3D)))

            interp_field_faces = downscale_L0_pickup_field_to_L1(var_name,
                                  L0_XC_faces, L0_YC_faces, var_grid_faces, L0_wet_grid_3D_faces,
                                  L1_XC_faces, L1_YC_faces, L1_wet_grid_3D_faces)

            interp_grid_compact = sf.sassie_n1_faces_to_compact(interp_field_faces, llc, n, levels=1)

            # ###################################################
            # # Make a plot if desired
            # plot_faces = {}
            # for face in range(1,6):
            #     plot_faces[face] = interp_field_faces[face][0,:,:]
            # plot_grid = sf.stitch_faces_to_single_grid(plot_faces, llc, rows=680)
            #
            # C = plt.imshow(plot_grid, origin='lower')
            # plt.colorbar(C)
            # plt.title(var_name)
            # plt.show()

        else:
            interp_grid_compact = np.zeros((1,n_rows_L1,n_cols_L1))

            print(np.shape(interp_grid_compact))

        interp_grids.append(interp_grid_compact)

    pickup_grid = stack_grids_to_pickup(interp_grids)
    print(np.shape(pickup_grid))

    output_dir = os.path.join('..', 'input')
    output_file = os.path.join(output_dir, 'pickup_seaice.0005785920')
    dtype = '>f8'
    global_metadata['nrecords'] = [np.shape(pickup_grid)[0]]
    write_seaice_pickup_file(output_file, dtype, pickup_grid, global_metadata)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_pickup_file(ecco_path)