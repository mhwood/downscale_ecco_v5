
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces
import argparse
import ast
import sys
# sys.path.insert(1,'/Users/michwood/Documents/Research/Projects/Ocean_Modeling/Downscale_Baffin_Bay/MITgcm/configurations/baffin_bay_ECCOv4r4_forcing/utils')
sys.path.insert(1,os.path.join('..','..','utils'))
import downscale_functions as df
import sassie_functions as sf


def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[mask_name]
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()
    return(faces,rows,cols)

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

def read_L1_faces_geometry(llc,rows,cols):
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

def read_mask_faces_from_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['wet_grid_'+hFac][:,:,:]
    # print(np.shape(mask))
    mask = mask[0,:,:]
    # print(np.shape(mask))
    mask_faces = sf.sassie_n1_compact_to_faces(mask,dim=2,Nr=90)
    ds.close()
    return(mask_faces)

def read_surface_variables_onto_L0_grid(L0_run_dir, var_name, nTimesteps, nTimestepsOut, faces, rows, cols, llc):

    print('       - Reading L0 exf field ' + var_name +' output onto N0 domain ')

    # make a blank grid of zeros
    grid_faces = {}
    for i in [1,2,3,4,5]:
        if i < 3:
            grid_face = np.zeros((nTimestepsOut,llc*3,llc))
        if i == 3:
            grid_face = np.zeros((nTimestepsOut,llc,llc))
        if i > 3:
            grid_face = np.zeros((nTimestepsOut,llc,llc*3))
        grid_faces[i] = grid_face

    var_file = os.path.join(L0_run_dir, 'mask_arctic_surface_' + var_name + '.bin')
    var_grid = np.fromfile(var_file, dtype='>f4')
    N = int(np.size(var_grid)/nTimesteps)
    var_grid = np.reshape(var_grid,(nTimesteps,N))
    var_grid = var_grid[:nTimestepsOut,:]

    for n in range(N):
        grid_faces[faces[n]][:,rows[n],cols[n]] = var_grid[:,n]

    return(grid_faces)

def downscale_L0_exf_field_to_L1(surface_var_name, nTimestepsOut,
                                  L0_XC_faces, L0_YC_faces, L0_surface_var_faces,
                                  L1_XC_faces, L1_YC_faces, L1_wet_grid_3D_faces):

    # make a grid of zeros to fill in for each face
    L1_surface_var_faces = {}
    for face in range(1,6):
        L1_surface_var_faces[face] = np.zeros((nTimestepsOut,
                                               np.shape(L1_XC_faces[face])[0],np.shape(L1_XC_faces[face])[1]))

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

    for timestep in range(nTimestepsOut):
        # if timestep % 50 == 0:
        print('          Working on timestep ' + str(timestep+1) + ' of ' + str(nTimestepsOut))
        for face in range(1,6):
            # print('            Working on face ' + str(face) + ' of 5')
            L0_XC = L0_XC_faces[face]
            L0_YC = L0_YC_faces[face]
            L0_surface_var = L0_surface_var_faces[face][timestep,:,:]
            L1_XC = L1_XC_faces[face]
            L1_YC = L1_YC_faces[face]
            L1_wet_grid = L1_wet_grid_3D_faces[face]
            downscaled_field = df.downscale_2D_field(L0_XC, L0_YC, L0_surface_var,
                                                     L0_wet_grid, L0_wet_grid_on_L1,
                                                     L1_XC, L1_YC, L1_wet_grid)
            L1_surface_var_faces[face][timestep,:,:] = downscaled_field

    return(L1_surface_var_faces)

def output_exf_variable(output_dir, var_grid_faces, var_name):

    var_grid_compact = sf.sassie_n1_faces_to_compact(var_grid_faces)

    print('   Creating exf field '+var_name+' of shape '+str(np.shape(var_grid_compact)))
    output_file = os.path.join(output_dir,'L1_exf_' + var_name + '.bin')
    var_grid_compact.ravel('C').astype('>f4').tofile(output_file)


########################################################################################################################

def create_exf_fields_via_interpolation(nTimesteps,nTimestepsOut,ecco_path):

    config_dir = '../..'
    LLC = 270
    llc = 1080

    # these are outlined stats of the different domains
    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L1_size = size_dict['N1_'+str(llc)]
    n_rows_L1 = L1_size[0]
    n_cols_L1 = L1_size[1]

    # this is the dir where the exf output will be stored
    if 'exf' not in os.listdir(os.path.join(config_dir,'N1_'+str(llc),'input')):
        os.mkdir(os.path.join(config_dir,'N1_'+str(llc),'input','exf'))
    output_dir = os.path.join(config_dir,'N1_'+str(llc),'input','exf')

    print('Creating the exf fields for the N1_'+str(llc)+' model from the output of the N0_'+str(LLC)+' model')

    print('    - Reading the mask to reference the variable to the LLC grid')
    nc_dict_file = os.path.join(config_dir,'N0_'+str(LLC),'input','dv_mask_reference_dict.nc')
    faces, rows, cols = read_mask_reference_from_nc_dict(nc_dict_file, 'mask_arctic_surface')

    print('    - Reading in the N0_'+str(LLC)+' exf diagnostics_vec output')
    L0_run_dir = os.path.join(config_dir, 'N0_'+str(LLC), 'run', 'dv')
    surface_var_name_set = ['APRESS', 'AQH', 'ATEMP', 'LWDOWN', 'PRECIP', 'SWDOWN',
                            'USTRESS', 'VSTRESS', 'WSPEED', 'RUNOFF']
    # surface_var_name_set = ['WSPEED']

    print('    - Reading in the geometry of the L0_'+str(LLC)+' domain')
    L0_XC_faces, L0_YC_faces = read_L0_faces_geometry(ecco_path,LLC)

    print('    - Reading in the geometry of the L1_'+str(llc)+' domain')
    L1_XC_faces, L1_YC_faces, N = read_L1_faces_geometry(llc,n_rows_L1,n_cols_L1)

    print('    - Limiting geometry of the N0_' + str(LLC) + ' domain')
    for face in range(1,6):
        if face<3:
            L0_XC_faces[face] = L0_XC_faces[face][-int(N / 4):, :]
            L0_YC_faces[face] = L0_YC_faces[face][-int(N / 4):, :]
        if face>3:
            L0_XC_faces[face] = L0_XC_faces[face][:,:int(N / 4):]
            L0_YC_faces[face] = L0_YC_faces[face][:,:int(N / 4):]

    print('    - Downscaling the exf fields')
    for ss in range(len(surface_var_name_set)):

        L0_surface_var_faces = read_surface_variables_onto_L0_grid(L0_run_dir, surface_var_name_set[ss], nTimesteps,nTimestepsOut,
                                                             faces,rows, cols, LLC)

        for face in range(1, 6):
            if face < 3:
                L0_surface_var_faces[face] = L0_surface_var_faces[face][:,-int(N / 4):, :]
            if face > 3:
                L0_surface_var_faces[face] = L0_surface_var_faces[face][:,:, :int(N / 4):]

        if surface_var_name_set[ss] == 'VSTRESS':
            L1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..','input', 'N1_wet_grids.nc'), hFac='S')
        elif surface_var_name_set[ss] == 'USTRESS':
            L1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..','input', 'N1_wet_grids.nc'), hFac='W')
        else:
            L1_wet_grid_3D_faces = read_mask_faces_from_nc(os.path.join('..', 'input', 'N1_wet_grids.nc'), hFac='C')

        surface_var_name = surface_var_name_set[ss]
        print('        - Creating the '+surface_var_name+' conditions')

        L1_surface_var_faces = downscale_L0_exf_field_to_L1(surface_var_name, nTimestepsOut,
                                                      L0_XC_faces, L0_YC_faces, L0_surface_var_faces,
                                                      L1_XC_faces, L1_YC_faces, L1_wet_grid_3D_faces)

        # ###################################################
        # # Make a plot if desired
        # plot_timestep = 2
        # plot_faces = {}
        # for face in range(1,6):
        #     plot_faces[face] = L1_surface_var_faces[face][plot_timestep,:,:]
        # plot_grid = sf.stitch_faces_to_single_grid(plot_faces, llc, rows=680)
        #
        # C = plt.imshow(plot_grid, origin='lower')
        # plt.colorbar(C)
        # plt.title(surface_var_name)
        # plt.show()

        output_exf_variable(output_dir, L1_surface_var_faces, surface_var_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--nTimesteps", action="store",
                        help="The number of timesteps in the diagnostics_vec output from the N0_270 domain.",
                        dest="nTimesteps",
                        type=int, required=True)

    parser.add_argument("-o", "--nTimestepsOut", action="store",
                        help="The number of timesteps to output for the N1_1080 domain (counting by steps in the N0_270 domain).", dest="nTimestepsOut",
                        type=int, required=False, default=-1)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default='../../ECCO')

    args = parser.parse_args()
    nTimesteps = args.nTimesteps
    nTimestepsOut = args.nTimestepsOut
    ecco_path = args.ecco_path

    if nTimestepsOut==-1:
        nTimestepsOut = nTimesteps

    create_exf_fields_via_interpolation(nTimesteps,nTimestepsOut,ecco_path)