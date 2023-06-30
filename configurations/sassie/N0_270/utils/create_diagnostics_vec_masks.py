import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import simplegrid as sg
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_faces_to_compact, llc_compact_to_faces
import argparse
import ast

def llc_270_blank_list():
    blank_list = [ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
                    25,26,27,29,30,31,32,33,34,35,36,39,40,41,42,43,44,45,54,133,142,150,151,159,160,161,166,
                    167,168,169,175,176,178,189,198,207,214,215,216,224,225,234,244,245,246,247,248,249,
                    250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,
                    270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,292,293,294,
                    295,296,297,368,369,417,418,419,420,424,425,426,427,428,429,433,434,435,436,437,438,442,
                    443,444,445,446,447,448,449,451,452,453,454,455,456,457,458,460,461,462,463,464,
                    465,466,467,469,470,471,472,473,474,475,476,477,478,481,482,483,484,485,486,495,
                    543,551,559,568,569,589,590,591,592,593,594,595,596,617,618,619,620,621,645,646,647,648,
                    673,674,675,700,701,702,727,728,729,731,754,755,756,758,780,781,782,783,784,785,807,
                    808,809,810,811,812,813,814,834,835,836,837,839,840,841,842,843,861,862,863,864,
                    866,867,868,869,870,871,888,889,
                    890,891,896,897,898,914,915,916,917,918,924,942,943,944,945,968,969,970,971,972,
                    983,984,985,986,987,996,997,998,999,1011,1012,1023,1024,1025,1026,1051,1052,1053]
    return(blank_list)

def read_global_XC_YC(ecco_dir,llc):

    grid_file_dir = os.path.join(ecco_dir,'mitgrid_tiles')
    # grid_file_dir = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in range(1,6):
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

    input_init_dir = os.path.join(ecco_dir, 'input_init')
    bathy_compact = np.fromfile(os.path.join(input_init_dir, 'bathy_llc' + str(llc)), '>f4')
    bathy_compact = np.reshape(bathy_compact, (13 * llc, llc))
    bathy_faces = llc_compact_to_faces(bathy_compact, less_output=True)

    return(XC_faces,YC_faces, bathy_faces)

def find_locations_in_faces(XC_faces, YC_faces,latitude):

    boundary_indices = {}

    for face in [1,2]:
        total_err = 1e10
        row_to_use = 0
        for row in range(np.shape(XC_faces[face])[0]):
            err = np.sum(np.abs((YC_faces[face][row,:]-latitude)))
            if err<total_err:
                total_err=err
                row_to_use = row
        boundary_indices[face] = row_to_use-1

    for face in [4,5]:
        total_err = 1e10
        col_to_use = 0
        for col in range(np.shape(XC_faces[face])[1]):
            err = np.sum(np.abs((YC_faces[face][:,col] - latitude)))
            if err < total_err:
                total_err = err
                col_to_use = col
        boundary_indices[face] = col_to_use+1

    return(boundary_indices)

def create_blank_cell_mask(bathy_faces):

    blank_cell_list = llc_270_blank_list()

    blank_faces = {}

    sNx = sNy = 30
    cell_counter = 1

    for face in range(1,6):
        bathy_face = bathy_faces[face]
        blank_face = np.ones_like(bathy_face).astype(int)
        nPx = int(np.shape(blank_face)[1] / sNx)
        nPy = int(np.shape(blank_face)[0] / sNy)

        for j in range(nPy):
            for i in range(nPx):
                if cell_counter in blank_cell_list:
                    blank_face[j  * sNy:(j  + 1) * sNy, i  * sNx:(i  + 1) * sNx] = 0
                cell_counter += 1

        # plt.subplot(1,2,1)
        # plt.imshow(bathy_face,origin='lower',cmap='Blues_r')
        # plt.title(str(face))
        # plt.subplot(1,2,2)
        # plt.imshow(blank_face,origin='lower',cmap='Blues_r')
        # plt.show()

        blank_faces[face] = blank_face
    return(blank_faces)

def create_boundary_mask(boundary_indices,layer,llc,bathy_faces,blank_faces):

    mask_faces = {}
    mask_faces[1] = np.zeros((3 * llc, llc))
    mask_faces[2] = np.zeros((3 * llc, llc))
    mask_faces[3] = np.zeros((llc, llc))
    mask_faces[4] = np.zeros((llc,3 * llc))
    mask_faces[5] = np.zeros((llc,3 * llc))

    mask_dict = np.zeros((llc*4, 3)).astype(int)

    transect_line = np.zeros((llc*4, ))

    counter = 1
    full_counter = 0
    for face in [1,2,4,5]:
        if face in [1]:# this face is special because we dont need the mediterranean
            for col in range(np.shape(mask_faces[face])[1]):
                if col<100:
                    if blank_faces[face][boundary_indices[face]+layer,col]==1 and bathy_faces[face][boundary_indices[face]+layer,col]<0:
                        mask_faces[face][boundary_indices[face]+layer,col] = counter
                        mask_dict[counter-1,0] = face
                        mask_dict[counter-1,1] = boundary_indices[face]+layer
                        mask_dict[counter-1,2] = col
                        transect_line[full_counter] = bathy_faces[face][boundary_indices[face]+layer,col]
                        # print(counter,face,boundary_indices[face]+layer,col)
                        counter += 1
                full_counter+=1
        if face in [2]:# this face is special because there is a sneaky blank cell
            for col in range(np.shape(mask_faces[face])[1]):
                if col>220:
                    if blank_faces[face][boundary_indices[face]+layer,col]==1 and bathy_faces[face][boundary_indices[face]+layer,col]<0:
                        mask_faces[face][boundary_indices[face]+layer,col] = counter
                        mask_dict[counter-1,0] = face
                        mask_dict[counter-1,1] = boundary_indices[face]+layer
                        mask_dict[counter-1,2] = col
                        transect_line[full_counter] = bathy_faces[face][boundary_indices[face]+layer,col]
                        # print(counter,face,boundary_indices[face]+layer,col)
                        counter += 1
                full_counter+=1
        if face in [4,5]:
            for row in range(np.shape(mask_faces[face])[0]):
                if blank_faces[face][row,boundary_indices[face]-layer]==1 and bathy_faces[face][row,boundary_indices[face]-layer]<0:
                    mask_faces[face][row,boundary_indices[face]-layer] = counter
                    # print(counter, face, row, boundary_indices[face]-layer)
                    mask_dict[counter-1, 0] = face
                    mask_dict[counter-1, 1] = row
                    mask_dict[counter-1, 2] = boundary_indices[face]-layer
                    transect_line[full_counter] = bathy_faces[face][row,boundary_indices[face]-layer]
                    counter += 1
                full_counter +=1

    mask_dict = mask_dict[:np.sum(mask_dict[:, 0] != 0), :]

    print('   This mask has '+str(counter-1)+' points')

    return(mask_faces,mask_dict,transect_line)


def create_surface_mask(boundary_indices,llc,bathy_faces,blank_faces):

    mask_faces = {}
    mask_faces[1] = np.zeros((3 * llc, llc))
    mask_faces[2] = np.zeros((3 * llc, llc))
    mask_faces[3] = np.zeros((llc, llc))
    mask_faces[4] = np.zeros((llc,3 * llc))
    mask_faces[5] = np.zeros((llc,3 * llc))

    mask_dict = np.zeros((llc*llc*5, 3))

    counter = 1
    for face in [1,2,3,4,5]:
        if face in [1,2]:
            for row in range(boundary_indices[face],np.shape(mask_faces[face])[0]):
                for col in range(np.shape(mask_faces[face])[1]):
                    if blank_faces[face][row,col]==1 and bathy_faces[face][row,col]<0:
                        mask_faces[face][row,col] = counter
                        mask_dict[counter-1,0] = face
                        mask_dict[counter-1,1] = row
                        mask_dict[counter-1,2] = col
                        counter += 1

        if face in [3]:
            for row in range(np.shape(mask_faces[face])[0]):
                for col in range(np.shape(mask_faces[face])[1]):
                    if blank_faces[face][row,col]==1 and bathy_faces[face][row,col]<0:
                        mask_faces[face][row,col] = counter
                        mask_dict[counter-1, 0] = face
                        mask_dict[counter-1, 1] = row
                        mask_dict[counter-1, 2] = col
                        counter += 1
        if face in [4,5]:
            for row in range(np.shape(mask_faces[face])[0]):
                for col in range(boundary_indices[face]):
                    if blank_faces[face][row,col]==1 and bathy_faces[face][row,col]<0:
                        mask_faces[face][row,col] = counter
                        mask_dict[counter-1, 0] = face
                        mask_dict[counter-1, 1] = row
                        mask_dict[counter-1, 2] = col
                        counter += 1

    mask_dict = mask_dict[:np.sum(mask_dict[:, 0] != 0), :]

    # print('   This mask has '+str(counter-1)+' points')

    return(mask_faces,mask_dict)


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


def plot_transect_lines(transect_lines,mask_names_list):

    fig = plt.figure(figsize=(8, 4))
    plt.style.use('dark_background')

    mask_colors = ['red','orange','green','purple']

    plt.plot([0, 0], [-9000, 1000],'w--',linewidth=0.6)
    plt.plot([270, 270], [-9000, 1000], 'w--',linewidth=0.6)
    plt.plot([540, 540], [-9000, 1000], 'w--',linewidth=0.6)
    plt.plot([810, 810], [-9000, 1000], 'w--',linewidth=0.6)
    plt.plot([1080, 1080], [-9000, 1000], 'w--',linewidth=0.6)

    for mn in range(len(mask_names_list)):
        plt.plot(transect_lines[mn],label=mask_names_list[mn],color=mask_colors[mn])

    plt.gca().set_xticks(np.arange(0,1080,135))


    plt.gca().set_ylim([-8000,500])

    plt.ylabel('Depth (m)')
    plt.xlabel('Grid cells along transect (including removed 0\'s)')
    plt.title('Bathymetry Profiles Along diagnostic_vec Masks')

    output_file = os.path.join('..','plots','dv_mask_bathymetry_profile')
    plt.savefig(output_file)
    plt.close(fig)


########################################################################################################################

def create_dv_masks(ecco_path,print_status):

    if print_status:
        print('Creating the diagnostics_vec masks to use in the L0 domain')

    if 'dv' not in os.listdir(os.path.join('..','input')):
        os.mkdir(os.path.join('..','input','dv'))

    llc = 270
    latitude = 41

    ###############################################################################################
    # Read in the grids

    if print_status:
        print('    Reading in the L0 domain files')

    # read the mitgrids to faces
    ecco_dir = os.path.join(ecco_path,'LLC'+str(llc)+'_Files')
    XC_faces, YC_faces, bathy_faces = read_global_XC_YC(ecco_dir, llc)

    blank_faces = create_blank_cell_mask(bathy_faces)

    # find the row (col) corresponding to 45N in faces 1 + 2 (4+5)
    boundary_indices = find_locations_in_faces(XC_faces, YC_faces,latitude)

    # add a buffer to include the boundary set up in Ians domain
    buffer_cells = 1
    boundary_indices[1] -= buffer_cells
    boundary_indices[2] -= buffer_cells
    boundary_indices[4] += buffer_cells
    boundary_indices[5] += buffer_cells

    all_mask_dicts = []
    mask_names_list = []
    transect_lines = []

    # make 4 layers of round-the-globe masks on faces
    for layer in range(4):
        if layer==0:
            mask_name = 'mask_'+str(latitude)+'N'
        else:
            mask_name = 'mask_'+str(latitude)+'N_i'+str(layer)
        print('  - Working on '+mask_name)
        mask_faces,mask_dict,transect_line = create_boundary_mask(boundary_indices, layer, llc, bathy_faces, blank_faces)

        # transect_lines.append(transect_line)
        # plt.plot(transect_line)
        # plt.show()

        print('     There are '+str(np.sum(mask_dict[:,0]==0))+' zero points out of '+str(np.shape(mask_dict)[0])+' points (should be none)')

        # convert the masks to compact and output as a binary file
        compact_mask = llc_faces_to_compact(mask_faces)
        # print(np.sum(compact_mask!=0),np.min(compact_mask[compact_mask!=0]),np.max(compact_mask))
        output_file = os.path.join('..','input','dv',mask_name+'.bin')
        compact_mask.ravel('C').astype('>f4').tofile(output_file)

        mask_names_list.append(mask_name)
        all_mask_dicts.append(mask_dict)

    print('  - Working on mask_arctic_surface')
    mask_faces, mask_dict = create_surface_mask(boundary_indices, llc, bathy_faces, blank_faces)

    print('     There are ' + str(np.sum(mask_dict[:, 0] == 0)) + ' zero points out of ' + str(
        np.shape(mask_dict)[0]) + ' points (should be none)')

    # convert the masks to compact and output as a binary file
    compact_mask = llc_faces_to_compact(mask_faces)
    output_file = os.path.join('..', 'input', 'dv', 'mask_arctic_surface.bin')
    compact_mask.ravel('C').astype('>f4').tofile(output_file)

    mask_names_list.append('mask_arctic_surface')
    all_mask_dicts.append(mask_dict)

    output_dir = os.path.join('..','input')
    output_file_name = 'dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list)

    # plot_transect_lines(transect_lines,mask_names_list)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    parser.add_argument("-p", "--print_status", action="store",
                        help="Print status of routine (1 for True, 0 for False).", dest="print_status",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    print_status = args.print_status

    if print_status>0:
        print_status=True
    else:
        print_status=False

    create_dv_masks(ecco_path,print_status)
