
import os
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_compact_to_faces

def sassie_n0_compact_to_faces(sassie_n0_compact, llc, n):
    sassie_faces = dict()

    dim = len(np.shape(sassie_n0_compact))
    Nr = 50

    # Face 1
    start_row = 0
    end_row = n
    if dim==2:
        sassie_faces[1] = sassie_n0_compact[start_row:end_row, :]
    if dim==3:
        sassie_faces[1] = sassie_n0_compact[:, start_row:end_row, :]

    # Face 2
    start_row = end_row
    end_row = start_row + n
    if dim == 2:
        sassie_faces[2] = sassie_n0_compact[start_row:end_row, :]
    if dim == 3:
        sassie_faces[2] = sassie_n0_compact[:, start_row:end_row, :]

    # Face 3
    start_row = end_row
    end_row = start_row + llc
    if dim == 2:
        sassie_faces[3] = sassie_n0_compact[start_row:end_row:, :]
    if dim == 3:
        sassie_faces[3] = sassie_n0_compact[:,start_row:end_row:, :]

    # Face 4
    start_row = end_row
    end_row = end_row + n
    if dim == 2:
        sassie_faces[4] = sassie_n0_compact[start_row:end_row,:].reshape(llc, n)
    if dim == 3:
        sassie_faces[4] = sassie_n0_compact[:,start_row:end_row,:].reshape(Nr,llc, n)

    # Face 5
    start_row = end_row
    end_row = end_row + n
    if dim == 2:
        sassie_faces[5] = sassie_n0_compact[start_row:end_row,:].reshape(llc, n)
    if dim == 3:
        sassie_faces[5] = sassie_n0_compact[:,start_row:end_row,:].reshape(Nr,llc, n)

    return sassie_faces


def sassie_n0_faces_to_compact(sassie_n0_faces, llc, n):
    Nr = 50

    if len(np.shape(sassie_n0_faces[1]))==2:
        sassie_n0_compact = np.vstack([sassie_n0_faces[1],
                                    sassie_n0_faces[2],
                                    sassie_n0_faces[3],
                                    sassie_n0_faces[4].reshape(n,llc),
                                    sassie_n0_faces[5].reshape(n,llc)])
    if len(np.shape(sassie_n0_faces[1]))==3:
        sassie_n0_compact = np.concatenate([sassie_n0_faces[1],
                                    sassie_n0_faces[2],
                                    sassie_n0_faces[3],
                                    sassie_n0_faces[4].reshape(Nr,n,llc),
                                    sassie_n0_faces[5].reshape(Nr,n,llc)],axis = 1)
    return sassie_n0_compact


def read_sassie_n0_grid(ecco_path):

    llc = 270
    n = 180

    XC_faces = {}
    YC_faces = {}

    for face in range(1, 6):

        llc_mitgrid_file = os.path.join(ecco_path, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles',
                                        'tile00' + str(face) + '.mitgrid')
        llc_mitgrid_grid = np.fromfile(llc_mitgrid_file, '>f8')
        if face in [1, 2]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16, 3 * llc + 1, llc + 1))
            llc_mitgrid_grid = llc_mitgrid_grid[:, -n - 1:, :]
        if face in [3]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16, llc + 1, llc + 1))
        if face in [4, 5]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16, llc + 1, 3 * llc + 1))
            llc_mitgrid_grid = llc_mitgrid_grid[:, :, :n + 1]

        XC_faces[face] = llc_mitgrid_grid[0,:-1,:-1]
        YC_faces[face] = llc_mitgrid_grid[1,:-1,:-1]

    return(XC_faces,YC_faces)


def read_sassie_n0_bathymetry(ecco_path):

    llc = 270
    n = 180

    bathy_file = os.path.join(ecco_path, 'LLC' + str(llc) + '_Files', 'input_init', 'bathy_llc' + str(llc))
    bathy_compact = np.fromfile(bathy_file, '>f4')
    bathy_compact = np.reshape(bathy_compact, (13 * llc, llc))
    bathy_faces = llc_compact_to_faces(bathy_compact, less_output=True)

    N0_bathy_faces = {}

    for face in range(1, 6):
        if face in [1, 2]:
            N0_bathy_faces[face] = bathy_faces[face][-n:, :]
        if face in [3]:
            N0_bathy_faces[face] = bathy_faces[face]
        if face in [4, 5]:
            N0_bathy_faces[face] = bathy_faces[face][:, :n]

    return(N0_bathy_faces)


def read_sassie_n1_grid(input_dir):
    llc = 1080
    n = 680

    XC_faces = {}
    YC_faces = {}

    for face in range(1, 6):

        llc_mitgrid_file = os.path.join(input_dir, 'tile00' + str(face) + '.mitgrid')
        llc_mitgrid_grid = np.fromfile(llc_mitgrid_file, '>f8')
        if face in [1, 2]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16, n + 1, llc + 1))
        if face in [3]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16, llc + 1, llc + 1))
        if face in [4, 5]:
            llc_mitgrid_grid = np.reshape(llc_mitgrid_grid, (16, llc + 1, n + 1))

        XC_faces[face] = llc_mitgrid_grid[0, :-1,:-1]
        YC_faces[face] = llc_mitgrid_grid[1, :-1,:-1]

    return (XC_faces, YC_faces)


def read_sassie_delR(domain_level):
    if domain_level==0:
        delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                         10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                         31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                         93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                         139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                         341.50, 364.50, 387.50, 410.50, 433.50, 456.50])
    else:
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
                         334.64, 365.93, 400.38, 438.23, 479.74])
    return(delR)


def sassie_n1_compact_to_faces(sassie_n1_compact, llc, n):
    sassie_faces = dict()

    dim = len(np.shape(sassie_n1_compact))
    Nr = 90

    # Face 1
    start_row = 0
    end_row = n #int(n*llc/tile_size)#
    # print(start_row, end_row, np.shape(sassie_n1_compact))
    if dim==2:
        sassie_faces[1] = sassie_n1_compact[start_row:end_row, :].reshape((n,llc))
    if dim==3:
        sassie_faces[1] = sassie_n1_compact[:, start_row:end_row, :]

    # Face 2
    start_row = end_row
    end_row = start_row + n
    # print(start_row, end_row, np.shape(sassie_n1_compact))
    if dim == 2:
        sassie_faces[2] = sassie_n1_compact[start_row:end_row, :].reshape((n,llc))
    if dim == 3:
        sassie_faces[2] = sassie_n1_compact[:, start_row:end_row, :]

    # Face 3
    start_row = end_row
    end_row = start_row + llc
    # print(start_row, end_row, np.shape(sassie_n1_compact))
    if dim == 2:
        sassie_faces[3] = sassie_n1_compact[start_row:end_row:, :].reshape((llc,llc))
    if dim == 3:
        sassie_faces[3] = sassie_n1_compact[:,start_row:end_row:, :]

    # Face 4
    start_row = end_row
    end_row = start_row + n
    # print(start_row, end_row, np.shape(sassie_n1_compact))
    if dim == 2:
        sassie_faces[4] = sassie_n1_compact[start_row:end_row,:].reshape((llc, n))
    if dim == 3:
        sassie_faces[4] = sassie_n1_compact[:,start_row:end_row,:].reshape(Nr,llc, n)

    # Face 5
    start_row = end_row
    end_row = end_row + n
    # print(start_row, end_row, np.shape(sassie_n1_compact))
    if dim == 2:
        sassie_faces[5] = sassie_n1_compact[start_row:end_row,:].reshape((llc, n))
    if dim == 3:
        sassie_faces[5] = sassie_n1_compact[:,start_row:end_row,:].reshape((Nr,llc, n))

    return sassie_faces


def sassie_n1_faces_to_compact(sassie_n1_faces, llc, n,levels=0):

    Nr = 90
    if levels!=0:
        Nr=levels

    if len(np.shape(sassie_n1_faces[1]))==2:
        sassie_n1_compact = np.vstack([sassie_n1_faces[1], \
                                    sassie_n1_faces[2], \
                                    sassie_n1_faces[3], \
                                    sassie_n1_faces[4].reshape(n,llc),\
                                    sassie_n1_faces[5].reshape(n,llc)])
    if len(np.shape(sassie_n1_faces[1]))==3:
        sassie_n1_compact = np.concatenate([sassie_n1_faces[1], \
                                    sassie_n1_faces[2], \
                                    sassie_n1_faces[3], \
                                    sassie_n1_faces[4].reshape(Nr,n,llc),\
                                    sassie_n1_faces[5].reshape(Nr,n,llc)],axis = 1)
    return sassie_n1_compact


def stitch_faces_to_single_grid(faces, llc, n):
    dim = len(np.shape(faces[1]))

    if dim == 2:
        stitched_grid = np.nan * np.ones((2* llc + n, 2* llc + n))

        stitched_grid[:n, :llc] = faces[1][-n:, :]
        stitched_grid[:n, llc:2*llc] = faces[2][-n:, :]
        stitched_grid[n:n+llc, llc:2*llc] = faces[3][:, :]
        stitched_grid[n:n+llc, 2*llc:] = faces[4][:, :n]
        stitched_grid[n+llc:n+2*llc, 2*llc:] = faces[5][:, :n]
    if dim == 3:
        stitched_grid = np.nan * np.ones((np.shape(faces[1])[0], 2 * llc + n, 2* llc + n))

        stitched_grid[:, :n, :llc] = faces[1][-n:, :]
        stitched_grid[:, :n, llc:2*llc] = faces[2][:, -n:, :]
        stitched_grid[:, n:n+llc, llc:2*llc] = faces[3][:, :, :]
        stitched_grid[:, n:n+llc, 2*llc:] = faces[4][:, :, :n]
        stitched_grid[:, n+llc:n+2*llc, 2*llc:] = faces[5][:, :, :n]

    return (stitched_grid)


def get_extended_var_grid_on_face(var_faces, face_number):

    var_face = var_faces[face_number]

    if face_number == 1 or face_number == 2:
        extended_var_face = np.zeros((np.shape(var_face)[0]+1,np.shape(var_face)[1]+2))
        extended_var_face[:-1, 1:-1] = var_face
        if face_number == 1:
            top = np.rot90(var_faces[3],k=3)
            left = np.rot90(var_faces[5])
            right = var_faces[2]
        if face_number == 2:
            top = var_faces[3]
            left = var_faces[1]
            right = np.rot90(var_faces[4])

        extended_var_face[:-1,0] = left[:,-1]
        extended_var_face[:-1, -1] = right[:, -1]
        extended_var_face[-1,1:-1] = top[0,:]

        # plt.subplot(2,3,2)
        # plt.imshow(top,origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(2, 3, 4)
        # plt.imshow(left, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(2, 3, 5)
        # plt.imshow(var_face, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(2, 3, 6)
        # plt.imshow(right, origin='lower',vmin=-10,vmax = 0)
        # plt.show()

    if face_number == 3:
        extended_var_face = np.zeros((np.shape(var_face)[0]+2,np.shape(var_face)[1]+2))
        extended_var_face[1:-1, 1:-1] = var_face

        top = np.rot90(var_faces[5],k=3)
        right = var_faces[4]
        left = np.rot90(var_faces[1])
        bottom = var_faces[2]

        extended_var_face[1:-1,0] = left[:,-1]
        extended_var_face[0, 1:-1] = bottom[-1, :]
        extended_var_face[-1,1:-1] = top[0,:]
        extended_var_face[1:-1, -1] = right[:, 0]

        # plt.subplot(3,3,2)
        # plt.imshow(top,origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 4)
        # plt.imshow(left, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 5)
        # plt.imshow(var_face, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 8)
        # plt.imshow(bottom, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 3, 6)
        # plt.imshow(right, origin='lower', vmin=-10, vmax=0)
        # plt.show()

    if face_number == 4 or face_number == 5:
        extended_var_face = np.zeros((np.shape(var_face)[0]+2,np.shape(var_face)[1]+1))
        extended_var_face[1:-1, 1:] = var_face
        if face_number == 4:
            top = var_faces[5]
            left = var_faces[3]
            bottom = np.rot90(var_faces[2],k=3)
        if face_number == 5:
            top = np.rot90(var_faces[1],k=3)
            left = np.rot90(var_faces[3])
            bottom = var_faces[4]

        extended_var_face[1:-1,0] = left[:,-1]
        extended_var_face[0, 1:] = bottom[-1, :]
        extended_var_face[-1,1:] = top[0,:]

        # plt.subplot(3,2,2)
        # plt.imshow(top,origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 2, 3)
        # plt.imshow(left, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 2, 4)
        # plt.imshow(var_face, origin='lower',vmin=-10,vmax = 0)
        # plt.subplot(3, 2, 6)
        # plt.imshow(bottom, origin='lower',vmin=-10,vmax = 0)
        # plt.show()
    return(extended_var_face)


def get_extended_var_grid_on_face_3D(var_faces, face_number, levels):

    var_face = var_faces[face_number]

    if face_number == 1 or face_number == 2:
        extended_var_face = np.zeros((levels,np.shape(var_face)[1]+1,np.shape(var_face)[2]+2))
    if face_number == 3:
        extended_var_face = np.zeros((levels,np.shape(var_face)[1]+2,np.shape(var_face)[2]+2))
    if face_number == 4 or face_number == 5:
        extended_var_face = np.zeros((levels,np.shape(var_face)[1]+2,np.shape(var_face)[2]+1))

    for ll in range(levels):
        level_var_faces = {}
        for face in range(1,6):
            level_var_faces[face] = var_faces[face][ll, :, :]
        extended_var_face[ll, :, :] = get_extended_var_grid_on_face(level_var_faces, face_number)

    return(extended_var_face)


def get_extended_var_grid_on_face_4D(var_faces, face_number, levels, timesteps):

    var_face = var_faces[face_number]

    if face_number == 1 or face_number == 2:
        extended_var_face = np.zeros((timesteps, levels,np.shape(var_face)[1]+1,np.shape(var_face)[2]+2))
    if face_number == 3:
        extended_var_face = np.zeros((timesteps, levels,np.shape(var_face)[1]+2,np.shape(var_face)[2]+2))
    if face_number == 4 or face_number == 5:
        extended_var_face = np.zeros((timesteps, levels,np.shape(var_face)[1]+2,np.shape(var_face)[2]+1))

    for t in range(timesteps):
        for ll in range(levels):
            level_var_faces = {}
            for face in range(1,6):
                level_var_faces[face] = var_faces[face][t, ll, :, :]
            extended_var_face[t, ll, :, :] = get_extended_var_grid_on_face(level_var_faces, face_number)

    return(extended_var_face)
