import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import argparse


def write_grid_to_nc(config_dir, var_name, full_grid):

    if 'beaufort' not in os.listdir(os.path.join(config_dir,'N1_1080','results')):
        os.mkdir(os.path.join(config_dir,'N1_1080','results','beaufort'))
    if var_name not in os.listdir(os.path.join(config_dir,'N1_1080','results','beaufort')):
        os.mkdir(os.path.join(config_dir,'N1_1080','results','beaufort',var_name))

    ds = nc4.Dataset(os.path.join(config_dir,'N1_1080','results','beaufort',var_name,'beaufort_'+var_name+'.nc'),'w')

    ds.createDimension('time',np.shape(full_grid)[0])
    ds.createDimension('depth', np.shape(full_grid)[1])
    ds.createDimension('row', np.shape(full_grid)[2])
    ds.createDimension('col', np.shape(full_grid)[3])

    var = ds.createVariable(var_name,'f4',('time','depth','row','col'))
    var[:,:,:,:] = full_grid
    ds.close()

    a=1

def create_beaufort_file(config_dir):

    var_names = ['THETA']

    mask_file = os.path.join(config_dir,'N1_1080','input','N1_1080_dv_mask_ref.nc')
    ds = nc4.Dataset(mask_file)
    grp = ds.groups['beaufort']
    faces = grp.variables['source_faces'][:]
    cols = grp.variables['source_cols'][:]
    rows = grp.variables['source_rows'][:]
    ds.close()

    Nr = 90
    llc = 1080
    n = 680

    for var_name in var_names:
        file_name = 'beaufort_mask_'+var_name+'.0005785921.bin'
        file_path = os.path.join(config_dir,'N1_1080','run','dv',file_name)
        file_grid = np.fromfile(file_path,'>f4')
        n_timesteps = int(np.size(file_grid)/(len(rows)*Nr))
        file_grid = np.reshape(file_grid,(n_timesteps,Nr,len(faces)))

        grid_started = False

        for timestep in range(n_timesteps):
            grid_faces = {}
            grid_faces[3] = np.zeros((Nr, llc, llc))
            grid_faces[4] = np.zeros((Nr, llc, n))

            for i in range(len(faces)):
                grid_faces[faces[i]][:, rows[i],cols[i]] = file_grid[timestep, :, i]

            file_grid = np.concatenate([grid_faces[3],grid_faces[4]],axis=2)

            if not grid_started:
                nonzero_rows, nonzero_cols = np.where(file_grid[0,:,:]!=0)
                min_row = np.min(nonzero_rows)
                max_row = np.max(nonzero_rows)
                min_col = np.min(nonzero_cols)
                max_col = np.max(nonzero_cols)

            file_grid = file_grid[:, min_row:max_row+1, min_col:max_col+1]

            if not grid_started:
                max_depth_level = 1
                for k in range(np.shape(file_grid)[0]):
                    if np.any(file_grid[k,:,:]!=0):
                        max_depth_level = k

            print(max_depth_level)
            file_grid = file_grid[:max_depth_level+1,:,:]

            file_grid = np.rot90(file_grid,axes=(1,2))
            file_grid = file_grid[np.newaxis, ...]

            if not grid_started:
                grid_started = True
                full_grid = file_grid
            else:
                full_grid = np.concatenate([full_grid,file_grid])

            print(np.shape(full_grid))
            # plt.imshow(file_grid[0,10,:,:],origin='lower')
            # plt.show()

        write_grid_to_nc(config_dir, var_name, full_grid)





########################################################################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="Path to the configurations.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_beaufort_file(config_dir)