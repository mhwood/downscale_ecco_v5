
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse

########################################################################################################################

def read_reference_dict(mask_name):
    ds = nc4.Dataset(os.path.join('..','input','dv_mask_reference_dict.nc'))
    grp = ds.groups[mask_name]
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    return(faces,rows,cols)

def read_grid_to_entire_globe(grid,llc,faces,rows,cols):

    full_grid = np.zeros((np.shape(grid)[0],4*llc))

    for i in range(np.shape(grid)[1]):
        if faces[i]==1:
            full_grid[:,cols[i]] = grid[:,i]
        if faces[i]==2:
            full_grid[:,llc+cols[i]] = grid[:,i]
        if faces[i]==4:
            full_grid[:,2*llc+rows[i]] = grid[:,i]
        if faces[i]==5:
            full_grid[:,3*llc+rows[i]] = grid[:,i]
    return(full_grid)


def plot_dv_output(mask_name,var_name,n_timesteps,plot_timestep,mask_len):

    llc = 270
    Nr = 50

    if 'plots' not in os.listdir(os.path.join('..')):
        os.mkdir(os.path.join('..','plots'))

    file_path = os.path.join('..','run','dv',mask_name+'_'+var_name+'.bin')
    grid = np.fromfile(file_path,'>f4')

    grid = np.reshape(grid,(n_timesteps,Nr,mask_len))
    grid = grid[plot_timestep,:,:]

    faces, rows, cols = read_reference_dict(mask_name)
    grid = read_grid_to_entire_globe(grid, llc, faces, rows, cols)

    fig = plt.figure(figsize=(15, 4))
    plt.style.use('dark_background')

    aspect = mask_len / (3*Nr)
    C = plt.imshow(grid,aspect = aspect)
    plt.colorbar(C)

    plt.ylabel('Depth Levels')
    plt.xlabel('Grid Points Along Mask')
    plt.title(var_name+' on mask '+mask_name+' (timestep '+str(plot_timestep+1)+' of '+str(n_timesteps)+')')

    output_file = os.path.join('..', 'plots', mask_name+'_'+var_name+'.png')
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--mask_name", action="store",
                        help="Name of the mask to plot.", dest="mask_name",
                        type=str, required=True)

    parser.add_argument("-v", "--var_name", action="store",
                        help="Name of the variable to plot.", dest="var_name",
                        type=str, required=True)

    parser.add_argument("-n", "--n_timesteps", action="store", help="The number of timesteps in the output.",
                        dest="n_timesteps", type=int, required=True)

    parser.add_argument("-p", "--plot_timestep", action="store", help="The timestep to plot.",
                        dest="plot_timestep", type=int, required=True)

    parser.add_argument("-l", "--mask_len", action="store", help="The length of the mask.",
                        dest="mask_len", type=int, required=True)

    args = parser.parse_args()
    mask_name = args.mask_name
    var_name = args.var_name
    n_timesteps = args.n_timesteps
    plot_timestep = args.plot_timestep
    mask_len = args.mask_len


    plot_dv_output(mask_name,var_name,n_timesteps,plot_timestep,mask_len)