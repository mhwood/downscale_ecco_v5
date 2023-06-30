
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

def read_grid_to_faces(grid,llc,faces,rows,cols):

    surface_faces = {}
    for face in range(1,6):
        if face<3:
            surface_faces[face] = np.zeros((3*llc,llc))
        if face==3:
            surface_faces[face] = np.zeros((llc,llc))
        if face>3:
            surface_faces[face] = np.zeros((llc,3*llc))

    for i in range(np.shape(grid)[0]):
        surface_faces[faces[i]][rows[i],cols[i]] = grid[i]

    return(surface_faces)

def stitch_together_arctic_faces(surface_faces,llc,extra_rows = 100):


    stitched_grid = np.nan*np.ones((llc+2*extra_rows,llc+2*extra_rows))

    stitched_grid[extra_rows:-extra_rows,:extra_rows] = np.rot90(surface_faces[1][-extra_rows:,:])
    stitched_grid[:extra_rows, extra_rows:-extra_rows] = surface_faces[2][-extra_rows:, :]
    stitched_grid[extra_rows:-extra_rows,extra_rows:-extra_rows] = surface_faces[3]
    stitched_grid[extra_rows:-extra_rows, -extra_rows:] = surface_faces[4][:, :extra_rows]
    stitched_grid[-extra_rows:, extra_rows:-extra_rows] = np.rot90(surface_faces[5][:, :extra_rows],k=3)

    return(stitched_grid)

def plot_dv_output(var_name,n_timesteps,plot_timestep):

    mask_name = 'mask_arctic_surface'

    llc = 270
    Nr = 50

    if 'plots' not in os.listdir(os.path.join('..')):
        os.mkdir(os.path.join('..','plots'))

    file_path = os.path.join('..','run','dv',mask_name+'_'+var_name+'.bin')
    grid = np.fromfile(file_path,'>f4')

    mask_len = int(np.size(grid)/n_timesteps)

    grid = np.reshape(grid,(n_timesteps,mask_len))

    grid = grid[plot_timestep,:]

    faces, rows, cols = read_reference_dict(mask_name)
    surface_faces = read_grid_to_faces(grid, llc, faces, rows, cols)

    fig = plt.figure(figsize=(15, 4))
    plt.style.use('dark_background')

    arctic_surface = stitch_together_arctic_faces(surface_faces,llc,extra_rows = 100)

    # aspect = mask_len / (3*Nr)
    C = plt.imshow(arctic_surface,origin='lower')#,aspect = aspect)
    plt.colorbar(C)

    # plt.ylabel('Depth Levels')
    # plt.xlabel('Grid Points Along Mask')
    plt.title(var_name+' on mask '+mask_name+' (timestep '+str(plot_timestep+1)+' of '+str(n_timesteps)+')')

    output_file = os.path.join('..', 'plots', mask_name+'_'+var_name+'.png')
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--var_name", action="store",
                        help="Name of the variable to plot.", dest="var_name",
                        type=str, required=True)

    parser.add_argument("-n", "--n_timesteps", action="store", help="The number of timesteps in the output.",
                        dest="n_timesteps", type=int, required=True)

    parser.add_argument("-p", "--plot_timestep", action="store", help="The timestep to plot.",
                        dest="plot_timestep", type=int, required=True)

    args = parser.parse_args()
    var_name = args.var_name
    n_timesteps = args.n_timesteps
    plot_timestep = args.plot_timestep


    plot_dv_output(var_name,n_timesteps,plot_timestep)