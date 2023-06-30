import os
import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cm
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf


def plot_3D_boundary_conditions(output_dir, input_dir, field_name, time_step, Nr, n, llc):

    N = 3*llc+2*n

    # read in the field grids
    south_file_name = os.path.join(input_dir,'obcs', 'N1_1080_BC_south_' + field_name + '.bin')
    south_field_grid = np.fromfile(south_file_name, '>f4')
    n_timesteps = int(np.size(south_field_grid)/(Nr*N))
    south_field_grid = np.reshape(south_field_grid, (n_timesteps, Nr, N))
    south_field_grid = south_field_grid[time_step, :, :]
    # south_field_grid = south_field_grid[:,:2*llc]

    east_file_name = os.path.join(input_dir,'obcs', 'N1_1080_BC_east_' + field_name + '.bin')
    east_field_grid = np.fromfile(east_file_name, '>f4')
    n_timesteps = int(np.size(east_field_grid) / (Nr * N))
    east_field_grid = np.reshape(east_field_grid, (n_timesteps, Nr, N))
    east_field_grid = east_field_grid[time_step, :, :]
    # east_field_grid = east_field_grid[:, -2 * llc:]

    n_cols = llc
    n_rows = llc + 4 * n
    bathy_file = os.path.join(input_dir, 'bathy_N1_no_walls')
    bathy_grid = np.fromfile(bathy_file, '>f4')
    bathy_grid_compact = np.reshape(bathy_grid, (n_rows, n_cols))
    bathy_grid_faces = sf.sassie_n1_compact_to_faces(bathy_grid_compact, llc, n)

    south_field_grid = np.ma.masked_where(south_field_grid==0, south_field_grid)
    east_field_grid = np.ma.masked_where(east_field_grid == 0, east_field_grid)

    fig = plt.figure(figsize=(12, 7))

    plt.style.use("dark_background")

    vmin = 0
    vmax = 0

    if field_name=='THETA':
        vmin = 0
        vmax = 17
        cmap = cm.thermal

    if field_name=='SALT':
        vmin = 32.5
        vmax = 36
        cmap = cm.haline

    if field_name in ['UVEL', 'VVEL']:
        vmin = -0.5
        vmax = 0.5
        cmap='seismic'

    # on the actual boundary
    plt.subplot(3,1,1)
    aspect = np.shape(south_field_grid)[1] / np.shape(south_field_grid)[0]
    aspect=10
    C = plt.imshow(south_field_grid, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax)
    plt.plot([1080,1080],[0,89],'w--',linewidth=0.6)
    plt.plot([2*1080, 2*1080], [0, 89], 'w--', linewidth=0.6)
    plt.plot([3*1080, 3*1080], [0, 89], 'w--', linewidth=0.6)
    plt.plot([3*1080+680, 3*1080+680], [0, 89], 'w--', linewidth=0.6)
    plt.text(1080 / 2, 5, 'Face 1', ha='center', va='top')
    plt.text(1080 + 1080 / 2, 5, 'Face 2', ha='center', va='top')
    plt.text(2*1080 + 1080 / 2, 5, 'Face 3', ha='center', va='top')
    plt.text(3*1080 + 680 / 2, 5, 'Face 4', ha='center', va='top')
    plt.text(3*1080+680 + 680 / 2, 5, 'Face 5', ha='center', va='top')
    # plt.text(1080, 85, 'South Boundary', ha='center', va='bottom')
    plt.gca().set_xticks([])
    plt.colorbar(C)
    # plt.xlabel('points along boundary')
    plt.ylabel('south boundary\ndepth cells')
    plt.title('BCs for ' + field_name + ' (time step ' + str(time_step)+')')

    plt.subplot(3, 1, 2)
    aspect = np.shape(east_field_grid)[1] / np.shape(east_field_grid)[0]
    aspect = 10
    C = plt.imshow(east_field_grid, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax)
    plt.plot([680, 680], [0, 89], 'w--',linewidth=0.6)
    plt.plot([2 * 680, 2 * 680], [0, 89], 'w--', linewidth=0.6)
    plt.plot([2 * 680 + 1080, 2 * 680+1080], [0, 89], 'w--', linewidth=0.6)
    plt.plot([2 * 680 + 2*1080, 2 * 1080 + 2*680], [0, 89], 'w--', linewidth=0.6)
    plt.text(680 / 2, 5, 'Face 1', ha='center', va='top')
    plt.text(680 + 680 / 2, 5, 'Face 2', ha='center', va='top')
    plt.text(2 * 680 + 1080 / 2, 5, 'Face 3', ha='center', va='top')
    plt.text(2 * 680 + 1080 + 1080 / 2, 5, 'Face 4', ha='center', va='top')
    plt.text(2 * 680 + 2*1080 + 1080 / 2, 5, 'Face 5', ha='center', va='top')
    # plt.text(1400, 85, 'East Boundary', ha='center', va='bottom')
    plt.colorbar(C)
    plt.xlabel('points along boundary')
    plt.ylabel('east boundary\ndepth cells')
    # plt.title('Boundary Conditions for ' + field_name + ' at time step ' + str(time_step))

    plt.subplot(3,1,3)
    plt.plot(np.arange(1080), bathy_grid_faces[1][0,:], 'w-')
    plt.plot(np.arange(1080,2*1080), bathy_grid_faces[2][0, :], 'w-')
    plt.plot(np.arange(2*1080,3*1080), bathy_grid_faces[4][:,-1], 'w-')
    plt.plot(np.arange(3*1080, 4*1080), bathy_grid_faces[5][:,-1], 'w-')

    plt.savefig(os.path.join(output_dir, 'bcs_' + field_name + '_' + str(time_step) + '.png'),bbox_inches='tight')
    plt.close(fig)


def plot_boundary_conditions(field_name,time_step):
    print('Plotting the boundary conditions for '+field_name+' at time step '+str(time_step))

    input_dir = os.path.join('..', 'input')
    output_dir = os.path.join('..','plots','init_files')

    if 'plots' not in os.listdir('..'):
        os.mkdir(os.path.join('..','plots'))
    if 'init_files' not in os.listdir(os.path.join('..','plots')):
        os.mkdir(os.path.join('..','plots','init_files'))

    Nr = 90
    llc = 1080
    n = 680

    plot_3D_boundary_conditions(output_dir, input_dir, field_name, time_step, Nr, n, llc)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--field_name", action="store",
                        help="Name of the boundary field to plot (THETA, SALT).", dest="field_name",
                        type=str, required=True)

    parser.add_argument("-t", "--time_step", action="store",
                        help="Timestep at which to plot the boundary field.", dest="time_step",
                        type=int, required=True)

    args = parser.parse_args()
    field_name = args.field_name
    time_step = args.time_step

    plot_boundary_conditions(field_name,time_step)