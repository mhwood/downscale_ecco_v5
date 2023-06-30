
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

def plot_exf_fields(plot_timestep,nTimesteps):

    if 'plots' not in os.listdir('..'):
        os.mkdir(os.path.join('..','plots'))
    if 'init_files' not in os.listdir(os.path.join('..','plots')):
        os.mkdir(os.path.join('..','plots','init_files'))

    llc = 1080
    n = 720

    n_rows = n*4 + llc
    n_cols = llc

    # these are the list of external forcing fields
    var_names = ['APRESS' ,'AQH' ,'ATEMP' ,'LWDOWN' ,'PRECIP' ,'SWDOWN' ,'USTRESS' ,'VSTRESS' ,'WSPEED','RUNOFF']

    print('Plotting the exf fields at timestep '+str(plot_timestep))

    fig = plt.figure(figsize=(17, 18))

    plt.style.use("dark_background")

    for vn in range(len(var_names)):

        var_name = var_names[vn]
        print('Reading in '+var_name)
        var_grid_file = os.path.join('..','input','exf','N1_exf_'+var_name+'.bin')
        var_grid = np.fromfile(var_grid_file,'>f4')
        var_grid_compact = np.reshape(var_grid,(nTimesteps,n_rows,n_cols))
        var_grid_compact = var_grid_compact[plot_timestep,:,:]
        var_grid_faces = sf.sassie_n1_compact_to_faces(var_grid_compact,llc=llc, n=720)
        var_grid = sf.stitch_faces_to_single_grid(var_grid_faces, llc=llc, n=720)

        if var_name in ['USTRESS' ,'VSTRESS']:
            vmin = -0.1
            vmax = 0.1
            cmap = 'RdBu'
        else:
            cmap='viridis'
            if np.any(var_grid>0):
                vmin = np.min(var_grid[np.logical_and(var_grid != 0,~np.isnan(var_grid))])
                vmax = np.max(var_grid[np.logical_and(var_grid != 0,~np.isnan(var_grid))])
            else:
                vmin = -0.1
                vmax = 0.1

        plt.subplot(4,3,vn+1)
        C = plt.imshow(var_grid, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)
        plt.title(var_name)
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
        plt.colorbar(C)

    plt.suptitle('External Forcing Fields (timestep = '+str(plot_timestep)+')')

    plt.savefig(os.path.join('..', 'plots','init_files','exf_fields_'+str(plot_timestep)+'.png'), bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--nTimesteps", action="store",
                        help="The number of timesteps stored in the exf files.", dest="nTimesteps",
                        type=int, required=True)

    parser.add_argument("-t", "--time_step", action="store",
                        help="Timestep of the output to plot.", dest="time_step",
                        type=int, required=True)

    args = parser.parse_args()
    time_step = args.time_step
    nTimesteps = args.nTimesteps

    plot_exf_fields(time_step, nTimesteps)