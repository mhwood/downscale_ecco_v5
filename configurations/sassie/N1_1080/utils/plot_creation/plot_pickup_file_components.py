
import os
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import mds
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf


def read_pickup_file_to_faces(input_dir):

    Nr = 90
    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup.0000000001'), returnmeta=True)

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]

        level=0
        var_grid = var_grid[level,:,:]
        var_grid_faces = sf.sassie_n1_compact_to_faces(var_grid,llc=1080,n=720)
        var_grid = sf.stitch_faces_to_single_grid(var_grid_faces,llc=1080,n=720)

        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def plot_pickup_file_components(output_dir,var_names,var_grids):
    fig = plt.figure(figsize=(18, 12))

    plt.style.use("dark_background")

    counter = 1

    for ff in range(len(var_names)):

        plt.subplot(3, 4, counter)

        print(' Plotting '+var_names[ff])#+' (global min: '+str(np.max(var_grids[ff][var_grids[ff]!=0]))+
              # ', min: '+str(np.min(var_grids[ff][var_grids[ff]!=0]))+')')

        var_grid_subset = var_grids[ff]

        cmap = 'viridis'

        if var_names[ff] in ['Uvel', 'Vvel', 'GuNm1', 'GuNm2', 'GvNm1', 'GvNm2']:
            cmap = 'RdBu'
            if var_names[ff] in ['Uvel', 'Vvel']:
                vmin = -0.15
                vmax = 0.15
            else:
                vmin = -1.5e-5
                vmax = 1.5e-5
        else:
            vmin = np.min(var_grid_subset[np.logical_and(var_grids[-3]!= 0,~np.isnan(var_grids[-3]))])
            vmax = np.max(var_grid_subset[np.logical_and(var_grids[-3]!= 0,~np.isnan(var_grids[-3]))])

        C = plt.imshow(var_grid_subset, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C)
        plt.title(var_names[ff])
        counter += 1

    plt.savefig(os.path.join(output_dir,'pickup_file_components.png'),bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs

if 'plots' not in os.listdir('..'):
    os.mkdir(os.path.join('..','plots'))
if 'init_files' not in os.listdir(os.path.join('..','plots')):
    os.mkdir(os.path.join('..','plots','init_files'))

input_dir = os.path.join('..','input')
output_dir = os.path.join('..','plots','init_files')

var_names,row_bounds,var_grids,global_metadata = read_pickup_file_to_faces(input_dir)

plot_pickup_file_components(output_dir,var_names,var_grids)

