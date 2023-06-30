
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
from matplotlib.gridspec import GridSpec
import sassie_functions as sf

def plot_bathymetry_with_BCs():

    config_dir = '../..'
    llc = 1080
    n = 720
    var_name = 'THETA'
    Nr = 90
    n_points = 2*llc+2*680

    if 'plots' not in os.listdir('..'):
        os.mkdir(os.path.join('..','plots'))
    if 'init_files' not in os.listdir(os.path.join('..','plots')):
        os.mkdir(os.path.join('..','plots','init_files'))

    f = open(os.path.join('..', '..', 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L1_size = size_dict['N1_1080']
    n_rows = L1_size[0]
    n_cols = L1_size[1]
    n_rows = 1080+4*720

    bathy_file = os.path.join('..', 'input', 'bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file, '>f4')
    bathy_grid_compact = np.reshape(bathy_grid, (n_rows, n_cols))
    bathy_grid_faces = sf.sassie_n1_compact_to_faces(bathy_grid_compact,llc,n)

    random_timestep = 171
    south_BC_file_path = os.path.join(config_dir,'N1_1080','input','obcs','N1_1080_BC_south_' + var_name + '.bin')
    south_bc_grid = np.fromfile(south_BC_file_path,'>f4')
    n_timesteps = int(np.size(south_bc_grid)/(n_points*Nr))
    south_bc_grid = np.reshape(south_bc_grid,(n_timesteps,Nr,n_points))
    south_bc_grid = south_bc_grid[random_timestep,:,:]

    east_BC_file_path = os.path.join(config_dir, 'N1_1080', 'input', 'obcs',
                                      'N1_1080_BC_east_' + var_name + '.bin')
    east_bc_grid = np.fromfile(east_BC_file_path, '>f4')
    n_timesteps = int(np.size(east_bc_grid) / (n_points * Nr))
    east_bc_grid = np.reshape(east_bc_grid, (n_timesteps, Nr, n_points))
    east_bc_grid = east_bc_grid[random_timestep,:,:]


    print('Plotting the bathymetry')

    fig = plt.figure(figsize=(9,12))
    gs = GridSpec(25, 12, left=0.05, right=0.95, hspace=0.2, wspace = 0.05)

    plt.style.use("dark_background")

    vmin = -10
    vmax = 0

    sNx = 60
    sNy = 60

    tile_counter = 0
    blank_list = []

    for face in range(1,6):
        if face==1:
            ax = fig.add_subplot(gs[10:14, 1:4])
        if face==2:
            ax = fig.add_subplot(gs[10:14, 5:8])
        if face==3:
            ax = fig.add_subplot(gs[5:9, 4:9])
        if face==4:
            ax = fig.add_subplot(gs[5:9, 8:12])
        if face==5:
            ax = fig.add_subplot(gs[:4, 8:12])

        if face in [1,2]:
            C = plt.imshow(bathy_grid_faces[face][40:,:], origin='lower', vmin=vmin, vmax=vmax, cmap='Blues_r')
            plt.title('Face '+str(face))
        elif face in [4,5]:
            C = plt.imshow(bathy_grid_faces[face][:,:-40], origin='lower', vmin=vmin, vmax=vmax, cmap='Blues_r')
            plt.title('Face '+str(face))
        else:
            C = plt.imshow(bathy_grid_faces[face], origin='lower', vmin=vmin, vmax=vmax, cmap='Blues_r')
            plt.title('Face '+str(face))


    ax = fig.add_subplot(gs[15:19, :])
    C = ax.imshow(south_bc_grid,aspect=10)
    plt.colorbar(C, fraction=0.046, pad=0.04)
    plt.title(var_name+' on South Boundary (Randomly Chosen Time Level = '+str(random_timestep)+')')
    plt.plot([1080,1080], [0,89],'w--',linewidth=0.6)
    plt.plot([2*1080, 2*1080], [0, 89], 'w--',linewidth=0.6)
    plt.plot([2*1080+680, 2*1080+680], [0, 89], 'w--',linewidth=0.6)
    ax.set_xticklabels([])
    plt.text(1080/2,5,'Face 1',ha='center',va='top')
    plt.text(1080 + 1080 / 2, 5, 'Face 2', ha='center', va='top')
    plt.text(2*1080 + 680 / 2, 5, 'Face 4', ha='center', va='top')
    plt.text(2*1080 + 680 + 680 / 2, 5, 'Face 5', ha='center', va='top')
    plt.text(2*1080+680,45,'Filled with zeros\nb/c I believe obcs expects\nvalues along these faces\n- not a real boundary',ha='center',va='center')
    plt.ylabel('Depth Cells')

    ax = fig.add_subplot(gs[20:24, :])
    C = ax.imshow(east_bc_grid,aspect=10)
    plt.colorbar(C, fraction=0.046, pad=0.04)
    plt.title(var_name + ' on East Boundary (Randomly Chosen Time Level = ' + str(random_timestep) + ')')
    plt.plot([680, 680], [0, 89], 'w--',linewidth=0.6)
    plt.plot([2 * 680, 2 * 680], [0, 89], 'w--',linewidth=0.6)
    plt.plot([2 * 680 + 1080, 2 * 680 + 1080], [0, 89], 'w--',linewidth=0.6)
    plt.text(680 / 2, 5, 'Face 1', ha='center', va='top')
    plt.text(680 + 680 / 2, 5, 'Face 2', ha='center', va='top')
    plt.text(2* 680 + 1080 / 2, 5, 'Face 4', ha='center', va='top')
    plt.text(2* 680 + 1080 + 1080 / 2, 5, 'Face 5', ha='center', va='top')
    plt.text(680, 45, 'Filled with zeros\nb/c I believe obcs expects\nvalues along these faces\n - not a real boundary', ha='center',va='center')
    plt.ylabel('Depth Cells')
    plt.xlabel('Points Along Boundary')

        # print('Face '+str(face),np.shape(bathy_grid_faces[face]))
        #
        # nPx = int(np.shape(bathy_grid_faces[face])[1] / sNx)
        # nPy = int(np.shape(bathy_grid_faces[face])[0] / sNy)
        #
        # for j in range(nPy):
        #     for i in range(nPx):
        #         tile_counter += 1
        #         bathy_subset = bathy_grid_faces[face][sNy * j:sNy * (j + 1), sNx * i:sNx * (i + 1)]
        #         if np.any(bathy_subset<0) and tile_counter not in manual_blank_list:
        #             rect = Rectangle((sNx*i,sNy*j),sNx,sNy,edgecolor='k',facecolor='none',alpha = 0.3)
        #         else:
        #             rect = Rectangle((sNx * i, sNy * j), sNx, sNy, edgecolor='k', facecolor='none', alpha=0.3, hatch='///')
        #             blank_list.append(tile_counter)
        #
        #         plt.gca().add_patch(rect)
        #         # plt.text(sNx * (i+0.5),sNy*(j+0.5),str(tile_counter),va='center',ha='center',color='k')

    ax = fig.add_subplot(gs[20:24, :])
    plt.text(0,0,'N0_1080 Bathymetry\nTile Size: '+str(sNx)+' x '+str(sNy)+
             '\nTotal Tiles: '+str(tile_counter)+'\nBlank Tiles: '+str(len(blank_list))+'\nRun Tiles: '+str(tile_counter-len(blank_list)),
             color='w',ha='center',va='center',fontsize = 16)
    plt.gca().set_xlim([-1, 1])
    plt.gca().set_ylim([-1, 1])
    plt.axis('off')

    # print('  blankList = ')
    # for tile in blank_list:
    #     print(' '+str(tile)+',')

    output_file_name = 'N1_1080_'+var_name+'_BC_with_bathy_ref.png'

    plt.savefig(os.path.join('..', 'plots','init_files',output_file_name))#, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    plot_bathymetry_with_BCs()