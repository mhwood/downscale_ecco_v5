
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse
import ast
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

def plot_bathymetry_with_procs():

    llc = 1080
    n = 680

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
    n_rows = 1080+4*680

    bathy_file = os.path.join('..', 'input', 'bathy_N1_walls')
    bathy_grid = np.fromfile(bathy_file, '>f4')
    bathy_grid_compact = np.reshape(bathy_grid, (n_rows, n_cols))
    bathy_grid_faces = sf.sassie_n1_compact_to_faces(bathy_grid_compact,llc,n)

    mask_name = 'N2_boundary'
    plot_as_points = True

    mask_file = os.path.join('..', 'input','dv', mask_name+'_mask.bin')
    mask_grid = np.fromfile(mask_file, '>f4')
    mask_grid_compact = np.reshape(mask_grid, (n_rows, n_cols))
    mask_grid_faces = sf.sassie_n1_compact_to_faces(mask_grid_compact, llc, n)

    print('Plotting the bathymetry')

    fig = plt.figure(figsize=(10,10))

    plt.style.use("dark_background")

    vmin = -10
    vmax = 0

    sNx = 40
    sNy = 40

    tile_counter = 0
    blank_list = []

    if sNx==180 and sNy==180:
        manual_blank_list = [11,17,18,29,43,74,79,105,120]
    else:
        manual_blank_list = []

    for face in range(1,6):
        if face==1:
            plt.subplot(3,3,7)
        if face==2:
            plt.subplot(3,3,8)
        if face==3:
            plt.subplot(3,3,5)
        if face==4:
            plt.subplot(3,3,6)
        if face==5:
            plt.subplot(3,3,3)
        C = plt.imshow(bathy_grid_faces[face], origin='lower', vmin=vmin, vmax=vmax, cmap='Blues_r')
        if np.any(mask_grid_faces[face])>0:
            if plot_as_points:
                rows,cols = np.where(mask_grid_faces[face]!=0)
                if len(rows)>0:
                    plt.plot(cols,rows,'r.')
            else:
                masked_mask = np.ma.masked_where(mask_grid_faces[face]==0,mask_grid_faces[face])
                C = plt.imshow(masked_mask, origin='lower', vmin=1, vmax=np.max(mask_grid_faces[face]))
        plt.title('Face '+str(face))

        print('Face '+str(face),np.shape(bathy_grid_faces[face]))

        nPx = int(np.shape(bathy_grid_faces[face])[1] / sNx)
        nPy = int(np.shape(bathy_grid_faces[face])[0] / sNy)

        for j in range(nPy):
            for i in range(nPx):
                tile_counter += 1
                bathy_subset = bathy_grid_faces[face][sNy * j:sNy * (j + 1), sNx * i:sNx * (i + 1)]
                if np.any(bathy_subset<0) and tile_counter not in manual_blank_list:
                    rect = Rectangle((sNx*i,sNy*j),sNx,sNy,edgecolor='k',facecolor='none',alpha = 0.3)
                else:
                    rect = Rectangle((sNx * i, sNy * j), sNx, sNy, edgecolor='k', facecolor='none', alpha=0.3, hatch='///')
                    blank_list.append(tile_counter)

                plt.gca().add_patch(rect)
                # plt.text(sNx * (i+0.5),sNy*(j+0.5),str(tile_counter),va='center',ha='center',color='k')

    plt.subplot(3,3,1)
    plt.text(0,0,'N0_1080 Bathymetry\nTile Size: '+str(sNx)+' x '+str(sNy)+
             '\nTotal Tiles: '+str(tile_counter)+'\nBlank Tiles: '+str(len(blank_list))+'\nRun Tiles: '+str(tile_counter-len(blank_list)),
             color='w',ha='center',va='center',fontsize = 16)
    plt.gca().set_xlim([-1, 1])
    plt.gca().set_ylim([-1, 1])
    plt.axis('off')

    # print('  blankList = ')
    # for tile in blank_list:
    #     print(' '+str(tile)+',')

    output_file_name = mask_name+'_mask.png'

    plt.savefig(os.path.join('..', 'plots','init_files',output_file_name), bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    plot_bathymetry_with_procs()