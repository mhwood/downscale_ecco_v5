
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse
import ast
import sys

from matplotlib.gridspec import GridSpec


def plot_test(config_dir):

    diag_file = os.path.join(config_dir,'N1_1080', 'run', 'T.0011577600.data')

    diag_grid = np.fromfile(diag_file, '>f4')
    diag_grid = np.reshape(diag_grid, (1, 90, 102600, 40))
    # diag_grid = np.reshape(diag_grid, (1, 1, 102600, 40))

    diag_grid = diag_grid[0, :, :, :]

    face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
    face_1 = face_1.reshape((90, 680, 1080))

    face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
    face_2 = face_2.reshape((90, 680, 1080))

    face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
    face_3 = face_3.reshape((90, 1080, 1080))

    face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
    face_4 = face_4.reshape((90, 1080, 680))

    face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
    face_5 = face_5.reshape((90, 1080, 680))

    plt.subplot(3,3,7)
    face_to_plot = face_1[0,:,:]
    face_to_plot = np.ma.masked_where(face_to_plot==0,face_to_plot)
    C = plt.imshow(face_to_plot)
    plt.colorbar(C)

    plt.subplot(3, 3, 8)
    face_to_plot = face_2[0, :, :]
    face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    C = plt.imshow(face_to_plot)
    plt.colorbar(C)

    plt.subplot(3, 3, 5)
    face_to_plot = face_3[0, :, :]
    face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    C = plt.imshow(face_to_plot)
    plt.colorbar(C)

    plt.subplot(3, 3, 6)
    face_to_plot = face_4[0, :, :]
    face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    C = plt.imshow(face_to_plot)
    plt.colorbar(C)

    plt.subplot(3, 3, 3)
    face_to_plot = face_5[0, :, :]
    face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    C = plt.imshow(face_to_plot)
    plt.colorbar(C)
    plt.show()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_test(config_dir)