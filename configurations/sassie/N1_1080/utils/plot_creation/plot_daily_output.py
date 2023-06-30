
import os
import time
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import cartopy.feature as cf
import simplegrid as sg
import matplotlib.path as mpath

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean.cm as cm
import argparse
from matplotlib.patches import Rectangle


def field_number_to_set_and_name(field_number,dim):

    # this is a dictionary referencing whats in the diagnostics file

    #2D files
    field_number_dict_2D = {1: ['ocean_state_2D_day_mean','ETAN',0, False],
                            2: ['ocean_state_2D_day_mean', 'PHIBOT', 1, False],
                            3: ['ocean_state_2D_day_mean', 'sIceLoad', 2, False],
                            4: ['seaice_state_day_mean','SIarea',0, False],
                            5: ['seaice_state_day_mean','SIheff',1, False],
                            6: ['seaice_state_day_mean','SIhsnow',2, False],
                            7: ['tr_adv_x_2D_day_mean', 'ADVxHEFF', 0, True],
                            8: ['tr_adv_x_2D_day_mean', 'ADVxSNOW', 1, True],
                            9: ['tr_adv_x_2D_day_mean', 'ADVyHEFF', 2, True],
                            10: ['tr_adv_x_2D_day_mean', 'ADVySNOW', 3, True],
                            11: ['EXF_day_mean', 'EXFaqh', 0, False],
                            12: ['EXF_day_mean', 'EXFatemp', 1, False],
                            13: ['EXF_day_mean', 'EXFempmr', 2, False],
                            14: ['EXF_day_mean', 'EXFevap', 3, False],
                            15: ['EXF_day_mean', 'EXFpreci', 4, False],
                            16: ['EXF_day_mean', 'EXFroff', 5, False],
                            17: ['EXF_day_mean', 'EXFqnet', 6, False],
                            18: ['EXF_day_mean', 'EXFhl', 7, False],
                            19: ['EXF_day_mean', 'EXFhs', 8, False],
                            20: ['EXF_day_mean', 'EXFlwdn', 9, False],
                            21: ['EXF_day_mean', 'EXFlwnet', 10, False],
                            22: ['EXF_day_mean', 'EXFswdn', 11, False],
                            23: ['EXF_day_mean', 'EXFswnet', 12, False],
                            24: ['EXF_day_mean', 'EXFuwind', 13, True],
                            25: ['EXF_day_mean', 'EXFtaux', 14, True],
                            26: ['EXF_day_mean', 'EXFvwind', 15, True],
                            27: ['EXF_day_mean', 'EXFtauy', 16, True],
                            28: ['oce_flux_day_mean', 'oceFWflx', 0, False],
                            29: ['oce_flux_day_mean', 'oceQnet', 1, False],
                            30: ['oce_flux_day_mean', 'oceQsw', 2, False],
                            31: ['oce_flux_day_mean', 'oceTAUX', 3, True],
                            32: ['oce_flux_day_mean', 'oceTAUY', 4, True],
                            33: ['oce_flux_day_mean', 'SFLUX', 5, False],
                            34: ['oce_flux_day_mean', 'TFLUX', 6, False],
                            35: ['seaice_flux_day_mean', 'SIatmFW', 0, False],
                            36: ['seaice_flux_day_mean', 'SIatmQnt', 1, False],
                            37: ['seaice_vel_day_mean', 'SIuice', 0, True],
                            38: ['seaice_vel_day_mean', 'SIvice', 1, True],
                            39: ['KPP_hbl_day_mean', 'KPPhbl', 0, False]}

    # 3D files
    field_number_dict_3D = {1: ['ocean_state_3D_day_mean', 'SALT', 0, False], # 1
                            2: ['ocean_state_3D_day_mean', 'THETA', 1, False],
                            3: ['tr_adv_r_day_mean', 'ADVr_SLT', 0, False],
                            4: ['tr_adv_r_day_mean', 'ADVr_TH', 1, False],
                            5: ['tr_adv_x_3D_day_mean', 'ADVx_SLT', 0, True],#2
                            6: ['tr_adv_x_3D_day_mean', 'ADVx_TH', 1, True],
                            7: ['tr_adv_x_3D_day_mean', 'ADVy_SLT', 2, True],
                            8: ['tr_adv_x_3D_day_mean', 'ADVy_TH', 3, True],
                            9: ['tr_diff_r_day_mean', 'DFrE_SLT', 0, False], # 3
                            10: ['tr_diff_r_day_mean', 'DFrE_TH', 1, False],
                            11: ['tr_diff_r_day_mean', 'DFrI_SLT', 2, False],
                            12: ['tr_diff_r_day_mean', 'DFrI_TH', 3, False],
                            13: ['tr_diff_x_day_mean', 'DFxE_SLT', 0, True], # 4
                            14: ['tr_diff_x_day_mean', 'DFxE_TH', 1, True],
                            15: ['tr_diff_x_day_mean', 'DFyE_SLT', 2, True],
                            16: ['tr_diff_x_day_mean', 'DFyE_TH', 3, True],
                            17: ['ocean_vel_day_mean', 'UVEL', 0, True], # 5
                            18: ['ocean_vel_day_mean', 'VVEL', 1, True],
                            19: ['ocean_vel_day_mean', 'WVEL', 2, False],
                            20: ['vol_adv_day_mean', 'UVELMASS', 0, True],
                            21: ['vol_adv_day_mean', 'VVELMASS', 1, True], # 6
                            22: ['vol_adv_day_mean', 'WVELMASS', 2, False],
                            23: ['phi_3D_day_mean', 'PHIHYD', 0, False],
                            24: ['phi_3D_day_mean', 'PHIHYDcR', 1, False],
                            25: ['phi_3D_day_mean', 'RHOAnoma', 2, False], # 7
                            26: ['KPP_mix_day_mean', 'KPPdiffS', 0, False],
                            #27: ['KPP_mix_day_mean', 'KPPdiffT', 1, False],
                            27: ['KPP_mix_day_mean', 'KPPviscA', 1, False]}

    if dim==2:
        field_set = field_number_dict_2D[field_number][0]
        field_name = field_number_dict_2D[field_number][1]
        field_set_index = field_number_dict_2D[field_number][2]
        is_vector = field_number_dict_2D[field_number][3]
    if dim==3:
        field_set = field_number_dict_3D[field_number][0]
        field_name = field_number_dict_3D[field_number][1]
        field_set_index = field_number_dict_3D[field_number][2]
        is_vector = field_number_dict_3D[field_number][3]

    return(field_set,field_name,field_set_index,is_vector)

def get_filename_date_pairs(config_dir, field_set):

    def file_iter_to_date(iter_number):
        timestep = 120
        seconds_elapsed = iter_number*timestep
        date = datetime(1992,1,1) + timedelta(seconds=seconds_elapsed)
        return(date.year, date.month, date.day)

    diag_dir = os.path.join(config_dir,'N1_1080','run','diags',field_set)
    all_file_names = os.listdir(diag_dir)

    file_prefixes = []
    dates = []

    for file_name in all_file_names:
        if file_name[0]!='.' and file_name[-5:]=='.data':
            file_prefix = file_name[:-5]
            file_prefixes.append(file_prefix)
            iter_number = int(file_prefix.split('.')[-1])
            year, month, day = file_iter_to_date(iter_number)
            dates.append(str(year)+'{:02d}'.format(month)+'{:02d}'.format(day))

    return(file_prefixes, dates)

def read_faces_geometry(config_dir, is_vector):
    llc = 1080
    N = 680

    grid_file_dir = os.path.join(config_dir,'N1_1080','input')
    XC_faces = {}
    YC_faces = {}
    for i in [1,2,3,4,5]:
        if i<3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, N)
        if i==3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, llc)
        if i>3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), N, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    AngleCS_faces = {}
    AngleSN_faces = {}
    if is_vector:
        AngleCS_grid = np.fromfile(os.path.join(grid_file_dir, 'AngleCS.data'), '>f4')
        N = int(np.size(AngleCS_grid) / (102600 * 40))
        AngleCS_grid = np.reshape(AngleCS_grid, (N, 102600, 40))

        face_1 = AngleCS_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1 = face_1[0, :, :]
        AngleCS_faces[1] = face_1

        face_2 = AngleCS_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2 = face_2[0, :, :]
        AngleCS_faces[2] = face_2

        face_3 = AngleCS_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3 = face_3[0, :, :]
        AngleCS_faces[3] = face_3

        face_4 = AngleCS_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4 = face_4[0, :, :]
        AngleCS_faces[4] = face_4

        face_5 = AngleCS_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5 = face_5[0, :, :]
        AngleCS_faces[5] = face_5

        AngleSN_grid = np.fromfile(os.path.join(grid_file_dir, 'AngleSN.data'), '>f4')
        N = int(np.size(AngleSN_grid) / (102600 * 40))
        AngleSN_grid = np.reshape(AngleSN_grid, (N, 102600, 40))

        face_1 = AngleSN_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1 = face_1[0, :, :]
        AngleSN_faces[1] = face_1

        face_2 = AngleSN_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2 = face_2[0, :, :]
        AngleSN_faces[2] = face_2

        face_3 = AngleSN_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3 = face_3[0, :, :]
        AngleSN_faces[3] = face_3

        face_4 = AngleSN_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4 = face_4[0, :, :]
        AngleSN_faces[4] = face_4

        face_5 = AngleSN_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5 = face_5[0, :, :]
        AngleSN_faces[5] = face_5


    return(XC_faces,YC_faces,AngleCS_faces,AngleSN_faces)

def read_field_to_faces(config_dir, dim, field_set, file_name, field_set_index):

    file_path = os.path.join(config_dir,'N1_1080','run','diags',field_set,file_name+'.data')

    if dim==2:
        diag_grid = np.fromfile(file_path, '>f4')
        N = int(np.size(diag_grid) / (102600 * 40))
        diag_grid = np.reshape(diag_grid, (N, 102600, 40))

        face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1 = face_1[field_set_index, :, :]

        face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2 = face_2[field_set_index, :, :]

        face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3 = face_3[field_set_index, :, :]

        face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4 = face_4[field_set_index, :, :]

        face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5 = face_5[field_set_index, :, :]

    if dim==3:
        depth_level = 10
        Nr = 90

        diag_grid = np.fromfile(file_path, '>f4')
        N = int(np.size(diag_grid) / (102600 * 40 * Nr))
        diag_grid = np.reshape(diag_grid, (N * Nr, 102600, 40))

        face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N* Nr, 680, 1080))
        face_1 = face_1[field_set_index*Nr + depth_level, :, :]

        face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N * Nr, 680, 1080))
        face_2 = face_2[field_set_index*Nr + depth_level, :, :]

        face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N * Nr, 1080, 1080))
        face_3 = face_3[field_set_index*Nr + depth_level, :, :]

        face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N * Nr, 1080, 680))
        face_4 = face_4[field_set_index*Nr + depth_level, :, :]

        face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N * Nr, 1080, 680))
        face_5 = face_5[field_set_index*Nr + depth_level, :, :]

    # plt.subplot(3,3,7)
    # face_to_plot = face_1
    # face_to_plot = np.ma.masked_where(face_to_plot==0,face_to_plot)
    # C = plt.imshow(face_to_plot)
    # plt.colorbar(C)
    #
    # plt.subplot(3, 3, 8)
    # face_to_plot = face_2
    # face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    # C = plt.imshow(face_to_plot)
    # plt.colorbar(C)
    #
    # plt.subplot(3, 3, 5)
    # face_to_plot = face_3
    # face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    # C = plt.imshow(face_to_plot)
    # plt.colorbar(C)
    #
    # plt.subplot(3, 3, 6)
    # face_to_plot = face_4
    # face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    # C = plt.imshow(face_to_plot)
    # plt.colorbar(C)
    #
    # plt.subplot(3, 3, 3)
    # face_to_plot = face_5
    # face_to_plot = np.ma.masked_where(face_to_plot == 0, face_to_plot)
    # C = plt.imshow(face_to_plot)
    # plt.colorbar(C)
    # plt.show()

    del diag_grid

    field_faces = {1: face_1, 2: face_2, 3: face_3, 4: face_4, 5: face_5}

    return (field_faces)

def read_vector_field_to_faces(config_dir, dim, field_set, file_name, field_set_index, AngleCS_faces, AngleSN_faces):

    def vector_companion_indices(field_set, field_set_index):
        if field_set == 'tr_adv_x_2D_day_mean':
            if field_set_index == 0:
                u_index = field_set_index
                v_index = 2
                output_type = 'u'
            if field_set_index == 2:
                v_index = field_set_index
                u_index = 0
                output_type = 'v'
            if field_set_index == 1:
                u_index = field_set_index
                v_index = 3
                output_type = 'u'
            if field_set_index == 3:
                v_index = field_set_index
                u_index = 1
                output_type = 'v'
        if field_set == 'EXF_day_mean':
            if field_set_index == 13:
                u_index = field_set_index
                v_index = 15
                output_type = 'u'
            if field_set_index == 15:
                v_index = field_set_index
                u_index = 13
                output_type = 'v'
            if field_set_index == 14:
                u_index = field_set_index
                v_index = 16
                output_type = 'u'
            if field_set_index == 16:
                v_index = field_set_index
                u_index = 14
                output_type = 'v'
        if field_set == 'oce_flux_day_mean':
            if field_set_index == 3:
                u_index = field_set_index
                v_index = 4
                output_type = 'u'
            if field_set_index == 4:
                v_index = field_set_index
                u_index = 3
                output_type = 'v'
        if field_set == 'seaice_vel_day_mean':
            if field_set_index == 0:
                u_index = field_set_index
                v_index = 1
                output_type = 'u'
            if field_set_index == 1:
                v_index = field_set_index
                u_index = 0
                output_type = 'v'

        ################################################################################
        # 3D vars

        if field_set == 'tr_adv_x_3D_day_mean':
            if field_set_index == 0:
                u_index = field_set_index
                v_index = 2
                output_type = 'u'
            if field_set_index == 2:
                v_index = field_set_index
                u_index = 0
                output_type = 'v'
            if field_set_index == 1:
                u_index = field_set_index
                v_index = 3
                output_type = 'u'
            if field_set_index == 3:
                v_index = field_set_index
                u_index = 1
                output_type = 'v'

        if field_set == 'tr_diff_x_day_mean':
            if field_set_index == 0:
                u_index = field_set_index
                v_index = 2
                output_type = 'u'
            if field_set_index == 2:
                v_index = field_set_index
                u_index = 0
                output_type = 'v'
            if field_set_index == 1:
                u_index = field_set_index
                v_index = 3
                output_type = 'u'
            if field_set_index == 3:
                v_index = field_set_index
                u_index = 1
                output_type = 'v'

        if field_set == 'ocean_vel_day_mean':
            if field_set_index == 0:
                u_index = field_set_index
                v_index = 1
                output_type = 'u'
            if field_set_index == 1:
                v_index = field_set_index
                u_index = 0
                output_type = 'v'

        if field_set == 'vol_adv_day_mean':
            if field_set_index == 0:
                u_index = field_set_index
                v_index = 1
                output_type = 'u'
            if field_set_index == 1:
                v_index = field_set_index
                u_index = 0
                output_type = 'v'

        return(u_index,v_index,output_type)

    file_path = os.path.join(config_dir,'N1_1080','run','diags',field_set,file_name+'.data')

    u_index, v_index, output_type = vector_companion_indices(field_set, field_set_index)

    if dim==2:
        diag_grid = np.fromfile(file_path, '>f4')
        N = int(np.size(diag_grid) / (102600 * 40))
        diag_grid = np.reshape(diag_grid, (N, 102600, 40))

        face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N, 680, 1080))
        face_1_u = face_1[u_index, :, :]
        face_1_v = face_1[v_index, :, :]
        if output_type=='u':
            face_1 = AngleCS_faces[1] * face_1_u - AngleSN_faces[1] * face_1_v
        else:
            face_1 = AngleSN_faces[1] * face_1_u + AngleCS_faces[1] * face_1_v

        face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N, 680, 1080))
        face_2_u = face_2[u_index, :, :]
        face_2_v = face_2[v_index, :, :]
        if output_type == 'u':
            face_2 = AngleCS_faces[2] * face_2_u - AngleSN_faces[2] * face_2_v
        else:
            face_2 = AngleSN_faces[2] * face_2_u + AngleCS_faces[2] * face_2_v

        face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N, 1080, 1080))
        face_3_u = face_3[u_index, :, :]
        face_3_v = face_3[v_index, :, :]
        if output_type == 'u':
            face_3 = AngleCS_faces[3] * face_3_u - AngleSN_faces[3] * face_3_v
        else:
            face_3 = AngleSN_faces[3] * face_3_u + AngleCS_faces[3] * face_3_v

        face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N, 1080, 680))
        face_4_u = face_4[u_index, :, :]
        face_4_v = face_4[v_index, :, :]
        if output_type == 'u':
            face_4 = AngleCS_faces[4] * face_4_u - AngleSN_faces[4] * face_4_v
        else:
            face_4 = AngleSN_faces[4] * face_4_u + AngleCS_faces[4] * face_4_v

        face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N, 1080, 680))
        face_5_u = face_5[u_index, :, :]
        face_5_v = face_5[v_index, :, :]
        if output_type == 'u':
            face_5 = AngleCS_faces[5] * face_5_u - AngleSN_faces[5] * face_5_v
        else:
            face_5 = AngleSN_faces[5] * face_5_u + AngleCS_faces[5] * face_5_v

    if dim==3:
        depth_level = 10
        Nr = 90

        diag_grid = np.fromfile(file_path, '>f4')
        N = int(np.size(diag_grid) / (102600 * 40 * Nr))
        diag_grid = np.reshape(diag_grid, (N * Nr, 102600, 40))

        face_1 = diag_grid[:, :int(1080 * 680 / 40), :]
        face_1 = face_1.reshape((N* Nr, 680, 1080))
        face_1_u = face_1[u_index * Nr + depth_level, :, :]
        face_1_v = face_1[v_index * Nr + depth_level, :, :]
        if output_type == 'u':
            face_1 = AngleCS_faces[1] * face_1_u - AngleSN_faces[1] * face_1_v
        else:
            face_1 = AngleSN_faces[1] * face_1_u + AngleCS_faces[1] * face_1_v

        face_2 = diag_grid[:, int(1080 * 680 / 40):int(2 * 1080 * 680 / 40), :]
        face_2 = face_2.reshape((N * Nr, 680, 1080))
        face_2_u = face_2[u_index * Nr + depth_level, :, :]
        face_2_v = face_2[v_index * Nr + depth_level, :, :]
        if output_type == 'u':
            face_2 = AngleCS_faces[2] * face_2_u - AngleSN_faces[2] * face_2_v
        else:
            face_2 = AngleSN_faces[2] * face_2_u + AngleCS_faces[2] * face_2_v

        face_3 = diag_grid[:, int(2 * 1080 * 680 / 40):int((1080 * 1080 + 2 * 1080 * 680) / 40), :]
        face_3 = face_3.reshape((N * Nr, 1080, 1080))
        face_3_u = face_3[u_index * Nr + depth_level, :, :]
        face_3_v = face_3[v_index * Nr + depth_level, :, :]
        if output_type == 'u':
            face_3 = AngleCS_faces[3] * face_3_u - AngleSN_faces[3] * face_3_v
        else:
            face_3 = AngleSN_faces[3] * face_3_u + AngleCS_faces[3] * face_3_v

        face_4 = diag_grid[:, int((1080 * 1080 + 2 * 1080 * 680) / 40):int((1080 * 1080 + 3 * 1080 * 680) / 40), :]
        face_4 = face_4.reshape((N * Nr, 1080, 680))
        face_4_u = face_4[u_index * Nr + depth_level, :, :]
        face_4_v = face_4[v_index * Nr + depth_level, :, :]
        if output_type == 'u':
            face_4 = AngleCS_faces[4] * face_4_u - AngleSN_faces[4] * face_4_v
        else:
            face_4 = AngleSN_faces[4] * face_4_u + AngleCS_faces[4] * face_4_v

        face_5 = diag_grid[:, int((1080 * 1080 + 3 * 1080 * 680) / 40):int((1080 * 1080 + 4 * 1080 * 680) / 40), :]
        face_5 = face_5.reshape((N * Nr, 1080, 680))
        face_5_u = face_5[u_index*Nr + depth_level, :, :]
        face_5_v = face_5[v_index*Nr + depth_level, :, :]
        if output_type == 'u':
            face_5 = AngleCS_faces[5] * face_5_u - AngleSN_faces[5] * face_5_v
        else:
            face_5 = AngleSN_faces[5] * face_5_u + AngleCS_faces[5] * face_5_v

    del diag_grid

    field_faces = {1: face_1, 2: face_2, 3: face_3, 4: face_4, 5: face_5}

    return (field_faces)

def field_name_to_plot_metadata(field_name, dim):
    metadata_dict_2D = {'ETAN':[-4,1,'viridis','m','Surface Height Anomaly'],
                        'PHIBOT':[-200,600,cm.haline,'m$^2$/s$^2$','Bottom Pressure Potential Anomaly'],
                        'sIceLoad':[0,5000,cm.ice,'kg/m$^2$','Sea Ice Loading'],
                        'SIarea':[0,1,cm.ice,'m$^2$/m$^2$','Sea Ice Area'],
                        'SIheff':[0,3,cm.ice,'m','Effective Sea Ice Thickness'],
                        'SIhsnow':[0,1,cm.ice,'m','Effective Snow Thickness on Sea Ice'],
                        'ADVxHEFF':[-2000,2000,cm.balance,'m$^2\\cdot$(m/s)','Zonal Advective Flux of Effective Sea Ice Thickness'],
                        'ADVxSNOW':[-400,400,cm.balance,'m$^2\\cdot$(m/s)','Zonal Advective Flux of Effective Snow Thickness on Sea Ice'],
                        'ADVyHEFF':[-2000,2000,cm.balance,'m$^2\\cdot$(m/s)','Meridional Advective Flux of Effective Sea Ice Thickness'],
                        'ADVySNOW':[-400,400,cm.balance,'m$^2\\cdot$(m/s)','Meridional Advective Flux of Effective Snow Thickness on Sea Ice'],
                        'EXFaqh':[0,0.01,cm.tempo,'kg/kg','Surface (2m) Specific Humidity'],
                        'EXFatemp': [250, 290, cm.matter_r, 'Kelvin', 'Surface (2m) Air Temperature'],
                        'EXFempmr': [-1e-7,1e-7, cm.curl, 'm/s', 'Net Upward Freshwater Flux'],
                        'EXFevap': [-5e-8,5e-8, cm.curl, 'm/s', 'Evaporation'],
                        'EXFpreci': [0, 5e-7, cm.rain, 'm/s', 'Precipitation'],
                        'EXFroff': [0, 1e-7, cm.rain, 'm/s', 'River Runoff'],
                        'EXFqnet': [-300, 300, cm.curl, 'W/m$^2$', 'Net Upward Heat Flux'],
                        'EXFhl': [-200, 200, cm.curl, 'W/m$^2$', 'Latent Heat Flux Into the Ocean'],
                        'EXFhs': [-150, 150, cm.curl, 'W/m$^2$', 'Sensible Heat Flux Into the Ocean'],
                        'EXFlwdn': [150, 350, cm.solar, 'W/m$^2$', 'Downward Longwave Radiation'],
                        'EXFlwnet': [-300, 300, cm.curl, 'W/m$^2$', 'Net Upward Longwave Radiation'],
                        'EXFswdn': [30, 300, cm.solar, 'W/m$^2$', 'Downward Shotwave Radiation'],
                        'EXFswnet': [-300, 300, cm.curl, 'W/m$^2$', 'Net Upward Shortwave Radiation'],
                        'EXFuwind': [-15, 15, cm.balance, 'm/s', 'Zonal 10m Wind Speed'],
                        'EXFtaux': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Zonal Surface Wind Stress'],
                        'EXFvwind': [-15, 15, cm.balance, 'm/s', 'Meridional 10m Wind Speed'],
                        'EXFtauy': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Meridional Surface Wind Stress'],
                        'oceFWflx': [-0.0005, 0.0005, cm.diff, 'kg/m$^2$/s', 'Net Surface Freshwater Flux Into the Ocean'],
                        'oceQnet': [-2000, 2000, cm.curl, 'W/m$^2$', 'Net Surface Heat Flux Into the Ocean'],
                        'oceQsw': [-300, 300, cm.curl, 'W/m$^2$', 'Net Shortwave Radiation Into the Ocean'],
                        'oceTAUX': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Zonal Surface Wind Stress'],
                        'oceTAUY': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Meridional Surface Wind Stress'],
                        'SFLUX': [-0.02, 0.02, cm.tarn, 'g/m^2/s', 'Total Salt Flux'],
                        'TFLUX': [-500,500, cm.curl, 'W/m$^2$', 'Total Heat Flux'],
                        'SIatmFW': [-0.0001, 0.0001, cm.diff_r, 'kg/m^2/s', 'Net Freshwater Flux from Atmosphere and Land'],
                        'SIatmQnt': [-300, 300, cm.curl, 'W/m$^2$', 'Net Atmospheric Heat Flux'],
                        'SIuice': [-1, 1, cm.balance, 'm/s', 'Zonal Sea Ice Velocity'],
                        'SIvice': [-1, 1, cm.balance, 'm/s', 'Meridional Sea Ice Velocity'],
                        'KPPhbl': [0,500, cm.deep, 'm', 'KPP Boundary Layer Depth'],
                        }

    metadata_dict_3D = {'THETA': [-2, 12, cm.thermal, '$^{\circ}$C', 'Potential Temperature'], #
                        'SALT': [25, 35, cm.haline, 'psu', 'Practical Salinity'], #
                        'ADVr_SLT': [-1e5, 1e5, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Vertical Advective Flux of Salinity'], #
                        'ADVr_TH': [-25000,25000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Vertical Advective Flux of Potential Temperature'], #
                        'ADVx_SLT': [-1e6, 1e6, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Zonal Advective Flux of Salinity'],
                        'ADVx_TH': [-50000, 50000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Zonal Advective Flux of Potential Temperature'],
                        'ADVy_SLT': [-1e6, 1e6, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Meridional Advective Flux of Salinity'],
                        'ADVy_TH': [-50000, 50000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Meridional Advective Flux of Potential Temperature'],
                        'DFrE_SLT': [-1, 1, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Salinity (Explicit Part)'],
                        'DFrE_TH': [-1, 1, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Potential Temperature (Explicit Part)'],
                        'DFrI_SLT': [-1000,1000, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Salinity (Implicit Part)'],
                        'DFrI_TH': [-2000,2000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Potential Temperature (Implicit Part)'],
                        'DFxE_SLT': [-1, 1, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Zonal Diffusive Flux of Salinity'],
                        'DFxE_TH': [-1, 1, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Zonal Diffusive Flux of Potential Temperature'],
                        'DFyE_SLT': [-1, 1, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Meridional Diffusive Flux of Salinity'],
                        'DFyE_TH': [-1, 1, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Meridional Diffusive Flux of Potential Temperature'],
                        'UVEL': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'], #
                        'VVEL': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity'], #
                        'WVEL': [-5e-4, 5e-4, cm.balance, 'm/s', 'Vertical Velocity'], #
                        'UVELMASS': [-1, 1, cm.balance, 'm/s', 'Zonal Mass-Weighted Component of Velocity'], #
                        'VVELMASS': [-1, 1, cm.balance, 'm/s', 'Meridional Mass-Weighted Component of Velocity'], #
                        'WVELMASS': [-5e-4, 5e-4, cm.balance, 'm/s', 'Vertical Mass-Weighted Component of Velocity'], #
                        'PHIHYD': [-10, 10, cm.curl, 'm$^2$/s$^2$', 'Hydrostatic Pressure Potential Anomaly'],
                        'PHIHYDcR': [-10, 10, cm.curl, 'm$^2$/s$^2$', 'Hydrostatic Pressure Potential Anomaly at Constant r'],
                        'RHOAnoma': [-10, 10, cm.curl, 'kg/m$^3$', 'Density Anomaly'],
                        'KPPdiffS': [0, 0.5, cm.dense, 'm$^2$/s', 'KPP Vertical Salt Diffusion Coefficient'],
                        'KPPdiffT': [0, 0.5, cm.dense, 'm$^2$/s', 'KPP Vertical Heat Diffusion Coefficient'],
                        'KPPviscA': [0, 0.25, cm.dense, 'm$^2$/s', 'KPP Vertical Eddy Viscosity Coefficient']}

    if dim==2:
        return(metadata_dict_2D[field_name])
    if dim==3:
        return(metadata_dict_3D[field_name])

def create_plot(output_dir, output_file, XC_faces, YC_faces, field_set, field_name, date, field_faces, dim):

    fontsize = 16

    fig = plt.figure(figsize=(10, 10))
    plt.style.use('dark_background')

    metadata = field_name_to_plot_metadata(field_name, dim)
    vmin = metadata[0]
    vmax = metadata[1]
    cmap = metadata[2]
    units = metadata[3]
    long_name = metadata[4]

    gs2 = GridSpec(15, 12, left=0.0, right=0.90, hspace=0.05)

    central_lon = -150

    ax = fig.add_subplot(gs2[:-2, :],projection=ccrs.NorthPolarStereo(central_longitude=central_lon))

    # ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)

    for face in [1,2,3,4,5]:
        XC = XC_faces[face]
        YC = YC_faces[face]
        field = field_faces[face]

        if np.any(field!=0):
            print('    - In face '+str(face)+' min: '+str(np.min(field[field!=0]))+', max: '+str(np.max(field[field!=0])))
        else:
            print('    - Face '+str(face)+' == 0')

        if face==4:
            XC_1 = np.copy(XC)
            XC_1[XC_1<0]+=360
            C = plt.pcolormesh(XC_1, YC, field, zorder=1, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
            cbar = plt.colorbar(C, fraction=0.026, pad=0.05)
            cbar.set_label(field_name+' ('+units+')', fontsize=fontsize)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fontsize)
        elif face==3:
            XC_1 = np.copy(XC[:540,:])
            XC_1[XC_1<0]+=360
            YC_1 = YC[:540, :]
            field_1 = field[:540,:]
            plt.pcolormesh(XC_1, YC_1, field_1, zorder=1, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax)

            XC_2 = np.copy(XC[540:, :])
            XC_2[XC_1 < 0] += 360
            YC_2 = YC[540:, :]
            field_2 = field[540:, :]
            plt.pcolormesh(XC_2, YC_2, field_2, zorder=1, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap,
                           vmin=vmin, vmax=vmax)

        #     XC_2 = np.copy(XC[:, :540])
        #     XC_2[XC_2 > 0] -= 360
        #     YC_2 = YC[:, :540]
        #     field_2 = field[:, :540]
        #     plt.pcolormesh(XC_2, YC_2, field_2, zorder=1, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap,
        #                    vmin=vmin, vmax=vmax)
        #     extent = [np.min(XC), np.max(XC), np.min(YC), np.max(YC)]
        #     plt.imshow(field, zorder=1, extent = extent,
        #                transform=ccrs.PlateCarree(), origin='lower', cmap=cmap)

        else:
            plt.pcolormesh(XC, YC, field, zorder=1, transform=ccrs.PlateCarree(), shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax)


    if dim==2:
        plt.title(long_name, fontsize=fontsize)
    if dim==3:
        plt.title(long_name + '\n(Depth Level = 10)', fontsize=fontsize)

    ax.add_feature(cf.LAND, zorder=2, edgecolor='black', facecolor='silver')

    ax.gridlines()

    extent = 4000000

    ax.set_extent((-extent,
                   extent,
                   -extent,
                   extent),
                  crs=ccrs.NorthPolarStereo(central_longitude=central_lon))

    # C = ax1.imshow(field_faces[3], origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)
    # plt.colorbar(C, fraction=0.031, pad=0.04)

    ax2 = fig.add_subplot(gs2[-1, 3:-2])
    year = int(date[:4])
    month = int(date[4:6])
    day = int(date[6:8])
    width = ((datetime(year,month,day)-datetime(year,1,1)).total_seconds())/((datetime(year+1,1,1)-datetime(year,1,1)).total_seconds())
    rect = Rectangle((year, 0), width, 1, fc='silver', ec='white')
    ax2.add_patch(rect)
    ax2.set_xlim([year, year+1])
    ax2.set_ylim([0, 1])
    ax2.set_xticks(np.arange(year, year+1, 1 / 12))
    ax2.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], fontsize=fontsize-2)
    plt.xlabel(str(year), fontsize=fontsize)
    ax2.set_yticks([])

    output_path = os.path.join(output_dir, output_file)
    plt.savefig(output_path,bbox_inches='tight')
    plt.close(fig)

def compile_panels_to_movie(config_dir,config_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L1',config_name,'plots','output',field_name)

    file_name = 'L1_'+field_name+'.mp4'

    os.chdir(panels_dir)
    os.system("ffmpeg -r 5 -i L1_CE_Greenland_"+field_name+"_%04d.png -vcodec mpeg4 -b 1M -y " + file_name)
    os.rename(file_name, os.path.join('..', file_name))

    os.chdir(pwd)

def plot_fields(config_dir, dim, field_number, remove_old, overwrite):

    if dim not in [2,3]:
        raise ValueError('Dimension should be 2 or 3')

    field_set, field_name, field_set_index, is_vector = field_number_to_set_and_name(field_number, dim)

    # time.sleep(field_number)
    if 'plots' not in os.listdir(os.path.join(config_dir, 'N1_1080')):
        os.mkdir(os.path.join(config_dir, 'N1_1080','plots'))
    if field_set not in os.listdir(os.path.join(config_dir, 'N1_1080','plots')):
        os.mkdir(os.path.join(config_dir, 'N1_1080','plots',field_set))
    if field_name not in os.listdir(os.path.join(config_dir, 'N1_1080','plots',field_set)):
        os.mkdir(os.path.join(config_dir, 'N1_1080','plots',field_set, field_name))
    output_dir = os.path.join(config_dir, 'N1_1080','plots',field_set, field_name)

    file_names, dates = get_filename_date_pairs(config_dir, field_set)

    XC_faces, YC_faces, AngleCS_faces, AngleSN_faces = read_faces_geometry(config_dir, is_vector)

    for f in range(len(file_names)):
        file_name = file_names[f]
        date = dates[f]
        print('        - Working on plot for '+field_name+' on '+date+' ('+file_name+')')

        output_file = 'N1_1080_' + field_name + '_' + date + '.png'

        if overwrite:
            continue_to_plot = True
        else:
            if output_file in os.listdir(output_dir):
                continue_to_plot = False
            else:
                continue_to_plot=True

        if continue_to_plot:

            print('    - Reading in the data')
            if is_vector:
                field_faces = read_vector_field_to_faces(config_dir, dim, field_set, file_name, field_set_index, AngleCS_faces, AngleSN_faces)
            else:
                field_faces = read_field_to_faces(config_dir, dim, field_set, file_name, field_set_index)

            print('    - Creating plot for ' + output_file)
            create_plot(output_dir, output_file, XC_faces, YC_faces, field_set, field_name, date, field_faces, dim)

            del field_faces




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config_dir", action="store",
                        help="The directory where the N1, N2, and N3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-d", "--dim", action="store",
                        help="The dimension of the set to plot.", dest="dim",
                        type=int, required=True)

    parser.add_argument("-f", "--field_number", action="store",
                        help="The number of the field to plot (1-24, used for parallel plotting on pleiades).", dest="field_number",
                        type=int, required=True)

    parser.add_argument("-r", "--remove_old", action="store",
                        help="Choose whether to remove all of the old files (1 is true, 0 is false).", dest="remove_old",
                        type=int, required=False, default = 0)

    parser.add_argument("-o", "--overwrite", action="store",
                        help="Choose whether to overwrite old files (1 is true, 0 is false).", dest="overwrite",
                        type=int, required=False, default=0)

    args = parser.parse_args()
    config_dir = args.config_dir
    dim = args.dim
    field_number = args.field_number
    remove_old = args.remove_old
    overwrite = args.overwrite

    if remove_old==0:
        remove_old = False
    else:
        remove_old = True

    if overwrite==0:
        overwrite = False
    else:
        overwrite = True

    plot_fields(config_dir, dim, field_number, remove_old, overwrite)
   

