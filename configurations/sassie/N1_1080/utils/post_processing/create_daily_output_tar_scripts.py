
import os
import numpy as np


numbers = np.arange(1508,1542)

var_names = ['EXF_day_mean','ocean_state_3D_day_mean','seaice_flux_day_mean','tr_adv_r_day_mean',#'tr_diff_x_day_mean',
             'KPP_hbl_day_mean','oce_flux_day_mean','seaice_state_day_mean','tr_adv_x_2D_day_mean','vol_adv_day_mean',
             'KPP_mix_day_mean','ocean_state_2D_day_mean','ocean_vel_day_mean','tr_adv_x_3D_day_mean',
             'phi_3D_day_mean','seaice_vel_day_mean','tr_diff_r_day_mean']
var_names = sorted(var_names)

print_tar_lines = True

print_Lou_copy_lines = True

if print_tar_lines:
    tar_dir = '/Volumes/zachariae/Research/Projects/Ocean_Modeling/Projects/Downscale_ECCOv5_Arctic/MITgcm/configurations/arctic_ECCOv5_forcing/N1_1080/results/N1_domain/tar_scripts'

    output_file = 'tar_diag_files_' + str(numbers[0])+'_'+str(numbers[-1]) + '.sh'
    output = ''

    for number in numbers:

        for var_name in var_names:
            command = 'echo "tarring '+var_name+' files for set '+str(number)+'"'
            output += command + '\n'

            if number >= 1000:
                command = 'cd ../../run_60/diags/'+var_name
            else:
                command = 'cd ../../run/diags/' + var_name
            # print(command)
            output += command + '\n'

            if number>=1000:
                command = 'tar -czvf ../../../results/N1_domain/'+var_name+'/'+var_name+'.00'+str(number)+'0000.tar.gz '+var_name+'.00'+str(number)+'*'
            else:
                command = 'tar -czvf ../../../results/N1_domain/' + var_name + '/' + var_name + '.00' + str(number) + '0000.tar.gz ' + var_name + '.00' + str(number) + '*'
            # print(command)
            output += command + '\n'

            command = 'cd ../../../results/N1_domain/'
            # print(command)
            output += command + '\n'

    f = open(os.path.join(tar_dir, output_file), 'w')
    f.write(output)
    f.close()

if print_Lou_copy_lines:

    lou_dir = '/Volumes/zachariae/Research/Projects/Ocean_Modeling/Projects/Downscale_ECCOv5_Arctic/MITgcm/' \
              'configurations/arctic_ECCOv5_forcing/N1_1080/results/N1_domain/Lou_scripts'

    path_name = '/nobackup/mwood7/Arctic/Nested_Models/MITgcm/configurations/N1_1080/results/N1_domain'

    for number in numbers:
        output_file = 'mv_zip_files_from_nobackupp_'+str(number)+'.sh'
        output = ''
        for var_name in var_names:
            if number >= 1000:
                command = 'mv '+path_name+'/'+var_name+'/'+var_name+'.00'+str(number)+'0000.tar.gz '+var_name+'/'
            else:
                command = 'mv ' + path_name + '/' + var_name + '/' + var_name + '.00' + str(number) + '0000.tar.gz ' + var_name + '/'
            output+=command+'\n'

        f = open(os.path.join(lou_dir,output_file),'w')
        f.write(output)
        f.close()