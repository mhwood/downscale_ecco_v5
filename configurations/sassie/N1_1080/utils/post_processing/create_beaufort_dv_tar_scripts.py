
import os
import numpy as np

# year = '2016'
# numbers = np.arange(639,654)

year = '2019'
numbers = np.arange(718,733)
subset = '2D'

year = '2020'
numbers = np.arange(743,755)
subset = '3D'

var_names_2D_1 = ['ETAN', 'AREA','EMPMR','HEFF','HSNOW']

var_names_2D_2 = ['QNET','SICELOAD', 'UICE','USTRESS','UWIND']

var_names_2D_3 = ['VICE','VSTRESS','VWIND','AQH', 'ATEMP']

var_names_2D_4 = ['EVAP', 'HL', 'HS', 'LWDOWN','LWFLUX','SWFLUX']

var_names_2D_5 = ['SWDOWN', 'PRECIP','RUNOFF', 'KPPFRAC', 'KPPHBL']

var_names_2D = var_names_2D_1 + var_names_2D_2 + var_names_2D_3 + var_names_2D_4 + var_names_2D_5

var_names_3D_6 = ['SALT','THETA']

var_names_3D_7 = ['UVEL','VVEL']

var_names_3D_8 = ['WVEL','KPPDIFFS']

var_names_3D_9 = ['KPPVISCA']

var_names_3D = var_names_3D_6 + var_names_3D_7 + var_names_3D_8 + var_names_3D_9

if subset=='2D':
    var_names = var_names_2D
if subset=='3D':
    var_names = var_names_3D


print_tar_lines = True
print_Lou_copy_lines = True

if print_tar_lines:
    tar_dir = '/Volumes/zachariae/Wood/Projects/Ocean_Modeling/Projects/Downscale_ECCOv5_Arctic/MITgcm/' \
              'configurations/arctic_ECCOv5_forcing/N1_1080/results/beaufort/tar_scripts'

    for var_name in var_names:
        command = 'mkdir ' + var_name
        # print(command)
        command = 'cd ../../../run/dv/'
        print(command)
        for number in numbers:
            # command = 'tar -czvf ../../results/beaufort/diagnostics_vec/'+var_name+'/beaufort_mask_'+var_name+'.000'+str(number)+'0000.tar.gz beaufort_mask_'+var_name+'.000'+str(number)+'*'
            # print(command)
            command = 'tar -czvf ../../results/beaufort/diagnostics_vec/' + var_name + '/beaufort_mask_' + var_name + '.000' + str(
                number) + '0000.tar.gz beaufort/'+var_name+'/beaufort_mask_' + var_name + '.000' + str(number) + '*'
            print(command)
        command = 'cd ../../results/beaufort/diagnostics_vec/'
        print(command)

if print_Lou_copy_lines:
    lou_dir = '/Volumes/zachariae/Wood/Projects/Ocean_Modeling/Projects/Downscale_ECCOv5_Arctic/MITgcm/' \
              'configurations/arctic_ECCOv5_forcing/N1_1080/results/beaufort/Lou_scripts'

    path_name = '/nobackup/mwood7/Arctic/Nested_Models/MITgcm/configurations/N1_1080/results/beaufort/diagnostics_vec'

    output = ''

    for var_name in var_names:
        command = 'mkdir '+var_name
        # print(command)
        output+=command+'\n'
        for number in numbers:
            command = 'echo "mv beaufort_mask_' + var_name + '.000' + str(number) + '0000.tar.gz ' + var_name + '"'
            # print(command)
            output += command + '\n'
            command = 'mv '+path_name+'/'+var_name+'/beaufort_mask_'+var_name+'.000'+str(number)+'0000.tar.gz '+var_name+'/'
            # print(command)
            output += command + '\n'

    output_file = 'mv_zip_files_from_nobackupp_'+year+'_'+subset+'.sh'
    print(output_file)
    f = open(os.path.join(lou_dir, output_file), 'w')
    f.write(output)
    f.close()