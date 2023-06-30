
import os
import numpy as np

def get_plot_var_list():
    plot_vars = ['EXF_day_mean', 'ocean_state_2D_day_mean', 'seaice_flux_day_mean',
                 'tr_adv_x_2D_day_mean', 'vol_adv_day_mean', 'KPP_hbl_day_mean',
                 'ocean_state_3D_day_mean', 'seaice_state_day_mean', 'tr_adv_x_3D_day_mean',
                 'KPP_mix_day_mean', 'ocean_vel_day_mean', 'seaice_vel_day_mean', 'tr_diff_r_day_mean',
                 'oce_flux_day_mean', 'phi_3D_day_mean', 'tr_adv_r_day_mean', 'tr_diff_x_day_mean']

    plot_var_subset_dict = {'EXF_day_mean': ['EXFaqh', 'EXFempmr', 'EXFhl', 'EXFlwdn', 'EXFpreci', 'EXFroff',
                                             'EXFswnet', 'EXFtauy', 'EXFvwind', 'EXFatemp', 'EXFevap', 'EXFhs',
                                             'EXFlwnet', 'EXFqnet', 'EXFswdn','EXFtaux', 'EXFuwind'],
                            'ocean_state_2D_day_mean': ['ETAN', 'PHIBOT', 'sIceLoad'],
                            'seaice_flux_day_mean': ['SIatmFW', 'SIatmQnt'],
                            'tr_adv_x_2D_day_mean': ['ADVxHEFF', 'ADVxSNOW', 'ADVyHEFF', 'ADVySNOW'],
                            'vol_adv_day_mean': ['UVELMASS', 'VVELMASS', 'WVELMASS'],
                            'KPP_hbl_day_mean': ['KPPhbl'],
                            'ocean_state_3D_day_mean': ['SALT', 'THETA'],
                            'seaice_state_day_mean': ['SIarea', 'SIheff', 'SIhsnow'],
                            'tr_adv_x_3D_day_mean': ['ADVx_SLT', 'ADVx_TH', 'ADVy_SLT', 'ADVy_TH'],
                            'KPP_mix_day_mean': ['KPPdiffS', 'KPPdiffT', 'KPPviscA'],
                            'ocean_vel_day_mean': ['UVEL', 'VVEL', 'WVEL'],
                            'seaice_vel_day_mean': ['SIuice', 'SIvice'],
                            'tr_diff_r_day_mean': ['DFrE_SLT', 'DFrE_TH', 'DFrI_SLT', 'DFrI_TH'],
                            'oce_flux_day_mean': ['SFLUX', 'TFLUX', 'oceFWflx', 'oceQnet', 'oceQsw', 'oceTAUX',
                                                  'oceTAUY'],
                            'phi_3D_day_mean': ['PHIHYD', 'PHIHYDcR', 'RHOAnoma'],
                            'tr_adv_r_day_mean': ['ADVr_SLT', 'ADVr_TH'],
                            'tr_diff_x_day_mean': ['DFxE_SLT', 'DFxE_TH', 'DFyE_SLT', 'DFyE_TH']}

    # plot_vars = ['EXF_day_mean']
    # plot_var_subset_dict = {'EXF_day_mean': ['EXFaqh', 'EXFpreci', 'EXFhs', 'EXFswdn', 'EXFtaux', 'EXFuwind']}

    return(plot_vars, plot_var_subset_dict)

def create_plot_tar_scripts(config_dir,years,plot_vars,plot_var_subset_dict):

    for year in years:
        output = 'cd ../'
        for var_set in plot_vars:
            output += '\ncd '+var_set
            for var_name in sorted(plot_var_subset_dict[var_set]):
                output += '\ncd '+var_name
                output += '\ntar -czvf ../../tar_files/'+var_name+'_'+str(year)+'.tar.gz *'+str(year)+'*'
                output += '\ncd ../'
            output += '\ncd ../'

        output_file = os.path.join(config_dir,'plots','tar_scripts','tar_'+str(year)+'.sh')
        f = open(output_file,'w')
        f.write(output)
        f.close()

def create_plot_mv_scripts(config_dir,years,plot_vars,plot_var_subset_dict):
    lou_dir = config_dir + '/plots/Lou_scripts'

    path_name = nobackup_dir + '/plots'

    for year in years:
        output = 'cd ..'
        for var_set in plot_vars:
            output += '\ncd ' + var_set
            for var_name in sorted(plot_var_subset_dict[var_set]):
                command = '\nmv ' + path_name + '/tar_files/' + var_name + '_' + str(year) + '.tar.gz ' + var_name
                output += command
            output += '\ncd ..'
        output += '\ncd mv_scripts'

        output_file = 'mv_files_from_nobackupp_' + str(year) + '.sh'
        f = open(os.path.join(lou_dir, output_file), 'w')
        f.write(output)
        f.close()

def create_N2_boundary_dv_tar_scripts(config_dir, N2_var_names, zip_numbers):
    tar_dir = config_dir + '/results/N2_boundary/tar_scripts'

    for number in zip_numbers:
        output = 'cd ../../../run/dv/N2_boundary'
        for var_name in N2_var_names:
            command = 'cd '+var_name
            output += '\n' + command
            command = 'tar -czvf ../../../../results/N2_boundary/'+var_name+'/N2_boundary_mask_'+\
                      var_name+'.000'+str(number)+'0000.tar.gz N2_boundary_mask_'+var_name+'.000'+str(number)+'*'
            output += '\n'+command
            command = 'cd ../'
            output += '\n'+command
        command = 'cd ../../../results/N2_boundary/tar_scripts'
        output += '\n' + command

        output_file = 'tar_N2_boundary_' + str(number)+'.sh'
        f = open(os.path.join(tar_dir, output_file), 'w')
        f.write(output)
        f.close()

def create_N2_boundary_dv_mv_scripts(config_dir, N2_var_names, zip_numbers):
    lou_dir = config_dir + '/results/N2_boundary/Lou_scripts'

    path_name = nobackup_dir + '/results/N2_boundary'

    for zip_number in zip_numbers:
        output = ''
        for var_name in N2_var_names:
            command = 'mv ' + path_name + '/' + var_name + '/N2_boundary_mask_' + var_name + '.000' + str(
                zip_number) + '0000.tar.gz ../' + var_name + '/'
            output += command + '\n'

        output_file = 'mv_zip_files_from_nobackupp_' + str(zip_number) + '.sh'
        f = open(os.path.join(lou_dir, output_file), 'w')
        f.write(output)
        f.close()

def create_beaufort_dv_tar_scripts():

    for var_name in var_names:
        command = 'mkdir ' + var_name
        # print(command)
        command = 'cd ../../../run/dv/'
        print(command)
        for number in numbers:
            command = 'tar -czvf ../../results/beaufort/diagnostics_vec/'+var_name+'/beaufort_mask_'+var_name+'.000'+str(number)+'0000.tar.gz beaufort_mask_'+var_name+'.000'+str(number)+'*'
            print(command)
        command = 'cd ../../results/beaufort/diagnostics_vec/'
        print(command)

def create_beaufort_dv_mv_scripts(config_dir, nobackup_dir):
    lou_dir = config_dir+'/results/beaufort/Lou_scripts'

    path_name = nobackup_dir+'/results/beaufort/diagnostics_vec'

    output = ''

    for var_name in var_names:
        command = 'mkdir ' + var_name
        # print(command)
        output += command + '\n'
        for number in numbers:
            command = 'mv ' + path_name + '/' + var_name + '/beaufort_mask_' + var_name + '.000' + str(
                number) + '0000.tar.gz ' + var_name + '/'
            # print(command)
            output += command + '\n'

    output_file = 'mv_zip_files_from_nobackupp_' + year + '_' + subset + '.sh'
    f = open(os.path.join(lou_dir, output_file), 'w')
    f.write(output)
    f.close()

config_dir = '/Volumes/ifenty/Wood/Projects/Ocean_Modeling/Projects/Downscale_ECCOv5_Arctic/MITgcm/' \
             'configurations/arctic_ECCOv5_forcing/N1_1080'

nobackup_dir = '/nobackup/mwood7/Arctic/Nested_Models/MITgcm/configurations/N1_1080'

subsets = ['plots']

years = np.arange(2014,2022)
zip_numbers = np.arange(579,789)

if 'plots' in subsets:
    plot_vars, plot_var_subset_dict = get_plot_var_list()
    create_plot_tar_scripts(config_dir,years,plot_vars,plot_var_subset_dict)
    create_plot_mv_scripts(config_dir, years, plot_vars, plot_var_subset_dict)

if 'N2_boundary' in subsets:
    N2_vars = ['AREA','ETAN','HEFF','HSNOW','SALT','THETA','UICE','UVEL','VICE','VVEL','WVEL']
    create_N2_boundary_dv_tar_scripts(config_dir, N2_vars, zip_numbers)
    create_N2_boundary_dv_mv_scripts(config_dir, N2_vars, zip_numbers)








