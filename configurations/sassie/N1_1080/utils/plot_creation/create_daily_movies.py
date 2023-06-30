
import os


def compile_panels_to_movie(plot_folder,var_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(plot_folder,'daily_movies','panels')

    file_name = 'N1_1080_'+var_name+'.mp4'

    os.chdir(panels_dir)
    os.system("ffmpeg -r 5 -i panel_%04d.png -vcodec mpeg4 -b 1M -y " + file_name)
    os.rename(file_name, os.path.join('..', file_name))

    os.chdir(pwd)

plot_folder = '/Volumes/ifenty/Wood/Projects/Ocean_Modeling/Projects/Downscale_ECCOv5_Arctic/MITgcm/' \
              'configurations/arctic_ECCOv5_forcing/N1_1080/plots'

var_sets = {'EXF_day_mean':['EXFaqh','EXFatemp','EXFempmr','EXFevap','EXFhl','EXFhs','EXFlwdn','EXFlwnet',
                            'EXFpreci','EXFqnet','EXFroff','EXFswdn','EXFswnet','EXFtaux','EXFtauy','EXFuwind','EXFvwind'],
            'KPP_hbl_day_mean':['KPPhbl'],
            'KPP_mix_day_mean':['KPPdiffS','KPPdiffT','KPPviscA'],
            'oce_flux_day_mean':['oceFWflx','oceQnet','oceQsw','oceTAUX','oceTAUY','SFLUX','TFLUX'],
            'ocean_state_2D_day_mean':['ETAN','PHIBOT','sIceLoad'],
            'ocean_state_3D_day_mean':['SALT','THETA'],
            'ocean_vel_day_mean':['UVEL','VVEL','WVEL'],
            'phi_3D_day_mean':['PHIHYD','PHIHYDcR','RHOAnoma'],
            'seaice_flux_day_mean':['SIatmFW','SIatmQnt'],
            'seaice_state_day_mean':['SIarea','SIheff','SIhsnow'],
            'seaice_vel_day_mean':['SIuice','SIvice'],
            'tr_adv_r_day_mean':['ADVr_SLT','ADVr_TH'],
            'tr_adv_x_2D_day_mean':['ADVxHEFF','ADVxSNOW','ADVyHEFF','ADVySNOW'],
            'tr_adv_x_3D_day_mean':['ADVx_SLT','ADVx_TH','ADVy_SLT','ADVy_TH'],
            'tr_diff_r_day_mean':['DFrE_SLT','DFrI_SLT','DFrE_TH','DFrI_TH'],
#            'tr_diff_x_day_mean':['DFxE_SLT'],
            'vol_adv_day_mean':['UVELMASS','VVELMASS','WVELMASS']}

# var_sets = {'seaice_state_day_mean':['SIarea','SIheff']}

if 'panels' in os.listdir(os.path.join(plot_folder,'daily_movies')):
    for file_name in os.listdir(os.path.join(plot_folder,'daily_movies','panels')):
        os.remove(os.path.join(plot_folder,'daily_movies','panels',file_name))
    os.rmdir(os.path.join(plot_folder,'daily_movies','panels'))

for var_set in list(var_sets.keys()):
    print('    - Generating movies for the '+var_set+' set')

    for var_name in var_sets[var_set]:

        print('    - Generating movies for '+var_name)

        # get an ordered list of all the files
        var_dir = os.path.join(plot_folder,var_set,var_name)
        file_names = []
        for file_name in os.listdir(var_dir):
            if file_name[0]!='.' and file_name[-4:]=='.png':
                file_names.append(file_name)
        file_names = sorted(file_names)
        print('             - Files '+file_names[0]+' to '+file_names[-1])

        # make a tmp dir
        if 'panels' not in os.listdir(os.path.join(plot_folder,'daily_movies')):
            os.mkdir(os.path.join(plot_folder,'daily_movies','panels'))

        # sym link all of the files in the temp dir with ordered numbers
        for f in range(len(file_names)):
            os.symlink(os.path.join(plot_folder,var_set,var_name,file_names[f]),
                       os.path.join(plot_folder,'daily_movies','panels','panel_'+'{:04d}'.format(f)+'.png'))

        # make the movie
        compile_panels_to_movie(plot_folder,var_name)

        # remove the tmp dir
        for f in range(len(file_names)):
            os.remove(os.path.join(plot_folder,'daily_movies','panels','panel_'+'{:04d}'.format(f)+'.png'))
        os.rmdir(os.path.join(plot_folder,'daily_movies','panels'))

