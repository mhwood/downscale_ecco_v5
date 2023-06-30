import os

dv_dir = os.path.join('..','run','dv')

boundary_mask_names = ['mask_41N_i1','mask_41N_i2','mask_41N_i3','mask_41N']
boundary_variables = ['THETA','SALT','UVEL','VVEL','ETAN']

surface_mask_names = []#['surface']
surface_variables = ['APRESS', 'AQH', 'ATEMP', 'LWDOWN', 'PRECIP', 'SWDOWN',
                    'USTRESS', 'VSTRESS', 'WSPEED', 'RUNOFF',
                     'ETAN','SALT','THETA','UVEL','VVEL']

for mask_name in boundary_mask_names:
    if mask_name not in os.listdir(dv_dir):
        os.mkdir(os.path.join(dv_dir,mask_name))
    for var_name in boundary_variables:
        if var_name not in os.listdir(os.path.join(dv_dir,mask_name)):
            os.mkdir(os.path.join(dv_dir,mask_name,var_name))

    for file_name in os.listdir(dv_dir):
        if file_name==mask_name+'.bin':
            os.rename(os.path.join(dv_dir,file_name),os.path.join(dv_dir,mask_name,file_name))
        if file_name[:len(mask_name)]==mask_name and file_name.split('.')[0].split('_')[-1] in boundary_variables:
            os.rename(os.path.join(dv_dir, file_name), os.path.join(dv_dir, mask_name,file_name.split('.')[0].split('_')[-1], file_name))


for mask_name in surface_mask_names:
    if mask_name not in os.listdir(dv_dir):
        os.mkdir(os.path.join(dv_dir,mask_name))
    for var_name in surface_variables:
        if var_name not in os.listdir(os.path.join(dv_dir,mask_name)):
            os.mkdir(os.path.join(dv_dir,mask_name,var_name))

    for file_name in os.listdir(dv_dir):
        if file_name==mask_name+'.bin':
            os.rename(os.path.join(dv_dir,file_name),os.path.join(dv_dir,mask_name,file_name))
        if file_name[:len(mask_name)]==mask_name and file_name.split('.')[0].split('_')[-1] in surface_variables:
            os.rename(os.path.join(dv_dir, file_name), os.path.join(dv_dir, mask_name,file_name.split('.')[0].split('_')[-1], file_name))
