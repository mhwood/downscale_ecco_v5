
import numpy as np


var_names = ['AREA','EMPMR','ETAN','HEFF','HSNOW',
             'QNET','SICELOAD', 'UICE','USTRESS','UWIND',
             'VICE','VSTRESS','VWIND','AQH', 'ATEMP',
             'EVAP', 'HL', 'HS', 'LWDOWN','LWFLUX','SWFLUX',
             'SWDOWN', 'PRECIP','RUNOFF', 'KPPFRAC', 'KPPHBL',
             'SALT','THETA',
             'UVEL','VVEL',
             'WVEL','KPPDIFFS',
             'KPPDIFFT', 'KPPVISCA']

for var_name in var_names:
    command = 'mkdir beaufort/' + var_name
    print(command)

    command = 'mv beaufort*' + var_name+'* beaufort/'+var_name
    print(command)



var_names = ['AREA','HEFF','HSNOW','UICE','VICE',
             'SALT','THETA', 'UVEL','VVEL','WVEL','ETAN']

for var_name in var_names:
    command = 'mkdir N2_boundary/' + var_name
    print(command)

    command = 'mv N2_boundary*' + var_name+'* N2_boundary/'+var_name
    print(command)

