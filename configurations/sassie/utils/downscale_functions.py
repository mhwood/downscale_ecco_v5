
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import netCDF4 as nc4
import xesmf as xe
import xarray as xr


def create_hFacC_grid(bathy,delR, hFacMin=0.2, hFacMinDr=5.0):
    # This is from MITgcm inside ini_masks_etc.F
    # Note: R_low is just the bathy
    # Note: drF is the grid spacing (provided)
    # Note: recip_drF is the reciprocal of drF
    # Note: Ro_surf is the z coord of the surface (essentially always 0 for the ocean?)
    # C--   Calculate lopping factor hFacC : over-estimate the part inside of the domain
    # C     taking into account the lower_R Boundary (Bathymetry / Top of Atmos)
    #         DO k=1, Nr
    #          hFacMnSz = MAX( hFacMin, MIN(hFacMinDr*recip_drF(k),oneRL) )
    #          DO j=1-OLy,sNy+OLy
    #           DO i=1-OLx,sNx+OLx
    # C      o Non-dimensional distance between grid bound. and domain lower_R bound.
    #            hFac_loc = (rF(k)-R_low(i,j,bi,bj))*recip_drF(k)
    # C      o Select between, closed, open or partial (0,1,0-1)
    #            hFac_loc = MIN( MAX( hFac_loc, zeroRL ) , oneRL )
    # C      o Impose minimum fraction and/or size (dimensional)
    #            IF ( hFac_loc.LT.hFacMnSz*halfRL .OR.
    #      &          R_low(i,j,bi,bj).GE.Ro_surf(i,j,bi,bj) ) THEN
    #              hFacC(i,j,k,bi,bj) = zeroRS
    #            ELSE
    #              hFacC(i,j,k,bi,bj) = MAX( hFac_loc, hFacMnSz )
    #            ENDIF
    #           ENDDO
    #          ENDDO
    #         ENDDO

    # Define grids with same names as those in MITgcm
    RU = -1 * np.cumsum(delR)
    RL = np.concatenate([[0], RU[:-1]])
    R_low = bathy
    drF = RL - RU
    recip_drF = 1 / drF

    # Pythonize the above loops
    hFacC = np.zeros((len(RL), np.shape(bathy)[0], np.shape(bathy)[1]))
    for k in range(len(RL)):
        print('       - Working on level '+str(k)+' of '+str(len(RL)))
        hFacMnSz = np.max([hFacMin, np.min([hFacMinDr * recip_drF[k], 1])])
        for i in range(np.shape(bathy)[0]):
            for j in range(np.shape(bathy)[1]):
                #      o Non-dimensional distance between grid bound. and domain lower_R bound.
                hFac_loc = (RL[k] - R_low[i, j]) * recip_drF[k]
                #      o Select between, closed, open or partial (0,1,0-1)
                hFac_loc = np.min([np.max([hFac_loc, 0]), 1])
                #      o Impose minimum fraction and/or size (dimensional)
                if hFac_loc <= hFacMnSz * 0.5 or R_low[i, j] >= 0:
                    hFacC[k, i, j] = 0
                else:
                    hFacC[k, i, j] = np.max([hFac_loc, hFacMnSz])

    return(hFacC)

def create_hFacS_grid(bathy,delR, hFacMin=0.2, hFacMinDr=5.0):
    # This is from MITgcm inside ini_masks_etc.F
    # Note: R_low is just the bathy
    # Note: drF is the grid spacing (provided)
    # Note: recip_drF is the reciprocal of drF
    # Note: Ro_surf is the z coord of the surface (essentially always 0 for the ocean?)

    RU = -1 * np.cumsum(delR)
    RL = np.concatenate([[0], RU[:-1]])
    R_low = bathy
    drF = RL - RU
    recip_drF = 1 / drF

    # rLowS = np.zeros((np.shape(bathy)[0], np.shape(bathy)[1]))
    rLowS = np.copy(bathy)
    for j in range(1,np.shape(rLowS)[0]):
        for i in range(np.shape(rLowS)[1]):
            rLowS[j,i] = np.max([R_low[j-1,i],R_low[j,i]])

    hFacS = np.zeros((len(delR), np.shape(bathy)[0], np.shape(bathy)[1]))
    for k in range(len(delR)):
        print('       - Working on level ' + str(k) + ' of ' + str(len(RL)))
        hFacMnSz = np.max([hFacMin, np.min([hFacMinDr * recip_drF[k], 1])])
        for j in range(np.shape(rLowS)[0]):
            for i in range(np.shape(rLowS)[1]):
                hFac1tmp = (RL[k] - rLowS[j,i]) * recip_drF[k]
                hFac_loc = np.min([hFac1tmp, 1])
    #      o Impose minimum fraction and/or size (dimensional)
                if hFac_loc<hFacMnSz*0.5 or rLowS[j,i]>=0:
                    hFac1tmp = 0
                else:
                    hFac1tmp = np.max([hFac_loc, hFacMnSz])
    #      o Reduce the previous fraction : substract the outside fraction
    #        (i.e., beyond reference (=at rest) surface position rSurfS)
                hFac2tmp = ( RL[k]-0 )*recip_drF[k]
                hFac_loc = hFac1tmp - np.max([hFac2tmp, 0])
    #      o Impose minimum fraction and/or size (dimensional)
                if hFac_loc<hFacMnSz*0.5:
                    hFacS[k,j,i]=0
                else:
                    hFacS[k,j,i]=np.max([hFac_loc, hFacMnSz])

    return(hFacS)

def create_hFacW_grid(bathy,delR, hFacMin=0.2, hFacMinDr=5.0):
    # This is from MITgcm inside ini_masks_etc.F

    RU = -1 * np.cumsum(delR)
    RL = np.concatenate([[0], RU[:-1]])
    R_low = bathy
    drF = RL - RU
    recip_drF = 1 / drF

    rLowW = np.copy(bathy)
    for j in range(np.shape(rLowW)[0]):
        for i in range(1,np.shape(rLowW)[1]):
            rLowW[j,i] = np.max([R_low[j,i-1],R_low[j,i]])

    hFacW = np.zeros((len(delR), np.shape(bathy)[0], np.shape(bathy)[1]))
    for k in range(len(delR)):
        print('       - Working on level ' + str(k) + ' of ' + str(len(RL)))
        hFacMnSz = np.max([hFacMin, np.min([hFacMinDr * recip_drF[k], 1])])
        for j in range(np.shape(rLowW)[0]):
            for i in range(np.shape(rLowW)[1]):
                hFac1tmp = (RL[k] - rLowW[j,i]) * recip_drF[k]
                hFac_loc = np.min([hFac1tmp, 1])
    #      o Impose minimum fraction and/or size (dimensional)
                if hFac_loc<hFacMnSz*0.5 or rLowW[j,i]>=0:
                    hFac1tmp = 0
                else:
                    hFac1tmp = np.max([hFac_loc, hFacMnSz])
    #      o Reduce the previous fraction : substract the outside fraction
    #        (i.e., beyond reference (=at rest) surface position rSurfS)
                hFac2tmp = ( RL[k]-0 )*recip_drF[k]
                hFac_loc = hFac1tmp - np.max([hFac2tmp, 0])
    #      o Impose minimum fraction and/or size (dimensional)
                if hFac_loc<hFacMnSz*0.5:
                    hFacW[k,j,i]=0
                else:
                    hFacW[k,j,i]=np.max([hFac_loc, hFacMnSz])

    return(hFacW)

def create_2D_wet_grid(bathy, delR0, hFac='C', hFacMin=0.2, hFacMinDr=5.0):

    delR = np.array([[delR0]])

    if hFac=='C':
        hFacGrid = create_hFacC_grid(bathy, delR, hFacMin, hFacMinDr)
    if hFac=='S':
        hFacGrid = create_hFacS_grid(bathy, delR, hFacMin, hFacMinDr)
    if hFac=='W':
        hFacGrid = create_hFacW_grid(bathy, delR, hFacMin, hFacMinDr)

    wet_grid = np.copy(hFacGrid)
    wet_grid[wet_grid>0]=1
    wet_grid=wet_grid[0,:,:]

    return(wet_grid)

def create_3D_wet_grid(bathy, delR, hFac='C', hFacMin=0.2, hFacMinDr=5.0):

    if hFac=='C':
        hFacGrid = create_hFacC_grid(bathy, delR, hFacMin, hFacMinDr)
    if hFac=='S':
        hFacGrid = create_hFacS_grid(bathy, delR, hFacMin, hFacMinDr)
    if hFac=='W':
        hFacGrid = create_hFacW_grid(bathy, delR, hFacMin, hFacMinDr)

    wet_grid = np.copy(hFacGrid)
    wet_grid[wet_grid>0]=1

    return(wet_grid)

def interpolate_var_grid_to_new_depth_levels(var_grid,wet_grid,delR_in,delR_out):

    Z_bottom_in = np.cumsum(delR_in)
    Z_top_in = np.concatenate([np.array([0]), Z_bottom_in[:-1]])
    Z_in = (Z_bottom_in + Z_top_in) / 2

    Z_bottom_out = np.cumsum(delR_out)
    Z_top_out = np.concatenate([np.array([0]), Z_bottom_out[:-1]])
    Z_out = (Z_bottom_out + Z_top_out) / 2


    if len(np.shape(var_grid))==3:
        new_var_grid = np.zeros((np.size(delR_out), np.shape(var_grid)[1],np.shape(var_grid)[2]))
        new_wet_grid = np.zeros((np.size(delR_out), np.shape(var_grid)[1], np.shape(var_grid)[2]))
        for i in range(np.shape(var_grid)[1]):
            for j in range(np.shape(var_grid)[2]):
                test_profile = var_grid[:, i, j]
                if np.sum(test_profile != 0) > 1:
                    set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                              bounds_error=False, fill_value=np.nan)
                    new_profile = set_int_linear(Z_out)

                    new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                    if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                        first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                        bottom_value = new_profile[~np.isnan(new_profile)][-1]
                        new_profile[np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                    new_var_grid[:, i, j] = new_profile

                if np.sum(test_profile == 0) == 1:
                    new_var_grid[0, i, j] = var_grid[0, i, j]

    elif len(np.shape(var_grid)) == 4:
        new_var_grid = np.zeros((np.shape(var_grid)[0],np.size(delR_out), np.shape(var_grid)[2], np.shape(var_grid)[3]))
        for t in range(np.shape(var_grid)[0]):
            for i in range(np.shape(var_grid)[2]):
                for j in range(np.shape(var_grid)[3]):

                    test_profile = var_grid[t, :, i, j]
                    if np.sum(test_profile != 0) > 1:
                        set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                                  bounds_error=False, fill_value=np.nan)
                        new_profile = set_int_linear(Z_out)

                        new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                        if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                            first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                            bottom_value = new_profile[~np.isnan(new_profile)][-1]
                            new_profile[
                                np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                        new_var_grid[t, :, i, j] = new_profile

                    if np.sum(test_profile == 0) == 1:
                        new_var_grid[t, 0, i, j] = var_grid[t, 0, i, j]
    else:
        raise ValueError('The input array should be dim 3 or 4')

    if np.shape(wet_grid)[0]!=len(delR_out):
        new_wet_grid = np.zeros((np.size(delR_out), np.shape(wet_grid)[1], np.shape(wet_grid)[2]))
        for i in range(np.shape(wet_grid)[1]):
            for j in range(np.shape(wet_grid)[2]):
                test_profile = wet_grid[:, i, j]
                if np.sum(test_profile != 0) > 1:
                    set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                              bounds_error=False, fill_value=np.nan)
                    new_profile = set_int_linear(Z_out)

                    new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                    if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                        first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                        bottom_value = new_profile[~np.isnan(new_profile)][-1]
                        new_profile[np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                    new_wet_grid[:, i, j] = new_profile

                if np.sum(test_profile == 0) == 1:
                    new_wet_grid[0, i, j] = wet_grid[0, i, j]
        new_wet_grid[np.isnan(new_wet_grid)] = 0
        new_wet_grid = np.round(new_wet_grid).astype(int)

        new_var_grid[np.isnan(new_var_grid)] = 0
    else:
        new_wet_grid = wet_grid


    return(new_var_grid,new_wet_grid)


def interpolate_var_grid_faces_to_new_depth_levels(var_grid_faces,wet_grid_faces,delR_in,delR_out):

    Z_bottom_in = np.cumsum(delR_in)
    Z_top_in = np.concatenate([np.array([0]), Z_bottom_in[:-1]])
    Z_in = (Z_bottom_in + Z_top_in) / 2

    Z_bottom_out = np.cumsum(delR_out)
    Z_top_out = np.concatenate([np.array([0]), Z_bottom_out[:-1]])
    Z_out = (Z_bottom_out + Z_top_out) / 2

    var_grid_faces_out = {}

    for face in range(1,6):
        print('          - Interpolating face '+str(face)+' to new depths')
        var_grid = var_grid_faces[face]
        wet_grid = wet_grid_faces[face]
        if len(np.shape(var_grid))==3:
            new_var_grid = np.zeros((np.size(delR_out),np.shape(var_grid)[1],np.shape(var_grid)[2]))
            for i in range(np.shape(var_grid)[1]):
                for j in range(np.shape(var_grid)[2]):
                    test_profile = var_grid[:, i, j]
                    if np.sum(test_profile != 0) > 1:
                        set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                                  bounds_error=False, fill_value=np.nan)
                        new_profile = set_int_linear(Z_out)

                        new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                        if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                            first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                            bottom_value = new_profile[~np.isnan(new_profile)][-1]
                            new_profile[
                                np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                        new_var_grid[:, i, j] = new_profile

                    if np.sum(test_profile == 0) == 1:
                        new_var_grid[0, i, j] = var_grid[0, i, j]

        elif len(np.shape(var_grid)) == 4:
            new_var_grid = np.zeros((np.shape(var_grid)[0],np.size(delR_out), np.shape(var_grid)[2], np.shape(var_grid)[3]))
            print(np.shape(var_grid),np.shape(new_var_grid))
            for t in range(np.shape(var_grid)[0]):
                for i in range(np.shape(var_grid)[2]):
                    for j in range(np.shape(var_grid)[3]):

                        test_profile = var_grid[t, :, i, j]
                        if np.sum(test_profile != 0) > 1:
                            set_int_linear = interp1d(Z_in[test_profile != 0], test_profile[test_profile != 0],
                                                      bounds_error=False, fill_value=np.nan)
                            new_profile = set_int_linear(Z_out)

                            new_profile[np.abs(Z_out) < np.abs(Z_in[0])] = new_profile[~np.isnan(new_profile)][0]
                            if np.size(np.abs(Z_in[test_profile == 0])) > 0:
                                first_zero_depth = np.abs(Z_in[test_profile == 0])[0]
                                bottom_value = new_profile[~np.isnan(new_profile)][-1]
                                new_profile[
                                    np.logical_and(np.isnan(new_profile), np.abs(Z_out) < first_zero_depth)] = bottom_value

                            new_var_grid[t, :, i, j] = new_profile

                        if np.sum(test_profile == 0) == 1:
                            new_var_grid[t, 0, i, j] = var_grid[t, 0, i, j]
        else:
            raise ValueError('The input array should be dim 3 or 4')



                    # plt.plot(Z_in, test_profile,'g-',linewidth=0.5)
                    # plt.plot(Z_out, new_profile,'b.')
                    # # plt.plot(Z_bottom_out, wet_grid[:,i,j])
                    # plt.show()
        var_grid_faces_out[face] = new_var_grid

    return(var_grid_faces_out)

def spread_var_horizontally_in_wet_grid(var_grid,wet_grid):
    rows = np.arange(np.shape(var_grid)[0])
    cols = np.arange(np.shape(var_grid)[1])
    Cols,Rows = np.meshgrid(cols,rows)

    is_remaining = np.logical_and(var_grid==0,wet_grid==1)
    n_remaining = np.sum(is_remaining)
    continue_iter = True
    for i in range(n_remaining):
        if continue_iter:
            Wet_Rows = Rows[wet_grid == 1]
            Wet_Cols = Cols[wet_grid == 1]
            Wet_Vals = var_grid[wet_grid == 1]
            Wet_Rows = Wet_Rows[Wet_Vals != 0]
            Wet_Cols = Wet_Cols[Wet_Vals != 0]
            Wet_Vals = Wet_Vals[Wet_Vals != 0]

            if len(Wet_Vals)>0:

                rows_remaining,cols_remaining = np.where(is_remaining)
                for ri in range(n_remaining):
                    row = rows_remaining[ri]
                    col = cols_remaining[ri]
                    row_col_dist = ((Wet_Rows.astype(float)-row)**2 + (Wet_Cols.astype(float)-col)**2)**0.5
                    closest_index = np.argmin(row_col_dist)
                    if row_col_dist[closest_index]<np.sqrt(2):
                        var_grid[row,col] = Wet_Vals[closest_index]

                is_remaining = np.logical_and(var_grid == 0, wet_grid == 1)
                n_remaining_now = np.sum(is_remaining)
                if n_remaining_now<n_remaining:
                    n_remaining = n_remaining_now
                else:
                    n_remaining = n_remaining_now
                    continue_iter=False

            else:
                continue_iter = False

    return(var_grid,n_remaining)

def spread_var_vertically_in_wet_grid(full_grid,level_grid,wet_grid,level,mean_vertical_difference):

    # if mean_vertical_difference!=0:
    #     print('Using a mean vertical difference of '+str(mean_vertical_difference))

    if level==0:
        bad_row, bad_col = np.where(np.logical_and(level_grid == 0, wet_grid == 1))
        plt.subplot(1, 2, 1)
        plt.imshow(level_grid,origin='lower')
        plt.subplot(1,2,2)
        plt.imshow(wet_grid,origin='lower')
        plt.show()
        raise ValueError('Cannot spread vertically in the surface layer e.g. at row='+str(bad_row[0])+', col='+str(bad_col[0]))

    is_remaining = np.logical_and(level_grid==0,wet_grid==1)
    rows_remaining, cols_remaining = np.where(is_remaining)
    for ri in range(len(rows_remaining)):
        row = rows_remaining[ri]
        col = cols_remaining[ri]
        level_grid[row,col] = full_grid[level-1,row,col]+mean_vertical_difference

    return(level_grid)

def read_regridder_functions_from_nc(input_dir, N0_XC_faces, N0_YC_faces, N1_XC_faces, N1_YC_faces):

    regridding_faces = {}

    for face in range(1,6):
        ds_in = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N0_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N0_YC_faces[face]).T),
            ))

        ds_out = xr.Dataset(
            coords=dict(
                lon=(["x", "y"], np.asfortranarray(N1_XC_faces[face]).T),
                lat=(["x", "y"], np.asfortranarray(N1_YC_faces[face]).T),
            ))

        filename = os.path.join(input_dir, 'regridding_weights_face_'+str(face)+'.nc')
        regridder = xe.Regridder(ds_in, ds_out, method='bilinear', reuse_weights=True, filename = filename)
        regridding_faces[face] = regridder

    return(regridding_faces)


def downscale_2D_field(L0_XC, L0_YC, L0_var, L0_wet_grid, L0_wet_grid_on_L1,
                       L1_XC, L1_YC, L1_wet_grid,spread_horizontally=True,remove_zeros=True):

    # try to interpolate everything using a linear interpolation first
    L0_points = np.hstack([np.reshape(L0_XC, (np.size(L0_XC), 1)),
                              np.reshape(L0_YC, (np.size(L0_YC), 1))])
    L0_values = np.reshape(L0_var, (np.size(L0_var), 1))


    # L0_wet_grid = np.reshape(L0_wet_grid, (np.size(L0_wet_grid), 1))

    if remove_zeros:
        L0_points = L0_points[L0_values[:, 0] != 0, :]
        L0_values = L0_values[L0_values[:, 0] != 0, :]

    grid = griddata(L0_points, L0_values, (L1_XC, L1_YC), method='linear',fill_value=0)
    grid = grid[:, :, 0]

    # mask out any values which should be 0'd based on the old bathy
    # grid[L0_wet_grid_on_L1 == 0] = 0

    # mask out any values which should be 0'd based on the new bathy
    grid[L1_wet_grid == 0] = 0

    # plt.imshow(grid, origin='lower')
    # plt.show()

    #spread the the variable outward to new wet cells
    if spread_horizontally:
        grid, _ = spread_var_horizontally_in_wet_grid(grid, L1_wet_grid)

    # C = plt.imshow(downscaled_grid,origin='lower',
    #                vmin=np.min(downscaled_grid[downscaled_grid!=0]),vmax=np.max(downscaled_grid[downscaled_grid!=0]))
    # plt.colorbar(C)
    # plt.show()
    
    return(grid)


def downscale_3D_field(L0_XC, L0_YC, L0_var, L0_wet_grid, L0_wet_grid_on_L1,
                       XC_subset, YC_subset, L1_wet_grid,
                       mean_vertical_difference=0,fill_downward=True,
                       printing=False,remove_zeros=True, testing = False):

    # full_grid = np.zeros_like(L1_wet_grid).astype(float)
    full_grid = np.zeros((np.shape(L0_wet_grid)[0],np.shape(XC_subset)[0],np.shape(XC_subset)[1]))

    if testing:
        K=1
    else:
        K=np.shape(L0_wet_grid)[0]

    for k in range(K):

        continue_to_interpolation = True

        if continue_to_interpolation:
            if printing:
                print('                - Working on level ' + str(k) + ' of ' + str(np.shape(L1_wet_grid)[0])+' ('+str(np.sum(L1_wet_grid[k, :, :] > 0))+' nonzero points found)')
            if np.any(L1_wet_grid[k, :, :] > 0):
                # take an initial stab at the interpolation
                L0_points = np.hstack([np.reshape(L0_XC, (np.size(L0_XC), 1)),
                                       np.reshape(L0_YC, (np.size(L0_YC), 1))])
                L0_values = np.reshape(L0_var[k, :, :], (np.size(L0_var[k, :, :]), 1))
                L0_wet_grid_vert = np.reshape(L0_wet_grid[k, :, :], (np.size(L0_wet_grid[k, :, :]), 1))
                L0_points = L0_points[L0_wet_grid_vert[:,0] != 0, :]
                L0_values = L0_values[L0_wet_grid_vert[:,0] != 0, :]
                if remove_zeros:
                    L0_points = L0_points[L0_values[:, 0] != 0, :]
                    L0_values = L0_values[L0_values[:, 0] != 0, :]

                if len(L0_points)>4:
                    grid = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='linear',fill_value=0)
                    grid = grid[:, :, 0]
                    # grid_nearest = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='nearest', fill_value=0)
                    # grid_nearest[:,:,0]
                    if not np.any(grid!=0):
                        grid = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='nearest', fill_value=0)
                        grid = grid[:, :, 0]
                else:
                    grid = np.zeros_like(XC_subset).astype(float)

                # if k==0:
                #     plt.imshow(grid,origin='lower')
                #     plt.show()

                # mask out any values which should be 0'd based on the old bathy
                grid[L0_wet_grid_on_L1[k, :, :] == 0] = 0

                # mask out any values which should be 0'd based on the new bathy
                grid[L1_wet_grid[k, :, :] == 0] = 0

                # spread the the variable outward to new wet cells
                grid, n_remaining = spread_var_horizontally_in_wet_grid(grid, L1_wet_grid[k, :, :])

                # if there are still values which need to be filled, spread downward
                if n_remaining > 0 and fill_downward and k>0:
                    grid = spread_var_vertically_in_wet_grid(full_grid, grid, L1_wet_grid[k, :, :], k, mean_vertical_difference)

                # if, for whatever reason, there are still values to be filled, then fill em with the nearest neighbor
                if n_remaining>0:
                    if len(L0_points) > 0:
                        grid_nearest = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='nearest', fill_value=0)
                        grid_nearest = grid_nearest[:, :, 0]
                        indices = np.logical_and(grid==0,L1_wet_grid[k, :, :] != 0)
                        grid[indices] = grid_nearest[indices]

                # # if there are any issues in the surface layer, then fill with a close point
                # if n_remaining > 0 and fill_downward and k == 0:
                #     bad_row, bad_col = np.where(np.logical_and(grid == 0, L1_wet_grid[0,:,:] == 1))
                #     good_row, good_col = np.where(np.logical_and(grid > 0, L1_wet_grid[0,:,:] == 1))
                #     for ri in range(len(bad_row)):
                #         dist = (bad_row[ri]-good_row)**2 + (bad_col[ri]-good_col)**2
                #         fill_index = np.argmin(dist)
                #         fill_val = grid[good_row[fill_index],good_col[fill_index]]
                #         # print('Filling row '+str(bad_row[ri])+', col '+str(bad_col[ri])+
                #         #       ' with value = '+str(fill_val)+' from row '+str(good_row[fill_index])+
                #         #       ', col '+str(good_col[fill_index]))
                #         grid[bad_row[ri],bad_col[ri]] = fill_val
            else:
                grid = np.zeros_like(XC_subset).astype(float)

        full_grid[k, :, :] = grid[:, :]

    return(full_grid)

def downscale_3D_field_with_zeros(L0_XC, L0_YC, L0_var, L0_wet_grid, L0_wet_grid_on_L1,
                                   XC_subset, YC_subset, L1_wet_grid,
                                   mean_vertical_difference=0,fill_downward=True,printing=False,remove_zeros=True):

    # full_grid = np.zeros_like(L1_wet_grid).astype(float)
    full_grid = np.zeros((np.shape(L0_wet_grid)[0],np.shape(XC_subset)[0],np.shape(XC_subset)[1]))

    for k in range(np.shape(L0_wet_grid)[0]):

        continue_to_interpolation = True

        tiny_value = 1e-14

        if continue_to_interpolation:
            if printing:
                print('                - Working on level ' + str(k) + ' of ' + str(np.shape(L1_wet_grid)[0])+' ('+str(np.sum(L1_wet_grid[k, :, :] > 0))+' nonzero points found)')
            if np.any(L1_wet_grid[k, :, :] > 0):
                # take an initial stab at the interpolation
                L0_points = np.hstack([np.reshape(L0_XC, (np.size(L0_XC), 1)),
                                       np.reshape(L0_YC, (np.size(L0_YC), 1))])
                L0_values = np.reshape(L0_var[k, :, :], (np.size(L0_var[k, :, :]), 1))
                L0_wet_grid_vert = np.reshape(L0_wet_grid[k, :, :], (np.size(L0_wet_grid[k, :, :]), 1))
                L0_points = L0_points[L0_wet_grid_vert[:,0] != 0, :]
                L0_values = L0_values[L0_wet_grid_vert[:,0] != 0, :]

                # fill the zeros with a very tiny value
                L0_values[L0_values[:, 0] == 0, :] = tiny_value

                if len(L0_points)>4:
                    grid = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='linear',fill_value=0)
                    grid = grid[:, :, 0]
                    # grid_nearest = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='nearest', fill_value=0)
                    # grid_nearest[:,:,0]
                    if not np.any(grid!=0):
                        grid = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='nearest', fill_value=0)
                        grid = grid[:, :, 0]
                else:
                    grid = np.zeros_like(XC_subset).astype(float)

                # if k==0:
                #     plt.imshow(grid,origin='lower')
                #     plt.show()

                # mask out any values which should be 0'd based on the old bathy
                grid[L0_wet_grid_on_L1[k, :, :] == 0] = 0

                # mask out any values which should be 0'd based on the new bathy
                grid[L1_wet_grid[k, :, :] == 0] = 0

                # spread the the variable outward to new wet cells
                grid, n_remaining = spread_var_horizontally_in_wet_grid(grid, L1_wet_grid[k, :, :])

                # if there are still values which need to be filled, spread downward
                if n_remaining > 0 and fill_downward and k>0:
                    grid = spread_var_vertically_in_wet_grid(full_grid, grid, L1_wet_grid[k, :, :], k, mean_vertical_difference)

                # if, for whatever reason, there are still values to be filled, then fill em with the nearest neighbor
                if n_remaining>0:
                    if len(L0_points) > 0:
                        grid_nearest = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='nearest', fill_value=0)
                        grid_nearest = grid_nearest[:, :, 0]
                        indices = np.logical_and(grid==0,L1_wet_grid[k, :, :] != 0)
                        grid[indices] = grid_nearest[indices]

                # now, add the end of it, make all of the tiny values actually 0
                grid[np.abs(grid)<=2*tiny_value] = 0

                # # if there are any issues in the surface layer, then fill with a close point
                # if n_remaining > 0 and fill_downward and k == 0:
                #     bad_row, bad_col = np.where(np.logical_and(grid == 0, L1_wet_grid[0,:,:] == 1))
                #     good_row, good_col = np.where(np.logical_and(grid > 0, L1_wet_grid[0,:,:] == 1))
                #     for ri in range(len(bad_row)):
                #         dist = (bad_row[ri]-good_row)**2 + (bad_col[ri]-good_col)**2
                #         fill_index = np.argmin(dist)
                #         fill_val = grid[good_row[fill_index],good_col[fill_index]]
                #         # print('Filling row '+str(bad_row[ri])+', col '+str(bad_col[ri])+
                #         #       ' with value = '+str(fill_val)+' from row '+str(good_row[fill_index])+
                #         #       ', col '+str(good_col[fill_index]))
                #         grid[bad_row[ri],bad_col[ri]] = fill_val
            else:
                grid = np.zeros_like(XC_subset).astype(float)

        full_grid[k, :, :] = grid[:, :]

    return(full_grid)


def downscale_2D_field_faces(L0_var_faces,
                             regridder_faces, L0_wet_grid_on_L1_faces,
                             L1_wet_grid_faces,
                             spread_horizontally=True):

    grid_faces = {}

    for face in range(1,6):

        # try to interpolate everything using a linear interpolation first
        regridder = regridder_faces[face]
        grid = regridder(np.asfortranarray(L0_var_faces[face]).T).T

        # mask out any values which should be 0'd based on the old bathy
        grid[L0_wet_grid_on_L1_faces[face] == 0] = 0

        # mask out any values which should be 0'd based on the new bathy
        grid[L1_wet_grid_faces[face] == 0] = 0

        # plt.imshow(grid, origin='lower')
        # plt.show()

        # spread the the variable outward to new wet cells
        if spread_horizontally:
            grid, _ = spread_var_horizontally_in_wet_grid(grid, L1_wet_grid_faces[face])

        # C = plt.imshow(downscaled_grid,origin='lower',
        #                vmin=np.min(downscaled_grid[downscaled_grid!=0]),vmax=np.max(downscaled_grid[downscaled_grid!=0]))
        # plt.colorbar(C)
        # plt.show()

        grid_faces[face] = grid

    return (grid_faces)

def downscale_3D_field_faces(L0_var_faces,
                             regridder_faces, L0_wet_grid_on_L1_faces,
                             L1_wet_grid_faces,
                             level_row_cols_per_face,
                             spread_horizontally=True):

    grid_faces = {}

    for face in range(1,6):

        regridder = regridder_faces[face]
        full_grid = np.zeros(level_row_cols_per_face[face]).astype(float)

        for k in range(1): #np.shape(L1_wet_grid)[0]):

            # try to interpolate everything using a linear interpolation first
            grid = regridder(np.asfortranarray(L0_var_faces[face][k,:,:]).T).T

            # mask out any values which should be 0'd based on the old bathy
            grid[L0_wet_grid_on_L1_faces[face][k,:,:] == 0] = 0

            # mask out any values which should be 0'd based on the new bathy
            grid[L1_wet_grid_faces[face][k,:,:] == 0] = 0

            # plt.imshow(grid, origin='lower')
            # plt.show()

            # spread the the variable outward to new wet cells
            if spread_horizontally:
                grid, _ = spread_var_horizontally_in_wet_grid(grid, L1_wet_grid)

            # C = plt.imshow(downscaled_grid,origin='lower',
            #                vmin=np.min(downscaled_grid[downscaled_grid!=0]),vmax=np.max(downscaled_grid[downscaled_grid!=0]))
            # plt.colorbar(C)
            # plt.show()

            full_grid[k,:,:] = grid

        grid_faces[face] = full_grid

    return (grid_faces)

def downscale_3D_boundary_points(L0_points, L0_values, L0_wet_grid,
                                XC_subset, YC_subset, L1_wet_grid,
                                mean_vertical_difference=0,fill_downward=True,
                                printing=False,remove_zeros=True):

    # full_grid = np.zeros_like(L1_wet_grid).astype(float)
    full_grid = np.zeros((np.shape(L1_wet_grid)[0],np.shape(XC_subset)[0], np.shape(XC_subset)[1]))
    full_interp_mask = np.zeros((np.shape(L1_wet_grid)[0], np.shape(XC_subset)[0], np.shape(XC_subset)[1]))

    for k in range(np.shape(L0_wet_grid)[0]):

        only_fill_downward = False
        continue_to_interpolation = True

        if continue_to_interpolation:
            if printing:
                print('    Working on level ' + str(k) + ' of ' + str(np.shape(L0_wet_grid)[0])+' ('+str(np.sum(L1_wet_grid[k, :, :] > 0))+' nonzero points found)')
            if np.sum(L1_wet_grid[k, :, :] > 0)>0:
                # take an initial stab at the interpolation
                L0_points_subset = np.copy(L0_points)
                L0_values_subset = L0_values[k,:]
                L0_wet_grid_vert = np.reshape(L0_wet_grid[k, :], (np.size(L0_wet_grid[k, :]), ))
                L0_points_subset = L0_points_subset[L0_wet_grid_vert != 0, :]
                L0_values_subset = L0_values_subset[L0_wet_grid_vert != 0]
                if remove_zeros:
                    L0_points_subset = L0_points_subset[L0_values_subset != 0, :]
                    L0_values_subset = L0_values_subset[L0_values_subset != 0]

                # 1 means the point has been assigned a value (or doesn't need one)
                # 0 means the point still needs a value
                interp_mask = np.copy(1 - L1_wet_grid[k, :, :])

                if len(L0_points_subset) > 0:

                    # print(len(L0_points_subset))

                    if len(L0_points_subset) > 4:
                        grid = griddata(L0_points_subset, L0_values_subset, (XC_subset, YC_subset), method='linear',fill_value=np.nan)

                        # mask out any values which should be 0'd based on the new bathy
                        grid[L1_wet_grid[k, :, :] == 0] = 0

                        # mask out values which became nans by the interpolation
                        grid[np.isnan(grid)] = 0
                        interp_mask[np.isnan(grid)] = 0
                    else:
                        grid = np.zeros(np.shape(XC_subset))

                    # check if any points should be filled by a nearest neighbor
                    if np.any(interp_mask==0):
                        # grid = griddata(L0_points, L0_values, (XC_subset, YC_subset), method='linear',fill_value=0)
                        # grid = grid[:, :, 0]
                        grid_nearest = griddata(L0_points_subset, L0_values_subset, (XC_subset, YC_subset), method='nearest', fill_value=np.nan)
                        grid[interp_mask==0] = grid_nearest[interp_mask==0]
                        interp_mask[interp_mask==0] = 1
                else:
                    grid = np.zeros(np.shape(XC_subset))
                    # print('No points found for this layer')

            else:
                grid = np.zeros((np.shape(XC_subset)[0], np.shape(XC_subset)[1])).astype(float)
                interp_mask = np.copy(1-L1_wet_grid[k, :, :])

        full_grid[k, :, :] = grid
        full_interp_mask[k,:,:] = interp_mask

    # plt.subplot(1,2,1)
    # plt.imshow(L1_wet_grid[:,0,:])
    # plt.subplot(1,2,2)
    # plt.imshow(full_interp_mask[:,0,:])
    # plt.show()

    if np.any(full_interp_mask==0):
        for k in range(1,np.shape(full_interp_mask)[0]):
            if np.any(full_interp_mask[k,:,:]==0):
                rows,cols = np.where(full_interp_mask[k, :, :]==0)
                for ri in range(len(rows)):
                    full_grid[k,rows[ri],cols[ri]] = full_grid[k-1,rows[ri],cols[ri]]

    return(full_grid)


