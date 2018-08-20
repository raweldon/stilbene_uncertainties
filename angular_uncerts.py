'''
    Systematic uncertainties from stilbene measurements at TUNL
    Relative uncertainties are large for small scatter angles (as expected)
    Numer to quote is the total angular uncertainty (sigma_tot)
'''
import numpy as np

# systematic uncertainties
# position
crys_beam_alignment = 1 # cm
laser_range_finder_centered = 0.1 
backing_det_pos = 0.2 # laser width = 4mm
support_align = 1 # backing det support alignment
det_interaction_location = np.sqrt(2.)/2. # corners of cube

sigma_pos = np.sqrt(crys_beam_alignment**2 + laser_range_finder_centered**2 + backing_det_pos**2 + support_align**2 + det_interaction_location**2)
print sigma_pos

# angle
axes_marks = 1 # degree
housing_mount = 1
#       BL 70                                                              BR 70
angles = [70., 60., 50., 40., 30., 20., 20., 30., 40., 50., 60., 70.]
dists = [66.4, 63.8, 63.2, 63.6, 64.9, 65.8, 66.0, 65.3, 63.4, 62.7, 64.3, 66.5]
det_size = [np.rad2deg(np.arctan((2.54+sigma_pos)/d)) for d in dists]
 
sigma_ang = [np.sqrt(axes_marks**2 + housing_mount**2 + d**2) for d in det_size]
sigma_ang2 = np.sqrt(axes_marks**2 + housing_mount**2)


# total uncert
En_11 = 11.325 # MeV
Ep_11 = [En_11 * np.sin(np.deg2rad(a))**2 for a in angles]
Ep_11_uncert = [En_11 * np.sin(2*np.deg2rad(a)) * np.deg2rad(sigma_ang[i]) for i, a in enumerate(angles)]

En_4 = 4.825
Ep_4 = [En_4 * np.sin(np.deg2rad(a))**2 for a in angles]
Ep_4_uncert = [En_4 * np.sin(2*np.deg2rad(a)) * np.deg2rad(sigma_ang[i]) for i, a in enumerate(angles)]

print [En_11 * np.sin(np.deg2rad(2*a)) * np.sqrt(np.tan(sigma_pos/dists[i])**2 + np.deg2rad(sigma_ang2)**2) for i, a in enumerate(angles)]

print '   angle (deg) |  uncert (deg) | uncert 11.3 (MeV) | rel uncert | uncert 4.8 (MeV) | rel uncert'
for i, sigma in enumerate(sigma_ang):
    if i == 0:
        print 'BL    ' + str(angles[i]) + '            ' + str(round(sigma,3)) + '            ' + str(round(Ep_11_uncert[i],3)
               ) + '           ' + str(round(Ep_11_uncert[i]/Ep_11[i],3)) + '           ' + str(round(Ep_4_uncert[i],3)
               ) + '           ' + str(round(Ep_4_uncert[i]/Ep_4[i],3))
    elif i == len(sigma_ang)-1:
        print 'BR    ' + str(angles[i]) + '            ' + str(round(sigma,3)) + '            ' + str(round(Ep_11_uncert[i],3)
               ) + '           ' + str(round(Ep_11_uncert[i]/Ep_11[i],3)) + '           ' + str(round(Ep_4_uncert[i],3)
               ) + '           ' + str(round(Ep_4_uncert[i]/Ep_4[i],3))
    else:
        print '      ' + str(angles[i]) + '            ' + str(round(sigma,3)) + '            ' + str(round(Ep_11_uncert[i],3)
               ) + '           ' + str(round(Ep_11_uncert[i]/Ep_11[i],3)) + '           ' + str(round(Ep_4_uncert[i],3)
               ) + '           ' + str(round(Ep_4_uncert[i]/Ep_4[i],3))


