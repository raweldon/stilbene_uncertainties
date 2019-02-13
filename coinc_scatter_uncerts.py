'''
Full uncertainty accounting
Split into sytematic and propagated uncertainties -- sytematic order of magnitude larger than prop
'''
import numpy as np

# measurements
#       BL 70                                                              BR 70
angles = [70., 60., 50., 40., 30., 20., 20., 30., 40., 50., 60., 70.]
dists = [66.4, 63.8, 63.2, 63.6, 64.9, 65.8, 66.0, 65.3, 63.4, 62.7, 64.3, 66.5]
Ens = [11.325, 4.825]
Ed = [8.209, 1.611] # from http://www.tunl.duke.edu/magnet.php , dcell length = 2.8575
Ed_loss_in_gas = [0.095, 0.344] # from http://www.tunl.duke.edu/magnet.php (pressure: 11.325 - 1.85psi, 4.825 - 1.78 psi)

# systematic uncertainties
## position
crys_beam_alignment = 1 # cm
laser_range_finder_centered = 0.1 
backing_det_pos = 0.2 # laser width = 4mm
support_align = 1 # backing det support alignment
det_interaction_location = np.sqrt(2.)/2. # corners of cube
sigma_pos = np.sqrt(crys_beam_alignment**2 + laser_range_finder_centered**2 + backing_det_pos**2 + support_align**2 + det_interaction_location**2)

## angle
axes_marks = 1 # degree
housing_mount = 1
sigma_ang = np.sqrt(axes_marks**2 + housing_mount**2)


for n, En in enumerate(Ens):
    Ep = [En*np.sin(np.deg2rad(ang))**2 for ang in angles]
    print '\n---------------------- En = ' + str(En) + ' ----------------------'
    print 'theta    sigma_sys    sigma_Ep    sigma_Ep_bar    sigma_tot    rel_uncert'
    ## systematics
    sigma_sys_degree = [np.sqrt(np.arctan(sigma_pos/d)**2 + np.deg2rad(sigma_ang)**2) for d in dists]
    sigma_sys = [En*np.sin(np.deg2rad(2*theta))*sigma_sys_degree[i] for i, theta in enumerate(angles)]

    # Propagated uncertainties (E_p = E_n * cos^2(theta))
    ## uncert in neutron energy
    md = 2.01410178    #D mass
    mn = 1.008665      #neutron mass
    m_he3 = 3.0160293  #He-3 mass
    Q = 3.26891        #d(d,n) Q-value
    sigma_Ed = 0.0005 # 500 eV (listed on tunl site)
    d_En_d_Ed = (md*mn/(mn + m_he3)**2 + (m_he3 - md)/(mn + m_he3))/(2*np.sqrt((Ed[n]*(m_he3 - md) + m_he3*Q)/(mn + m_he3) + mn*md*Ed[n]/(mn + m_he3)**2)
                ) + mn*md/(2*(mn + m_he3)*np.sqrt(mn*md*Ed[n]))  # from 8/20/2018 notes

    sigma_n_prod = Ed_loss_in_gas[n]/np.sqrt(3.) # 103 keV loss through d cell, uniform distribution
    sigma_En = np.sqrt(sigma_Ed**2 + sigma_n_prod**2)

    ## uncert in scatter angle
    sigma_theta = [np.arctan((2.54)/d) for d in dists]

    sigma_Ep = [np.sqrt((np.sin(np.deg2rad(theta))**2*sigma_En)**2 + (En*np.sin(np.deg2rad(2*theta))*sigma_theta[i])**2) for i, theta in enumerate(angles)]
    #print sigma_Ep
    N = 100 # counts
    sigma_Ep_bar = sigma_Ep/np.sqrt(N)

    # total uncertainty
    sigma_tot = [np.sqrt(sigma_sys[i]**2 + sigma_Ep_bar[i]**2) for i, sig in enumerate(sigma_sys)]
    for s, sig in enumerate(sigma_tot):
        print '{:^7} {:>8} {:>11} {:>13} {:>14} {:>12}'.format(angles[s], round(sigma_sys[s],3), round(sigma_Ep[s],3), round(sigma_Ep_bar[s],3), round(sigma_tot[s],3), round(sigma_Ep_bar[s]/Ep[s],3))

    print sigma_En
    print (np.sin(np.deg2rad(theta))**2*sigma_En)**2, (En*np.sin(np.deg2rad(2*theta))*sigma_theta[i])**2


