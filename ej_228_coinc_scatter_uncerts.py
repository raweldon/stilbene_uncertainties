'''
Note: This code was updated in the plastic final analysis on daq2

Full uncertainty accounting
Split into sytematic and propagated uncertainties -- sytematic order of magnitude larger than prop
'''
import numpy as np

def Ed_to_En(md, mn, m3He, Q, Ed):
    return ((np.sqrt(md*mn*Ed)/(mn + m3He)) + np.sqrt((md*mn*Ed)/(mn + m3He)**2 + (Ed*(m3He - md) + m3He*Q)/(mn + m3He)))**2

# measurements
#       BL 70                                                              BR 70
angles = [70., 60., 50., 40., 30., 20., 20., 30., 40., 50., 60., 70.]
dists = [66.4, 63.8, 63.2, 63.6, 64.9, 65.8, 66.0, 65.3, 63.4, 62.7, 64.3, 66.5]
Ens = [11.325, 4.825]
Ed = [8.209, 1.611] # from http://www.tunl.duke.edu/magnet.php , dcell length = 2.8575 cm
Ed_loss_in_gas = [0.095, 0.344] # from http://www.tunl.duke.edu/magnet.php (pressure: 11.325 - 1.85psi, 4.825 - 1.78 psi)
counts = ((165.0, 318.0, 375.0, 437.0, 430.0, 441.0, 416.0, 435.0, 420.0, 392.0, 335.0, 215.0),
          (140.0, 231.0, 278.0, 322.0, 318.0, 298.0, 285.0, 273.0, 288.0, 318.0, 233.0, 126.0)) # from /home/radians/raweldon/tunl.2018.1_analysis/plastic_analysis/final_analysis/analysis_lo/plastic_lo_vs_en.py

# systematic uncertainties
## position
crys_beam_alignment = 1 # cm
backing_det_pos = 0.5 # laser width = +-2mm, laser range finder center = +- 0.1, slop - rest 
support_align = 1 # backing det support alignment
det_interaction_location = np.sqrt(2.)/2. # corners of cube -- NOT A SYSTEMATIC
sigma_pos = np.sqrt(crys_beam_alignment**2 + backing_det_pos**2 + support_align**2)# + det_interaction_location**2)

# constants
md = 2.01410178    #D mass
mn = 1.008665      #neutron mass
m_he3 = 3.0160293  #He-3 mass
Q = 3.26891        #d(d,n) Q-value
sigma_Ed = 0.0005 # 500 eV (listed on tunl site)

for n, En in enumerate(Ens):
    Ep = [En*np.sin(np.deg2rad(ang))**2 for ang in angles]
    print '\n---------------------- En = ' + str(En) + ' ----------------------'
    print 'theta    sigma_sys    sigma_Ep    sigma_Ep_bar    sigma_tot    rel_uncert'
    ## systematics
    sigma_sys_degree = [np.arctan(sigma_pos/d) for d in dists]
    sigma_sys = [En*np.sin(np.deg2rad(2*theta))*sigma_sys_degree[i] for i, theta in enumerate(angles)]

    # Propagated uncertainties (E_p = E_n * cos^2(theta))
    #d_En_d_Ed = (md*mn/(mn + m_he3)**2 + (m_he3 - md)/(mn + m_he3))/(2*np.sqrt((Ed[n]*(m_he3 - md) + m_he3*Q)/(mn + m_he3) + mn*md*Ed[n]/(mn + m_he3)**2)
    #            ) + mn*md/(2*(mn + m_he3)*np.sqrt(mn*md*Ed[n]))  # from 8/20/2018 notes

    En_loss_in_gas = Ed_to_En(md, mn, m_he3, Q, Ed[n] + Ed_loss_in_gas[n]/2.) - Ed_to_En(md, mn, m_he3, Q, Ed[n] - Ed_loss_in_gas[n]/2)

    sigma_n_prod = En_loss_in_gas/np.sqrt(3.) # uniform distribution
    #simga_n_prod_ed = Ed_loss_in_gas[n]/np.sqrt(3.)

    ## uncert in neutron energy
    sigma_En = np.sqrt(sigma_Ed**2 + sigma_n_prod**2)

    ## uncert in scatter angle
    sigma_theta = [np.arctan((2.54)/d) for d in dists]

    sigma_Ep = [np.sqrt((np.sin(np.deg2rad(theta))**2*sigma_En)**2 + (En*np.sin(np.deg2rad(2*theta))*sigma_theta[i])**2) for i, theta in enumerate(angles)]
    #print sigma_Ep
    N = 100 # counts
    sigma_Ep_bar = [sigma_Ep[i]/np.sqrt(N) for i, N in enumerate(counts[n])]

    # total uncertainty
    sigma_tot = [np.sqrt(sigma_sys[i]**2 + sigma_Ep_bar[i]**2) for i, sig in enumerate(sigma_sys)]
    for s, sig in enumerate(sigma_tot):
        print '{:^7} {:>8} {:>11} {:>13} {:>14} {:>12}'.format(angles[s], round(sigma_sys[s],3), round(sigma_Ep[s],3), round(sigma_Ep_bar[s],3), round(sigma_tot[s],3), 
                                                               round(sigma_Ep_bar[s]/Ep[s],3))

    print En_loss_in_gas
    print sigma_n_prod
    print sigma_En
    print (np.sin(np.deg2rad(theta))**2*sigma_En)**2, (En*np.sin(np.deg2rad(2*theta))*sigma_theta[i])**2


