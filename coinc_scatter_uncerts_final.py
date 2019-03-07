#!/usr/bin/env python
'''
Original code on desktop at /raweldon/Research/TUNL/git_programs/stilbene_uncertainties and on github
Full uncertainty accounting
Split into sytematic and propagated uncertainties -- sytematic order of magnitude larger than prop
'''
import numpy as np

def Ed_to_En(md, mn, m3He, Q, Ed):
    return ((np.sqrt(md*mn*Ed)/(mn + m3He)) + np.sqrt((md*mn*Ed)/(mn + m3He)**2 + (Ed*(m3He - md) + m3He*Q)/(mn + m3He)))**2

def calc_ep_uncerts(counts, beam_4MeV, print_unc):
    # measurements
    #       BL 70                                                    BR 70
    angles = [70., 60., 50., 40., 30., 20., 20., 30., 40., 50., 60., 70.]
    dists = [66.4, 63.8, 63.2, 63.6, 64.9, 65.8, 66.0, 65.3, 63.4, 62.7, 64.3, 66.5]
    if beam_4MeV:
        En = 4.825
        Ed = 1.611
        Ed_loss_in_gas = 0.344
    else:
        En = 11.325
        Ed = 8.209 # from http://www.tunl.duke.edu/magnet.php , dcell length = 2.8575 cm
        Ed_loss_in_gas = 0.095 # from http://www.tunl.duke.edu/magnet.php (pressure: 11.325 - 1.85psi, 4.825 - 1.78 psi)

    # systematic uncertainties
    ## position
    backing_det_pos = 0.3 # laser width = +-2mm, laser range finder center = +- 0.1, slop - rest
    support_align = 0.5 # backing det support alignment
    sigma_pos = np.sqrt(backing_det_pos**2 + support_align**2)# + det_interaction_location**2)
    sigma_en_sys = 0.0397 # from pandas_ej228_plots.py (final paper plots)

    # constants
    md = 2.01410178    #D mass
    mn = 1.008665      #neutron mass
    m_he3 = 3.0160293  #He-3 mass
    Q = 3.26891        #d(d,n) Q-value
    sigma_Ed = 0.0005 # 500 eV (listed on tunl site)

    Ep = [En*np.sin(np.deg2rad(ang))**2 for ang in angles]

    # calc uncerts
    if print_unc:
        print '\n---------------------- En = ' + str(En) + ' ----------------------'
        print 'theta    sigma_sys    sigma_Ep    sigma_Ep_bar    rel_uncert    sigma_tot'
    ## systematics
    sigma_pos_degree = [np.arctan(sigma_pos/d) for d in dists]
    sigma_theta_sys = [En*np.sin(np.deg2rad(2*theta))*sigma_pos_degree[i] for i, theta in enumerate(angles)]
    sigma_sys = [np.sqrt(sigma_en_sys**2 + sig**2) for sig in sigma_theta_sys]

    En_loss_in_gas = Ed_to_En(md, mn, m_he3, Q, Ed + Ed_loss_in_gas/2.) - Ed_to_En(md, mn, m_he3, Q, Ed - Ed_loss_in_gas/2)

    sigma_n_prod = En_loss_in_gas/np.sqrt(3.) # uniform distribution
    #simga_n_prod_ed = Ed_loss_in_gas[n]/np.sqrt(3.)

    ## uncert in neutron energy
    sigma_En = np.sqrt(sigma_Ed**2 + sigma_n_prod**2)

    ## uncert in scatter angle
    sigma_theta = [np.arctan((2.54)/(d + 2.54)) for d in dists]

    sigma_Ep = [np.sqrt((np.sin(np.deg2rad(theta))**2*sigma_En)**2 + (En*np.sin(np.deg2rad(2*theta))*sigma_theta[i])**2)
                for i, theta in enumerate(angles)]
    #print sigma_Ep
    N = 100 # counts
    sigma_Ep_bar = [sigma_Ep[i]/np.sqrt(N) for i, N in enumerate(counts)]

    # total uncertainty
    sigma_tot = [np.sqrt(sigma_sys[i]**2 + sigma_Ep_bar[i]**2) for i, sig in enumerate(sigma_sys)]
    if print_unc:
        for s, sig in enumerate(sigma_tot):
            print '{:^7} {:>8} {:>11} {:>13} {:>14} {:>12}'.format(angles[s], round(sigma_sys[s],3), round(sigma_Ep[s],3), round(sigma_Ep_bar[s],3),
                                                                   round(sigma_Ep_bar[s]/Ep[s],3), round(sigma_tot[s],3))
    return sigma_Ep_bar, sigma_tot

if __name__ == '__main__':
    calc_ep_uncerts([100]*12, beam_4MeV=True, print_unc=True)
    calc_ep_uncerts([100]*12, beam_4MeV=False, print_unc=True)