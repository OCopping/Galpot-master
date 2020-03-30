import numpy as np
import matplotlib.pyplot as plt
import os

import astropy.units as u
from astropy.coordinates import SkyCoord

import numpy

import scipy.optimize as sopt

import gal_uvw

#import vcirc

import scipy.interpolate as sinterp

import astropy.io.fits as afits
import scipy

M_NFW = 80.
q_NFW = 1.
rs_NFW = 16.

phi1_prog = 0.
phi2_prog = 0.

mu_phi2_prog = 0.

dist_prog = 19.1
rv_prog = 50.

iteration = 0

R_phi12_radec = np.array([[0.5964467, 0.27151332, -0.75533559],
                         [-0.48595429, -0.62682316, -0.60904938],
                         [0.63882686, -0.73032406, 0.24192354]])

a_g = np.array([[-0.0548755604, +0.4941094279, -0.8676661490],
                [-0.8734370902, -0.4448296300, -0.1980763734], 
                [-0.4838350155, 0.7469822445, +0.4559837762]])

G = 43007.105731706317

l_lmc = 280.4652*np.pi/180.
b_lmc = -32.8884*np.pi/180.

# Here we draw the potential values from McMillan's posterior distribution

np.random.seed(732)
mcmillan = np.genfromtxt('ProductionRunBig2.tab')
m200s = mcmillan[1:,43]

#m200_file = open('m200_vals','w')
np.savetxt('m200_vals',m200s)

mcmillan_rows  = np.random.choice(len(mcmillan),1000,replace=False)
print("picked rows:", len(mcmillan_rows))

# chi2_eval is the main function
def chi2_eval(pid):
    print('[Run {}]'.format(pid))
    
    for i in range(100):
        # generate observable parameters for the LMC
        mu_alpha_lmc = np.random.normal(1.91,0.02)
        mu_delta_lmc = np.random.normal(0.229,0.047)
        rv_lmc = np.random.normal(262.2,3.4)
        dist_lmc = np.random.normal(49590.,540.)

        # this loops over the mass of the LMC between 5e10 Msun and 2.5e11 Msun
        for M_LMC in range(5,30,5):


            if M_LMC > 2.:
                rs_LMC = np.sqrt(G*M_LMC*8.7/91.7**2.)-8.7
            else:
                rs_LMC = np.sqrt(G*2.*8.7/91.7**2.)-8.7

            tmax = 13.8
            
            global iteration
            iteration = 0
            lhood = chi2_worker(mu_alpha_lmc, mu_delta_lmc, rv_lmc, dist_lmc, M_LMC,rs_LMC, tmax, pid)
        
    return lhood
def chi2_worker(mu_alpha_lmc, mu_delta_lmc, rv_lmc, dist_lmc, M_LMC, rs_LMC, tmax, pid):

    global iteration
    iteration += 1
    #print(iteration)

    # take 1 row from McMillan's potential
    mcm = mcmillan[mcmillan_rows[pid]]

    pid_pot = pid
    
    pid = 0
    
    # take this row and generate a potential file from it.
    fileout = open('pot/PJM17_{0}.Tpot'.format(pid),'w')

    print(4, file=fileout)
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4}'.format(mcm[0],mcm[1],mcm[2],mcm[3],mcm[4]), file=fileout) # thin disk
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4}'.format(mcm[5],mcm[6],mcm[7],mcm[8],mcm[9]), file=fileout) # thick disk
    print('5.31319e+07 7 -0.085 4 0', file=fileout) # HI disk
    print('2.17995e+09 1.5 -0.045 12 0', file=fileout) # H2 disk
    
    print(2, file=fileout)
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4} {5}'.format(mcm[20],mcm[21],mcm[22],mcm[23],mcm[24],mcm[25]), file=fileout) # Bulge
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4} {5}'.format(mcm[26],mcm[27],mcm[28],mcm[29],mcm[30],mcm[31]), file=fileout) # Halo

    Usun = mcm[32]*1000. # Motion of the sun in x in km/s
    Vsun = mcm[33]*1000. # Motion of the sun in y in km/s
    Wsun = mcm[34]*1000. # Motion of the sun in z in km/s
    R0 = mcm[-6] # The locaton of the sun in kpc
    V0 = mcm[-5] # The circular velocity at the Sun's location

    M200 = 4.*np.pi*mcm[26]*mcm[30]**3.*(np.log(1.+mcm[-9]/mcm[30])-mcm[-9]/(mcm[-9]+mcm[30]))/(1.e10) # This the mass of the NFW

    c200 = mcm[-9]/mcm[30] # Concentration of the halo
    rs = mcm[30] # Scale radius of the halo

    fileout.close()

    # Now we need the LMC's 3d position and velocity

    # This is the 3d velocity of the Sun
    vlsr = np.array([Usun,Vsun+V0,Wsun])
    
    k_mu = 4.74047

    uvw_stationary = -vlsr

    # This takes the LMC's angle on the sky
    gc = SkyCoord(b=b_lmc*u.radian,l=l_lmc*u.radian,frame='galactic')

    # This line computes the LMC's x,y,z position
    x_lmc,y_lmc,z_lmc = np.array([-R0,0.,0.])+dist_lmc/1000.*np.array([np.cos(l_lmc)*np.cos(b_lmc),np.sin(l_lmc)*np.cos(b_lmc),np.sin(b_lmc)])

    # This line computes the LMC's vx,vy,vz (https://github.com/segasai/astrolibpy)
    vx_lmc,vy_lmc,vz_lmc = gal_uvw.gal_uvw(distance=dist_lmc,ra=np.array(gc.icrs.ra),dec=np.array(gc.icrs.dec),lsr=np.array([-Usun,0,Wsun]),pmra=mu_alpha_lmc,pmdec=mu_delta_lmc,vrad=rv_lmc)


    # These two lines fix up the velocity since gal_uvw has slightly different coordinates
    vy_lmc += Vsun+V0 
    vx_lmc = -vx_lmc

    os.system('./a.out {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}'.format(M_LMC,rs_LMC,x_lmc,y_lmc,z_lmc,vx_lmc,vy_lmc,vz_lmc,tmax,rs,M200,c200,pid))

    # reads two values: no. of in-falls and max distance of LMC from MW
    data = np.genfromtxt('orbit_{0}.txt'.format(pid))

    fileout = open('nperi_vs_mass.txt','a')

    print('{0} {1} {2} {3} {4} {5}'.format(pid_pot, M200,rs,M_LMC,data[0],data[1]), file=fileout)

    fileout.close()
    
    return 0

for pid in range(10):
    chi2_eval(pid)

