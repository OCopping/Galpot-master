

import scipy as sci
import matplotlib.pyplot as plt
import astropy.units as u

#datfile = open('nperi_vs_mass.txt', 'r')

a = sci.genfromtxt('nperi_vs_mass.txt')

pids = a[:,0]
mw_masses = a[:,1]
rs_es = a[:,2]
lmc_masses = a[:,3]
infalls = a[:,4]
largest_radii = a[:,5]

print(len(mw_masses))

avg_nperi = sci.average(infalls)

print(avg_nperi)
