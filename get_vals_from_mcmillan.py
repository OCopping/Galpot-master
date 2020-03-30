import numpy as np
import subprocess

np.random.seed(732)
mcmillan = np.genfromtxt('ProductionRunBig2.tab')
#m200s = mcmillan[1:,43]

##m200_file = open('m200_vals','w')
#np.savetxt('m200_vals',m200s)

mcmillan_rows  = np.random.choice(len(mcmillan),1000,replace=False)
print("picked rows:", len(mcmillan_rows))

n_threads = 8
num_rows = len(mcmillan_rows) / n_threads
print(num_rows)

for i in range(n_threads):
    #subprocess.run(['mkdir', 'thread_{}'.format(i)])
    subprocess.run('touch', 'thread_{0}/data_{1}'.format(i,num_rows*(i+1)))
    
