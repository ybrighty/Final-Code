import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import structure as st
import calc
import os


ds_data = np.loadtxt("plot support files\\Ag-AlN Double SCC 9.txt", skiprows=7)
ss_data = np.loadtxt("plot support files\\Ag-H-AlN Single SCC 10.txt", skiprows=7)
    
wl = np.array(ds_data[:,0])
T = np.array(ds_data[:,1])
R = np.array(ds_data[:,2])
#for i in range(len(wl)):
 #   print wl[i],T[i],R[i]
plt.plot(wl, T, label='T')
plt.plot(wl, R, label='R')
AM1p5_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)
solar_wl = AM1p5_data[:,0]
solar_intens = AM1p5_data[:,3]
solar_intens /= max(solar_intens)
plt.plot(solar_wl, solar_intens, label="AM1.5", color=(0.55,0.55,0.55))
#photopic_data = np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
#photopic_wl = photopic_data[:,0]
#photopic_intens = photopic_data[:,1]
#plt.plot(photopic_wl, photopic_intens, label="photopic", color=(1,0,0))
plt.title('Reflectance and Transmittance plots')
plt.xlabel('Wavelength ({})'.format('nm'))
plt.xlim(xmin=0)
plt.xlim(xmax=4000)
plt.ylim(0, 1)
plt.legend(loc=0)

wl = ss_data[:,0]
T = ss_data[:,1]
R = ss_data[:2]
plt.plot(wl, T, label='T')
plt.plot(wl, R, label='R')
AM1p5_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)
solar_wl = AM1p5_data[:,0]
solar_intens = AM1p5_data[:,3]
solar_intens /= max(solar_intens)
plt.plot(solar_wl, solar_intens, label="AM1.5", color=(0.55,0.55,0.55))
#photopic_data = np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
#photopic_wl = photopic_data[:,0]
#photopic_intens = photopic_data[:,1]
plt.plot(photopic_wl, photopic_intens, label="photopic", color=(1,0,0))
plt.title('Reflectance and Transmittance plots')
plt.xlabel('Wavelength ({})'.format('nm'))
plt.xlim(xmin=0)
plt.xlim(xmax=4000)
plt.ylim(0, 1)
plt.legend(loc=0)
