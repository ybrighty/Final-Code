'''computes the color value given wavelength, T, R

can be made to read the data from a text file with wavelength / T / R

refer to https://windows.lbl.gov/software/window/Docs/WINDOW_4.0_Documentation_of_Calculation_Procedures.pdf
for calculation methods -- Bright Ye, yufeng.ye@mail.utoronto.ca
'''


import numpy as np
from opt_sim import calc
import os

#function used in conversion from CIEXYZ to CIELAB

def f(t):
        if t > (6/29)**3:
                return t**(1/float(3))
        else:
                return 1/float(3)*((29/float(6))**2)*t+4/float(29)


def color_calc(wl,T,R):
#        os.chdir(os.getcwd() + "\\plot support files")
        S_data = np.loadtxt("plot support files\\Standard_Illuminant_D65.txt")
        xyzdata = np.loadtxt("plot support files\\(photopic)ciexyz64.txt")

        s_wl = S_data[:,0]
        s_array = S_data[:,1]

        xyzwl = xyzdata[:,0]
        d_wl = (xyzwl[-1] - xyzwl[0]) / len(xyzwl)


        if wl[0] > xyzwl[0] or wl[-1] < xyzwl[-1]:
                raise Exception("The entire visible spectrum must be available in the structure simulation.")

        # Future person beware!! T & R must be given in %'s (i.e., 1-100 as supposed to 0.01 - 1)! -- Bright
        if T[int(len(T)/2)] < 1:
                T *= 100
                R *= 100

        xyz_sens = xyzdata[:,1:]

        T_interp = calc.interpolate(wl, T, xyzwl)
        R_interp = calc.interpolate(wl, R, xyzwl)
        S_interp = calc.interpolate(s_wl, s_array, xyzwl)
        

        # initializing
        T_xyz = [0, 0, 0]
        R_xyz = [0, 0, 0]
        
        for i in range(3):
                T_xyz[i] = np.trapz(T_interp * S_interp * xyz_sens[:,i], xyzwl) / float(np.trapz(S_interp * xyz_sens[:,1], xyzwl))
                R_xyz[i] = np.trapz(R_interp * S_interp * xyz_sens[:,i], xyzwl) / float(np.trapz(S_interp * xyz_sens[:,1], xyzwl))
        
        
        X,Y,Z = T_xyz[:]
        Xn, Yn, Zn = 95.047, 100.000, 108.883
        L = 116*f(Y/Yn)-16
        a = 500*(f(X/Xn)-f(Y/Yn)) + 0.4124  # this is how much "a" is off by, when we give 100% T or R -- Bright
        b = 200*(f(Y/Yn)-f(Z/Zn)) - 0.9589 # this is how much "a" is off by, when we give 100% T or R -- Bright

        # reset T / R to decimals (i.e. back to 0.01 - 1 as supposed to 1-100)
        if T[int(len(T)/2)] > 1:
                T /= 100
                R /= 100


        return L,a,b

        '''
        print T_color
        print R_color
        '''

'''
#uncomment below when not calling but directly using this file to compute color
input_data=np.loadtxt("plot support files\\SCC_ds_5.txt", skiprows=2)
wl = input_data[:,0]
R = input_data[:,1] / 100
T = input_data[:,2] / 100

### Use the code below if need to recompute TSER & VLT to match LBNL output
P_data=np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
S_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)

A = [1] * len(T) - T - R
index_lowerUV = 0
index_lowervis = 32
index_uppervis = 94
index_upper = 454
photopic_array = np.interp(wl[index_lowervis:index_uppervis],P_data[:, 0] , P_data[:, 1])
sol_array = np.interp(wl[index_lowerUV:index_upper],S_data[:, 0] , S_data[:, 3])
Tvis = sum(T[index_lowervis:index_uppervis]*photopic_array)/(sum(photopic_array))
#no -0.04
#TSER = 1 - sum(T[index_lowerUV:index_upper]*sol_array)/(sum(sol_array)) 
TSER = sum(R[index_lowerUV:index_upper]*sol_array)/(sum(sol_array))+0.85*sum(A[index_lowerUV:index_upper]*sol_array)/(sum(sol_array))
print ('VLT = %s' % float("{0:.4f}".format(Tvis)))
print ('TSER = %s' % float("{0:.4f}".format(TSER)))



print(color_calc(wl,T,R))
'''
