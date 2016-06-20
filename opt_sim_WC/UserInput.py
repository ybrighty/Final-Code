'''
UserInput.py generates the performances of a coating given 
thickness parameters by the user

edit values like 'AlN_seed = x' to change thickness parameters,
for double stack / inclusion of seeding layers, uncomment
lines like 'Ag_top_thickness = y' and include the corresponding
'GivenAgTop = nklib.Ag(Ag_top_thickness)' and
'MyStack = ML([GivenSeed, GivenAgBot, GivenAlNBot, GivenAgTop, GivenAlNTop])'
'''

import ga_2layers_aln_bright as bright
import opt_sim as opt
import numpy as np
import opt_sim.nklib as nklib
import opt_sim.structure
from opt_sim.structure import MultiLayer as ML
import matplotlib.pyplot as plt
import GenerateColor 

P_data=np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
S_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)
#S_data = np.loadtxt("plot support files\\Standard_Illuminant_D65.txt")

AlN_seed = 3
Ag_bot_thickness = 10 # 12.2 (-3.5) originally for scc 10
AlN_bot_thickness = 43 # 54.2 originally for scc 10
#Ag_top_thickness = 14
#AlN_top_thickness = 45.1

GivenSeed = nklib.AlN(AlN_seed)
GivenAgBot = nklib.Ag(Ag_bot_thickness)
GivenAlNBot = nklib.AlN(AlN_bot_thickness)

#GivenAgTop = nklib.Ag(Ag_top_thickness)
#GivenAlNTop = nklib.AlN(AlN_top_thickness)

#MyStack = ML([GivenAgBot, GivenAlNHBot])
MyStack = ML([GivenSeed, GivenAgBot, GivenAlNBot])
#MyStack = ML([GivenSeed, GivenAgBot, GivenAlNBot, GivenAgTop, GivenAlNTop])
MyStack.calculate_TR()
MyStack = bright.correct_TR(MyStack)
MyStack.A = [1] * len(MyStack.T) - MyStack.T - MyStack.R

L,a,b = GenerateColor.color_calc(MyStack.wl, MyStack.T, MyStack.R)

index_lowerUV = 0
index_lowervis = 32
index_uppervis = 94
index_upper = 454

photopic_array = np.interp(MyStack.wl[index_lowervis:index_uppervis],P_data[:, 0] , P_data[:, 1])
sol_array = np.interp(MyStack.wl[index_lowerUV:index_upper],S_data[:, 0] , S_data[:, 3])

Tvis = sum(MyStack.T[index_lowervis:index_uppervis]*photopic_array)/(sum(photopic_array))
#yes -0.04
#TSER = 1 - sum(MyStack.T[index_lowerUV:index_upper]*sol_array)/(sum(sol_array)) - 0.04
TSER = sum(MyStack.R[index_lowerUV:index_upper]*sol_array)/(sum(sol_array))+0.85*sum(MyStack.A[index_lowerUV:index_upper]*sol_array)/(sum(sol_array))

#print 'UV lower bound:', MyStack.wl[index_lowerUV]
#print 'Visible lower bound:', MyStack.wl[index_lowervis]
#print 'Visible upper bound:', MyStack.wl[index_uppervis]
#print 'IR upper bound:', MyStack.wl[index_upper]
#print MyStack.wl[:]


print ("==============================================================")
#print (GivenAgBot, GivenAlNBot)
print (GivenSeed, GivenAgBot, GivenAlNBot)
#print (GivenSeed, GivenAgBot, GivenAlNBot, GivenAgTop, GivenAlNTop) 
print ('VLT = %s' % float("{0:.4f}".format(Tvis)))
print ('TSER = %s' % float("{0:.4f}".format(TSER)))
#for i in range(len(MyStack.wl)):
#    print MyStack.wl[i],",", MyStack.T[i],",", MyStack.R[i]
#print ('Expected VLT = %s' % float(75.6))
#print ('Expected TSER = %s' % float(46))
print 'Color in L* a* b* space is:', L, a, b
opt.plot.TR([MyStack])
plt.show()

