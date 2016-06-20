'''
GenerateAll.py creates all possible combinations of thicknesses
of Single Stack coatings and evaluate their performances by
iterating over thicknesses and computing VLT, TSER, L, a, b

general: good results at Ag: 12-15 nm, AlN(H): 20-70 nm

but to choose your own iteration boundaries, just change the while loop
conditions and variable value resets

**alteratively, uncomment the 'criteria' variables and the if condition at
the end to only output those thickness combinations that satisfy your criteria
(ex. TSER + Tvis > 1.20) or use the 'best_criteria' and output only the
best match
'''


import ga_2layers_aln_bright as bright

import os
import numpy as np
import opt_sim.nklib as nklib
import opt_sim.structure
from opt_sim.structure import MultiLayer as ML
import GenerateColor

def output(write_to, Tvis, TSER, L, a, b, GivenAlNHSeed, GivenAgBot, GivenAlNHBot):

    write_to.write(str(Tvis) + "," + str(TSER) + "," + str(L) + "," + str(a) +"," + str(b) + "," + str(GivenAlNHSeed.d) + "," + str(GivenAgBot.d) + "," + str(GivenAlNHBot.d))
    write_to.write("\n")

P_data=np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
S_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)

write_to = open("All_SingleStack_Results.txt", 'w')
write_to.write("VLT TSER L a b Seed Layer1 Layer2")


Ag_thickness = 8.5
AlNH_seed_thickness = 2.8
#best_criteria = 0.05

while AlNH_seed_thickness < 15:
    Ag_thickness = 8.5
    AlNH_seed_thickness += 0.2
    while Ag_thickness < 14:
        Ag_thickness += 0.2
        AlNH_thickness = 40
        while AlNH_thickness < 55:
            #Ag_thickness=random.uniform(13, 16)
            #AlN_thickness=random.uniform(20,55)
            AlNH_thickness += 0.2
            
            GivenAgBot = nklib.Ag(Ag_thickness)
            GivenAlNHBot = nklib.AlNH(AlNH_thickness)
            GivenAlNHSeed = nklib.AlNH(AlNH_seed_thickness)
            #GivenAgTop = nklib.Ag(18.9+3.5)
            #GivenAlNTop = nklib.AlN(60)

            MyStack = ML([GivenAlNHSeed, GivenAgBot, GivenAlNHBot])
            #MyStack = ML([GivenAgBot, GivenAlNBot, GivenAgTop, GivenAlNTop])
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
            #TSER = 1 - sum(MyStack.T[index_lowerUV:index_upper]*sol_array)/(sum(sol_array)) - 0.04
            TSER = sum(MyStack.R[index_lowerUV:index_upper]*sol_array)/(sum(sol_array))+0.85*sum(MyStack.A[index_lowerUV:index_upper]*sol_array)/(sum(sol_array))

        #print 'UV lower bound:', index_lowerUV
        #print 'Visible lower bound:', index_lowervis
        #print 'Visible upper bound:', index_uppervis
        #print 'IR upper bound:', index_upper
        #print MyStack.wl[:]

            #criteria = abs(Tvis - 0.77) + abs(TSER - 0.47)

            output(write_to, Tvis, TSER, L, a, b, GivenAlNHSeed, GivenAgBot, GivenAlNHBot)

            '''
            if (criteria < best_criteria):
                best_criteria = criteria
                print ("==============================================================")
                print (GivenAlNSeed, GivenAgBot, GivenAlNBot)
                #print (GivenAgBot, GivenAlNBot, GivenAgTop, GivenAlNTop) 
                print ('VLT = %s' % float("{0:.4f}".format(Tvis)))
                print ('TSER = %s' % float("{0:.4f}".format(TSER)))
                print 'Color in L* a* b* space is:', L, a, b

            ''' 
write_to.close()
