from opt_sim.structure import MultiLayer as ML
from opt_sim import *
import random
from Queue import PriorityQueue
from scipy import stats
import GenerateColor

elements=PriorityQueue() 

populationSize=80

numberOfLoser=60 # The number of population dying after each generation

new_population=100


P_data=np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
S_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)
# V_data is visible of AM1.5
V_data = np.loadtxt("plot support files\\AM1.5_Vis.txt")
# U_data is UV of AM1.5
U_data = np.loadtxt("plot support files\\AM1.5_UV.txt", skiprows=2)
# I data is IR of AM1.5
I_data = np.loadtxt("plot support files\\AM1.5_IR.txt")



class constants(object):
    '''
Below, lower_limit and upper_limit arguments are used to provide the range of wavelength
for optimization. This values are particularly important when optimizing for a certain range of
wavelengths on the photopic curve. The default range is 390-700. Note: the upper_limit value should
always be less than 700 and that of the lower_limit should be greater than 390.

T_factor and R_factor coresspond to the percentage of importance of optimizing for transmittance
of VL and reflectance of Infrared, respectively.

max_DLC_thickness coressponds to the maximum thickness that can be assigned for a given DLC layer. 
'''
    def __init__(self,lower_limit, upper_limit, p1, p2, p3, wanted_TSER, max_DLC_thickness):
        self.lower_limit=lower_limit
        self.upper_limit=upper_limit
        self.p1=p1 #priority constant for VLT
        self.p2=p2 #priority constant for TSER
        self.p3=p3 #priority constant for color (L,a,b)
        self.wanted_TSER=wanted_TSER
        self.max_DLC_thickness=max_DLC_thickness


myConst=constants(lower_limit=390, upper_limit=700,
                  p1=12, p2=1,p3=0.9, wanted_TSER=.55, max_DLC_thickness=150.0)


## FOR correct_TR()
glass_ext = np.loadtxt("glass_extinct_assume_n1.5.txt")
##wl_A, A = glass_ext[2:-25,0], glass_ext[2:-25,1]    #210nm to 2475nm Ag/DLC
wl_A, A = glass_ext[6:-20,0], glass_ext[6:-20,1]    #230nm to 2500nm Ag/DLC/AlN
##wl_A, A = glass_ext[46:-214,0], glass_ext[46:-214,1]  #430nm to 1530nm Ag/TiO2
MIN_WL = 230
MAX_WL = 2500
Tg = 0.96
Rg = 0.04
    
def correct_TR(struct, incidence="top"):
    if not struct._TR_calculated:
        struct.calculate_TR()
        
    # Commented out because right now wl_A and struct.wl
    # have the same range and step
    T0 = calc.interpolate(struct.wl, struct.T, wl_A)
    R0 = calc.interpolate(struct.wl, struct.R, wl_A)
##    T0 = struct.T
##    R0 = struct.R

    _layers_list = struct._layers_list[:]
    _layers_list.reverse()
    struct_rev = ML(_layers_list, unit=struct.unit,
                    label=struct.label+"_rev", min_wl=struct._min_wl,
                    max_wl = struct._max_wl, wl_step=struct._wl_step,
                    n0=struct.ns, ns=struct.n0)
    struct_rev.calculate_TR()

    Tr0 = calc.interpolate(struct_rev.wl, struct_rev.T, wl_A)
    Rr0 = calc.interpolate(struct_rev.wl, struct_rev.R, wl_A)
##    Tr0 = struct_rev.T
##    Rr0 = struct_rev.R

    if incidence == "top":
        T = (1 - A) * Tg * T0 / (1 - (1 - A)**2 * Rg * Rr0)
        R = R0 + (1 - A)**2 * Rg * T0 * Tr0 / (1 - (1 - A)**2 * Rg * Rr0)
    elif incidence == "bot":
        T = (1 - A) * Tg * Tr0 / (1 - (1 - A)**2 * Rg * Rr0)
        R = Rg + (1 - A)**2 * Rr0 * Tg * Tg / (1 - (1 - A)**2 * Rg * Rr0)

    struct_rev.wl = wl_A
    struct_rev.T = T
    struct_rev.R = R

    return struct_rev


def create_structure():
    '''
  This function creates a structure composed of randomly generated thicknesses.
  It returns a tuple with a list of layer thicknesses and priority value. 
  
 '''
    L0_thickness = 15
    seedlayer = AlN(L0_thickness)

    L1_thickness=random.uniform(12, 20.)
    L1_thickness=float("{0:.1f}".format(L1_thickness))
    layer_one=Ag(L1_thickness)
        
    L2_thickness=random.uniform(20.0, myConst.max_DLC_thickness)
    L2_thickness=float("{0:.1f}".format(L2_thickness))
    layer_two=AlNH(L2_thickness)
        
    L3_thickness=random.uniform(12, 20.)
    L3_thickness=float("{0:.1f}".format(L3_thickness))
    layer_three=Ag(L3_thickness)
        
    L4_thickness=random.uniform(20.0, myConst.max_DLC_thickness)
##    L4_thickness=random.uniform(20.0, 60)
    L4_thickness=float("{0:.1f}".format(L4_thickness))
    layer_four=AlNH(L4_thickness)
        
    mystruct = ML([seedlayer, layer_one, layer_two, layer_three, layer_four],
                  min_wl=MIN_WL, max_wl=MAX_WL)
    ML.calculate_TR(mystruct)

    # Doing unchecked attribute mutation here...future person beware. This
    # only applies because of comment above in correct_TR
    mystruct = correct_TR(mystruct)

    index_lowerUV=0
    index_lowervis=0
    index_uppervis=0
    upper_middle=0
    index_upper=0

    if mystruct.wl[len(mystruct.wl)-1]<3000:
        index_upper=len(mystruct.wl)-1
    for index,  i in enumerate(mystruct.wl):
        if index_lowervis ==0 and i>=myConst.lower_limit:
            index_lowervis = index
        elif index_uppervis==0 and i>=myConst.upper_limit:
            index_uppervis=index
        elif myConst.upper_limit!= 700 and upper_middle==0:
            if i>=700:
                upper_middle=index
        elif index_upper==0 and i>=3000:
            index_upper=index
            break
    if myConst.upper_limit==700:
        upper_middle=index_uppervis
        
    

        photopic_array = np.interp(mystruct.wl[index_lowervis:index_uppervis],P_data[:, 0] , P_data[:, 1])
    V_normalized_data = np.interp(mystruct.wl[index_lowervis:index_uppervis],V_data[:, 0], V_data[:, 3]) / max(V_data[:, 3])
    T_ideal_array = photopic_array/V_normalized_data
    #T_one_array = [1]*len(mystruct.wl[index_lowervis:index_uppervis])
	
    sol_array = np.interp(mystruct.wl[index_lowerUV:index_upper],S_data[:, 0] , S_data[:, 3])
    sol_array_UV = np.interp(mystruct.wl[index_lowerUV:index_lowervis],U_data[:, 0] , U_data[:, 3])
    sol_array_vis = np.interp(mystruct.wl[index_lowervis:index_uppervis],V_data[:, 0] , V_data[:, 3])
    sol_array_IR = np.interp(mystruct.wl[index_uppervis:index_upper],I_data[:, 0] , I_data[:, 3])
    
    Tvis = sum(mystruct.T[index_lowervis:index_uppervis]*photopic_array)/(sum(photopic_array))
   
    #Tpriority = sum(mystruct.T[index_lowervis:index_uppervis]*T_ideal_array)/(sum(T_ideal_array))
    #Tpriority = sum(mystruct.T[index_lowervis:index_uppervis]*sol_array_vis*photopic_array)/(sum(sol_array_vis*photopic_array))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]*sol_array_vis/max(sol_array_vis)-photopic_array)**2))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]-sol_array_vis/max(sol_array_vis))**2))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]-T_one_array)**2))
    Tpriority = Tvis
    
    TSER_UV = 1 - sum(mystruct.T[index_lowerUV:index_lowervis]*sol_array_UV)/(sum(sol_array_UV))
    TSER_vis = 1 - sum(mystruct.T[index_lowervis:index_uppervis]*sol_array_vis)/(sum(sol_array_vis))
    TSER_IR = 1 - sum(mystruct.T[index_uppervis:index_upper]*sol_array_IR)/(sum(sol_array_IR))
    # 0.0342 = % of AM1.5 energy in UV
    # 0.426 = % of AM1.5 energy in visible
    # 0.54 = % of AM1.5 energy in IR (700-2500)
    #TSER = TSER_UV * 0.0342 + TSER_vis * 0.426 + TSER_IR * 0.54 - .04
    TSER = 1 - sum(mystruct.T[index_lowerUV:index_upper]*sol_array)/(sum(sol_array)) - 0.04
	
    #TSER_priority = TSER_UV * 0.0342 / (0.0342 + 0.54) + TSER_IR * 0.54 / (0.0342 + 0.54)
    TSER_priority = TSER_IR

    L,a,b = GenerateColor.color_calc(mystruct.wl, mystruct.T, mystruct.R)
    Cpriority = 1/ (abs(a) + abs(b))
	
    #priority = myConst.p1 * Tpriority + myConst.p2 * (TSER_IR - myConst.wanted_TSER)
    priority = myConst.p1 * Tpriority + myConst.p2 * TSER_priority + myConst.p3 * Cpriority 
    #priority = myConst.p1 * Tpriority * Tvis + myConst.p2 * TSER_priority
    #priority = -1 / (Tvis - myConst.p1) - \
    #           abs(myConst.p2 * (TSER - myConst.wanted_TSER)) 				  
													  
    return (priority, L0_thickness, L1_thickness, L2_thickness, L3_thickness, L4_thickness,
            Tvis, TSER)
        
def create_inital_population():
    '''
  The following function creates the inital population containing 80 elements in the priority queue.
  It determines the thicknesses based on random choices in range 20 to 300nm for DLC layers and from
  3 to 20nm for silver layers. The numbers selected randomly get rounded to one decimal place.
 '''
    counter=0

    print 'Creating the inital population'

    while counter < populationSize:

        elements.put(create_structure())
        counter+=1

def recalculate_priority(arg):
    '''
This function returns a tuple containing a list of layers with calculated priority.
The structure is given as an argument
'''
    seed=AlN(arg[1])
    first=Ag(arg[2])
    second=AlNH(arg[3])
    third=Ag(arg[4])
    fourth=AlNH(arg[5])
    mystruct = ML([seed, first, second, third, fourth],
                  min_wl=MIN_WL, max_wl=MAX_WL)
    ML.calculate_TR(mystruct)

    # Doing unchecked attribute mutation here...future person beware. This
    # only applies because of comment above in correct_TR
    mystruct = correct_TR(mystruct)
            
    index_lowerUV=0
    index_lowervis=0
    index_uppervis=0
    upper_middle=0
    index_upper=0

    if mystruct.wl[len(mystruct.wl)-1]<3000:
        index_upper=len(mystruct.wl)-1
    for index,  i in enumerate(mystruct.wl):
        if index_lowervis ==0 and i>=myConst.lower_limit:
            index_lowervis = index
        elif index_uppervis==0 and i>=myConst.upper_limit:
            index_uppervis=index
        elif myConst.upper_limit!= 700 and upper_middle==0:
            if i>=700:
                upper_middle=index
        elif index_upper==0 and i>=3000:
            index_upper=index
            break
    if myConst.upper_limit==700:
        upper_middle=index_uppervis
        

        photopic_array = np.interp(mystruct.wl[index_lowervis:index_uppervis],P_data[:, 0] , P_data[:, 1])
    V_normalized_data = np.interp(mystruct.wl[index_lowervis:index_uppervis],V_data[:, 0], V_data[:, 3]) / max(V_data[:, 3])
    T_ideal_array = photopic_array/V_normalized_data
    #T_one_array = [1]*len(mystruct.wl[index_lowervis:index_uppervis])
	
    sol_array = np.interp(mystruct.wl[index_lowerUV:index_upper],S_data[:, 0] , S_data[:, 3])
    sol_array_UV = np.interp(mystruct.wl[index_lowerUV:index_lowervis],U_data[:, 0] , U_data[:, 3])
    sol_array_vis = np.interp(mystruct.wl[index_lowervis:index_uppervis],V_data[:, 0] , V_data[:, 3])
    sol_array_IR = np.interp(mystruct.wl[index_uppervis:index_upper],I_data[:, 0] , I_data[:, 3])
    
    Tvis = sum(mystruct.T[index_lowervis:index_uppervis]*photopic_array)/(sum(photopic_array))
   
    #Tpriority = sum(mystruct.T[index_lowervis:index_uppervis]*T_ideal_array)/(sum(T_ideal_array))
    #Tpriority = sum(mystruct.T[index_lowervis:index_uppervis]*sol_array_vis*photopic_array)/(sum(sol_array_vis*photopic_array))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]*sol_array_vis/max(sol_array_vis)-photopic_array)**2))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]-sol_array_vis/max(sol_array_vis))**2))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]-T_one_array)**2))
    Tpriority = Tvis
    
    TSER_UV = 1 - sum(mystruct.T[index_lowerUV:index_lowervis]*sol_array_UV)/(sum(sol_array_UV))
    TSER_vis = 1 - sum(mystruct.T[index_lowervis:index_uppervis]*sol_array_vis)/(sum(sol_array_vis))
    TSER_IR = 1 - sum(mystruct.T[index_uppervis:index_upper]*sol_array_IR)/(sum(sol_array_IR))
    # 0.0342 = % of AM1.5 energy in UV
    # 0.426 = % of AM1.5 energy in visible
    # 0.54 = % of AM1.5 energy in IR (700-2500)
    #TSER = TSER_UV * 0.0342 + TSER_vis * 0.426 + TSER_IR * 0.54 - .04
    TSER = 1 - sum(mystruct.T[index_lowerUV:index_upper]*sol_array)/(sum(sol_array)) - 0.04
	
    #TSER_priority = TSER_UV * 0.0342 / (0.0342 + 0.54) + TSER_IR * 0.54 / (0.0342 + 0.54)
    TSER_priority = TSER_IR

    L,a,b = GenerateColor.color_calc(mystruct.wl, mystruct.T, mystruct.R)
    Cpriority = 1/ (abs(a) + abs(b))
	
    #priority = myConst.p1 * Tpriority + myConst.p2 * (TSER_IR - myConst.wanted_TSER)
    priority = myConst.p1 * Tpriority + myConst.p2 * TSER_priority + myConst.p3 * Cpriority 
    #priority = myConst.p1 * Tpriority * Tvis + myConst.p2 * TSER_priority
    #priority = -1 / (Tvis - myConst.p1) - \
    #           abs(myConst.p2 * (TSER - myConst.wanted_TSER))                            
    
    return (priority, arg[1], arg[2], arg[3], arg[4], arg[5], Tvis, TSER)

#the fitness is not recalculated 
def crossover_and_mutation(parent_one_, parent_two_):
    '''
 This function combines the selected parents to create two children and then,
 interchage the position of some elements to create two more new children.
 Finally, children get added to the prirority queue.

'''

    parent_one=elements.queue[parent_one_]  
    parent_two=elements.queue[parent_two_]
    
    LAg_thickness=random.uniform(12,20.0)
    LAgtop_thickness=random.uniform(11,20)

    LAg=float("{0:.1f}".format(LAg_thickness))
    LAgtop=float("{0:.1f}".format(LAgtop_thickness))
        
    LAlNH_thickness=random.uniform(20.0, myConst.max_DLC_thickness)
    LAlNH=float("{0:.1f}".format(LAlNH_thickness))

    # variables on RHS has to be tuples, so even for single layers you have
    # to index something like [3:4] rather than [3] (this is terrible I know;
    # I didn't write this -Remy)
    # The data structure of parent/children is (priority, layer1, layer2,
    # layer3, layer4) as a tuple, so the layer number enumerate starting
    # from 1.

    # future person beware -- parent_one[1:2] is for seeding layer, decrement all indices to exclude seed layer -- Bright
    child_one = (None,) + parent_one[1:2] + (LAg,) + parent_one[3:4] + (LAgtop,) + parent_two[5:6]
    child_two = (None,) + parent_one[1:2] + parent_two[2:4] + parent_one[4:5] + (LAlNH,)
    child_three = (None,) + parent_one[1:2] + parent_one[4:6] + parent_two[2:4]
    child_four = (None,) + parent_one[1:2] + parent_two[4:6] + parent_one[2:4]
    
#### Uncomment for 6-layer structure

##    # not worked out yet
##    child_one=parent_one[:3]+parent_two[3:6]+(L5_thickness, )
##    
##    child_two=parent_two[:3]+(L4_thickness, )+parent_one[4:]
##        
##    child_three=child_one[0:1]+child_one[1:2]+child_one[4:]+child_one[3:4]+child_one[2:3]
##    
##    child_four=child_two[0:1]+child_two[1:2]+child_two[4:]+child_two[3:4]+child_two[2:3]
    

    
    for list_ in [child_one, child_two, child_three, child_four]:

        elements.put(recalculate_priority(list_))
        elements.get()    #remove an element with the lowest priority value
        
def change_generation():
    '''
Remove the population dying after each generation and add a new population
to create a new generation containing the survivors. Only 20 population
survive at the end of each generation and the remaining 60 die. 
'''

#Remove the population dying after each generation
    temp_var=0
    while temp_var<numberOfLoser:
        elements.get()
        temp_var+=1
    
    counter=0
    while counter < new_population:
         
        elements.put(create_structure())
        
# start removing onece the population has reached 80
        if counter>=numberOfLoser:
            elements.get()
            
        counter+=1

def calculate_priority(arg):
    '''
This function returns the recalculated priority value by recalculate_priority() function. 
'''
    return recalculate_priority(arg)[0]


def plot_structure(arg, text_only=False):
    '''
Basically, this function plots the graph of reflectance and transmittance of a
given three layer structure versus wavelength. The argument is given as a tuple or a list of 3 elements.
In addition, it provides the values of reflectance, transmittance, R_color and T_color.
'''
    seed = AlN(arg[1])
    first = Ag(arg[2])
    second = AlNH(arg[3])
    third = Ag(arg[4])
    fourth = AlNH(arg[5])
    mystruct = ML([seed, first, second, third, fourth],
                  min_wl=MIN_WL, max_wl=MAX_WL)
    TR(mystruct)

    show()
    


#### Uncomment for 6-layer structure

##    fifth = Ag(arg[4])
##    sixth = DLC80WA(arg[5])
##    
##    mystruct = ML([first, second, third, fourth,
##                   fifth, sixth, AlN(20)],
##                  min_wl=230, max_wl=2476)
    print mystruct

    # Doing unchecked attribute mutation here...future person beware. This
    # only applies because of comment above in correct_TR
    mystruct = correct_TR(mystruct)
    
    index_lowerUV=0
    index_lowervis=0
    index_uppervis=0
    upper_middle=0
    index_upper=0

    if mystruct.wl[len(mystruct.wl)-1]<3000:
        index_upper=len(mystruct.wl)-1
    for index,  i in enumerate(mystruct.wl):
        if index_lowervis ==0 and i>=myConst.lower_limit:
            index_lowervis = index
        elif index_uppervis==0 and i>=myConst.upper_limit:
            index_uppervis=index
        elif myConst.upper_limit!= 700 and upper_middle==0:
            if i>=700:
                upper_middle=index
        elif index_upper==0 and i>=3000:
            index_upper=index
            break
    if myConst.upper_limit==700:
        upper_middle=index_uppervis
        
    photopic_array = np.interp(mystruct.wl[index_lowervis:index_uppervis],P_data[:, 0] , P_data[:, 1])
    V_normalized_data = np.interp(mystruct.wl[index_lowervis:index_uppervis],V_data[:, 0], V_data[:, 3]) / max(V_data[:, 3])
    T_ideal_array = photopic_array/V_normalized_data
    #T_one_array = [1]*len(mystruct.wl[index_lowervis:index_uppervis])
	
    sol_array = np.interp(mystruct.wl[index_lowerUV:index_upper],S_data[:, 0] , S_data[:, 3])
    sol_array_UV = np.interp(mystruct.wl[index_lowerUV:index_lowervis],U_data[:, 0] , U_data[:, 3])
    sol_array_vis = np.interp(mystruct.wl[index_lowervis:index_uppervis],V_data[:, 0] , V_data[:, 3])
    sol_array_IR = np.interp(mystruct.wl[index_uppervis:index_upper],I_data[:, 0] , I_data[:, 3])
    
    Tvis = sum(mystruct.T[index_lowervis:index_uppervis]*photopic_array)/(sum(photopic_array))
   
    #Tpriority = sum(mystruct.T[index_lowervis:index_uppervis]*T_ideal_array)/(sum(T_ideal_array))
    #Tpriority = sum(mystruct.T[index_lowervis:index_uppervis]*sol_array_vis*photopic_array)/(sum(sol_array_vis*photopic_array))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]*sol_array_vis/max(sol_array_vis)-photopic_array)**2))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]-sol_array_vis/max(sol_array_vis))**2))
    #Tpriority = 1/np.sqrt(sum((mystruct.T[index_lowervis:index_uppervis]-T_one_array)**2))
    Tpriority = Tvis
    
    TSER_UV = 1 - sum(mystruct.T[index_lowerUV:index_lowervis]*sol_array_UV)/(sum(sol_array_UV))
    TSER_vis = 1 - sum(mystruct.T[index_lowervis:index_uppervis]*sol_array_vis)/(sum(sol_array_vis))
    TSER_IR = 1 - sum(mystruct.T[index_uppervis:index_upper]*sol_array_IR)/(sum(sol_array_IR))
    # 0.0342 = % of AM1.5 energy in UV
    # 0.426 = % of AM1.5 energy in visible
    # 0.54 = % of AM1.5 energy in IR (700-2500)
    #TSER = TSER_UV * 0.0342 + TSER_vis * 0.426 + TSER_IR * 0.54 - .04
    TSER = 1 - sum(mystruct.T[index_lowerUV:index_upper]*sol_array)/(sum(sol_array)) - 0.04
	
    #TSER_priority = TSER_UV * 0.0342 / (0.0342 + 0.54) + TSER_IR * 0.54 / (0.0342 + 0.54)
    TSER_priority = TSER_IR

    L,a,b = GenerateColor.color_calc(mystruct.wl, mystruct.T, mystruct.R)
    Cpriority = 1/ (abs(a) + abs(b))
	
    #priority = myConst.p1 * Tpriority + myConst.p2 * (TSER_IR - myConst.wanted_TSER)
    priority = myConst.p1 * Tpriority + myConst.p2 * TSER_priority + myConst.p3 * Cpriority 
    #priority = myConst.p1 * Tpriority * Tvis + myConst.p2 * TSER_priority
    #priority = -1 / (Tvis - myConst.p1) - \
    #           abs(myConst.p2 * (TSER - myConst.wanted_TSER))                            

    print 'priorities are: ',myConst.p1, myConst.p2, myConst.p3
    print 'TSER = %s' % float("{0:.4f}".format(TSER))
    print 'TSER in UV = %s' % float("{0:.4f}".format(TSER_UV))
    print 'TSER in visible = %s' % float("{0:.4f}".format(TSER_vis))
    print 'TSER in IR = %s' % float("{0:.4f}".format(TSER_IR))

#    print 'Tpriority = %s' % float("{0:.4f}".format(Tpriority))

    print 'Tvis = %s' % float("{0:.4f}".format(Tvis)) 
    
    print 'T_color:', L,a,b

#    print 'UV lower bound:', index_lowerUV
 #   print 'Visible lower bound:', index_lowervis
  #  print 'Visible upper bound:', index_uppervis
   # print 'IR upper bound:', index_upper
    #print mystruct.wl[:]
    
    if not text_only:
        TR(mystruct, show_solar=True, min_wl=200, max_wl=2500)
        plot.plt.grid("on")
        plot.plt.plot([400,400], [0,1], "k--")
        plot.plt.plot([700,700], [0,1], "k--")
        plot.plt.text(425, 0.70, "visible", fontsize=16)

        show()

    return mystruct
