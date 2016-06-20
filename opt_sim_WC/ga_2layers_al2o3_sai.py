from opt_sim.structure import MultiLayer as ML
from opt_sim import *
import random
from Queue import PriorityQueue
from scipy import stats

elements=PriorityQueue() 

populationSize=80

numberOfLoser=60 # The number of population dying after each generation

new_population=100


P_data=np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
S_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)



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
    def __init__(self,lower_limit, upper_limit, p1, p2, wanted_TSER, max_DLC_thickness):
        self.lower_limit=lower_limit
        self.upper_limit=upper_limit
        self.p1=p1
        self.p2=p2
        self.wanted_TSER=wanted_TSER
        self.max_DLC_thickness=max_DLC_thickness


myConst=constants(lower_limit=400, upper_limit=700,
                  p1=1.2, p2=5, wanted_TSER=.50, max_DLC_thickness=100.0)


## FOR correct_TR()
glass_ext = np.loadtxt("glass_extinct_assume_n1.5.txt")
wl_A, A = glass_ext[2:-25,0], glass_ext[2:-25,1]    #210nm to 2475nm Ag/DLC
##wl_A, A = glass_ext[6:-25,0], glass_ext[6:-25,1]    #230nm to 2475nm Ag/DLC/AlN
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

    return T, R


def create_structure():
    '''
  This function creates a structure composed of randomly generated thicknesses.
  It returns a tuple with a list of layer thicknesses and priority value. 
  
 '''
       
    L1_thickness=random.uniform(5., 30.)
    L1_thickness=float("{0:.1f}".format(L1_thickness))
    layer_one=Ag(L1_thickness)
        
    L2_thickness=random.uniform(10.0, 100.0)
    L2_thickness=float("{0:.1f}".format(L2_thickness))
    layer_two=Al2O3(L2_thickness)
        
    mystruct = ML([layer_one, layer_two],
                  min_wl=MIN_WL, max_wl=MAX_WL)
    ML.calculate_TR(mystruct)

    # Doing unchecked attribute mutation here...future person beware. This
    # only applies because of comment above in correct_TR
    mystruct.T, mystruct.R = correct_TR(mystruct)

    index_lower=0
    index_middle=0
    upper_middle=0
    index_upper=0

    if mystruct.wl[len(mystruct.wl)-1]<3000:
        index_upper=len(mystruct.wl)-1
    for index,  i in enumerate(mystruct.wl):
        if index_lower ==0 and i>=myConst.lower_limit:
            index_lower = index
        elif index_middle==0 and i>=myConst.upper_limit:
            index_middle=index
        elif myConst.upper_limit!= 700 and upper_middle==0:
            if i>=700:
                upper_middle=index
        elif index_upper==0 and i>=3000:
            index_upper=index
            break
    if myConst.upper_limit==700:
        upper_middle=index_middle
    
    T_array = np.interp(mystruct.wl[index_lower:index_middle],P_data[:, 0] , P_data[:, 1])
    sol_array=  np.interp(mystruct.wl[index_lower:index_upper],S_data[:, 0] , S_data[:, 3])
        
    Tvis = sum(mystruct.T[index_lower:index_middle]*T_array)/(sum(T_array))
    TSER = 1 - sum(mystruct.T[index_lower:index_upper]*sol_array)/(sum(sol_array)) - 0.04
            
    priority = -1 / (Tvis - myConst.p1) - \
               abs(myConst.p2 * (TSER - myConst.wanted_TSER))

    return (priority, L1_thickness, L2_thickness, Tvis, TSER)
        
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
    first=Ag(arg[1])
    second=Al2O3(arg[2])
    mystruct = ML([first, second],
                  min_wl=MIN_WL, max_wl=MAX_WL)
    ML.calculate_TR(mystruct)

    # Doing unchecked attribute mutation here...future person beware. This
    # only applies because of comment above in correct_TR
    mystruct.T, mystruct.R = correct_TR(mystruct)
            

    index_lower=0
    index_middle=0
    upper_middle=0
    index_upper=0

    if mystruct.wl[len(mystruct.wl)-1]<3000:
        index_upper=len(mystruct.wl)-1
    for index,  i in enumerate(mystruct.wl):
        if index_lower ==0 and i>=myConst.lower_limit:
            index_lower = index
        elif index_middle==0 and i>=myConst.upper_limit:
            index_middle=index
        elif myConst.upper_limit!= 700 and upper_middle==0:
            if i>=700:
                upper_middle=index
        elif index_upper==0 and i>=3000:
            index_upper=index
            break
    if myConst.upper_limit==700:
        upper_middle=index_middle
        
    T_array = np.interp(mystruct.wl[index_lower:index_middle],P_data[:, 0] , P_data[:, 1])
    sol_array=  np.interp(mystruct.wl[index_lower:index_upper],S_data[:, 0] , S_data[:, 3])
        
    Tvis = sum(mystruct.T[index_lower:index_middle]*T_array)/(sum(T_array))
    TSER = 1 - sum(mystruct.T[index_lower:index_upper]*sol_array)/(sum(sol_array)) - 0.04
            
    priority = -1 / (Tvis - myConst.p1) - \
               abs(myConst.p2 * (TSER - myConst.wanted_TSER))

    return (priority, arg[1], arg[2], Tvis, TSER)

#the fitness is not recalculated 
def crossover_and_mutation(parent_one_, parent_two_):
    '''
 This function combines the selected parents to create two children and then,
 interchage the position of some elements to create two more new children.
 Finally, children get added to the prirority queue.

'''
    parent_one=elements.queue[parent_one_]  
    parent_two=elements.queue[parent_two_]
    
    Ag_thickness=random.uniform(5., 30.)
    Ag_thickness=float("{0:.1f}".format(Ag_thickness))
    
    Al2O3_thickness=random.uniform(10.0, 100.0)
    Al2O3_thickness=float("{0:.1f}".format(Al2O3_thickness))
    
    child_one=parent_one[0:2]+(Ag_thickness, )
    
    child_two=parent_two[0:2]+(Ag_thickness, )

    child_three=child_one[0:1]+(Al2O3_thickness, )+child_one[2:]
    
    child_four=child_two[0:1]+(Al2O3_thickness, )+child_two[2:]
        
    
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
    first = Ag(arg[1])
    second = Al2O3(arg[2])
    mystruct = ML([first, second],
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
    
##    mystruct.calculate_color()

    # Doing unchecked attribute mutation here...future person beware. This
    # only applies because of comment above in correct_TR
    mystruct.T, mystruct.R = correct_TR(mystruct)
    
    index_lower=0
    index_middle=0
    upper_middle=0
    index_upper=0

    if mystruct.wl[len(mystruct.wl)-1]<3000:
        index_upper=len(mystruct.wl)-1
    for index,  i in enumerate(mystruct.wl):
        if index_lower ==0 and i>=myConst.lower_limit:
            index_lower = index
        elif index_middle==0 and i>=myConst.upper_limit:
            index_middle=index
        elif myConst.upper_limit!= 700 and upper_middle==0:
            if i>=700:
                upper_middle=index
        elif index_upper==0 and i>=3000:
            index_upper=index
            break
    if myConst.upper_limit==700:
        upper_middle=index_middle
        
    T_array = np.interp(mystruct.wl[index_lower:index_middle],P_data[:, 0] , P_data[:, 1])
    sol_array=  np.interp(mystruct.wl[index_lower:index_upper],S_data[:, 0] , S_data[:, 3])
        
    Tvis = sum(mystruct.T[index_lower:index_middle]*T_array)/(sum(T_array))
    TSER = 1 - sum(mystruct.T[index_lower:index_upper]*sol_array)/(sum(sol_array)) - 0.04
            
    priority = -1 / (Tvis - myConst.p1) - \
               abs(myConst.p2 * (TSER - myConst.wanted_TSER))
    
    print 'TSER = %s' % float("{0:.4f}".format(TSER))

    print 'Tvis = %s' % float("{0:.4f}".format(Tvis))
    
##    print 'T_color:', mystruct.T_color
##    print 'R_color:', mystruct.R_color

    if not text_only:
        TR(mystruct, show_solar=True, min_wl=200, max_wl=2500)
        plot.plt.grid("on")
        plot.plt.plot([400,400], [0,1], "k--")
        plot.plt.plot([700,700], [0,1], "k--")
        plot.plt.text(425, 0.70, "visible", fontsize=16)

        show()

    return mystruct
