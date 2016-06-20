from opt_sim.structure import MultiLayer as ML
from opt_sim import *
import random
from Queue import PriorityQueue
from scipy import stats

elements=PriorityQueue() 

populationSize=80

numberOfLoser=60 # The number of population dying after each generation

new_population=100


P_data=np.loadtxt("plot support files\\photopic.txt", skiprows=1)
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
    def __init__(self,lower_limit, upper_limit, T_factor, R_factor, max_DLC_thickness):
        self.lower_limit=lower_limit
        self.upper_limit=upper_limit
        self.T_factor=T_factor
        self.R_factor=R_factor
        self.max_DLC_thickness=max_DLC_thickness


myConst=constants(lower_limit=390, upper_limit=700, T_factor=.55, R_factor=.45, max_DLC_thickness=150.0)



def create_structure():
    '''
  This function creates a structure composed of randomly generated thicknesses.
  It returns a tuple with a list of layer thicknesses and priority value. 
  
 '''
    L1_thickness=random.uniform(20.0, myConst.max_DLC_thickness)
    L1_thickness=float("{0:.1f}".format(L1_thickness))
    layer_one=DLC80WA(L1_thickness)
        
    L2_thickness=random.uniform(3.0, 20.0)
    L2_thickness=float("{0:.1f}".format(L2_thickness))
    layer_two=Ag(L2_thickness)
        
    L3_thickness=random.uniform(20.0, myConst.max_DLC_thickness)
    L3_thickness=float("{0:.1f}".format(L3_thickness))
    layer_three=DLC80WA(L3_thickness)
        
    L4_thickness=random.uniform(3.0, 20.0)
    L4_thickness=float("{0:.1f}".format(L4_thickness))
    layer_four=Ag(L4_thickness)
        
    L5_thickness=random.uniform(20.0, myConst.max_DLC_thickness)
    L5_thickness=float("{0:.1f}".format(L5_thickness))
    layer_five=DLC80WA(L5_thickness)
    
    L6_thickness = random.uniform(3.0, 20.0)
    L6_thickness=float("{0:.1f}".format(L6_thickness))
    layer_six = Ag(L6_thickness)
        
    L7_thickness = random.uniform(20.0, myConst.max_DLC_thickness)
    L7_thickness=float("{0:.1f}".format(L7_thickness))
    layer_seven = DLC80WA(L7_thickness)
        
    mystruct = ML([layer_one, layer_two, layer_three, layer_four, layer_five, layer_six, layer_seven])
    ML.calculate_TR(mystruct)
        
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
    R_array=  np.interp(mystruct.wl[upper_middle:index_upper],S_data[:, 0] , S_data[:, 3])
        
    Transmittance = sum(mystruct.T[index_lower:index_middle]*T_array)/(sum(T_array))
    Reflectance = sum(mystruct.R[upper_middle:index_upper]*R_array)/(sum(R_array))
            
    priority= myConst.R_factor*Reflectance + myConst.T_factor*Transmittance
    
    return (priority, L1_thickness, L2_thickness, L3_thickness, L4_thickness, L5_thickness, L6_thickness, L7_thickness)        

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
    first=DLC80WA(arg[1])
    second = Ag(arg[2])
    third=DLC80WA(arg[3])
    fourth=Ag(arg[4])
    fifth=DLC80WA(arg[5])
    sixth=Ag(arg[6])
    seventh=DLC80WA(arg[7])
    
    mystruct = ML([first, second, third, fourth, fifth, sixth, seventh])
    ML.calculate_TR(mystruct)
            

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
    R_array=  np.interp(mystruct.wl[upper_middle:index_upper],S_data[:, 0] , S_data[:, 3])
        
    Transmittance = sum(mystruct.T[index_lower:index_middle]*T_array)/(sum(T_array))
    Reflectance = sum(mystruct.R[upper_middle:index_upper]*R_array)/(sum(R_array))
    
    priority = myConst.R_factor*Reflectance + myConst.T_factor*Transmittance
    return (priority, arg[1], arg[2], arg[3], arg[4], arg[5], arg[6], arg[7])


#the fitness is not recalculated 
def crossover_and_mutation(parent_one_, parent_two_):
    '''
 This function combines the selected parents to create two children and then,
 interchage the position of some elements to create two more new children.
 Finally, children get added to the prirority queue.

'''

    parent_one=elements.queue[parent_one_]  
    parent_two=elements.queue[parent_two_]

    L6_thickness = random.uniform(3.0, 20.0)
    L6_thickness=float("{0:.1f}".format(L6_thickness))
        
    L7_thickness = random.uniform(20.0, myConst.max_DLC_thickness)
    L7_thickness=float("{0:.1f}".format(L7_thickness))
    
    
    child_one=parent_one[0:4]+parent_two[4:6] +(L6_thickness, L7_thickness)
    
    child_two=parent_two[0:4]+parent_one[4:6] +(L6_thickness, L7_thickness)

        
    child_three=child_one[:1]+child_one[7:]+child_one[4:5]+child_one[3:4]+child_one[2:3]+child_one[5:7]+child_one[1:2]
    
    child_four=child_two[:1]+child_two[1:2]+child_two[4:5]+child_two[3:4]+child_two[2:3]+child_two[7:]+child_two[6:7]+child_two[5:6]
     

    
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

            
def plot_structure(arg):
    '''
Basically, this function plots the graph of reflectance and transmittance of a
given three layer structure versus wavelength. The argument is given as a tuple or a list of 3 elements.
In addition, it provides the values of reflectance, transmittance, R_color and T_color.
'''
    first = DLC80WA(arg[0])
    second = Ag(arg[1])
    third = DLC80WA(arg[2])
    fourth = Ag(arg[3])
    fifth = DLC80WA(arg[4])
    sixth=Ag(arg[5])
    seventh=DLC80WA(arg[6])
    mystruct = ML([first, second, third, fourth, fifth, sixth, seventh])
    print mystruct
    
    mystruct.calculate_color()
    TR(mystruct)
    
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
    R_array=  np.interp(mystruct.wl[upper_middle:index_upper],S_data[:, 0] , S_data[:, 3])
        
    Transmittance = sum(mystruct.T[index_lower:index_middle]*T_array)/(sum(T_array))
    Reflectance = sum(mystruct.R[upper_middle:index_upper]*R_array)/(sum(R_array))

    priority_value = myConst.R_factor*Reflectance + myConst.T_factor*Transmittance
    
    print 'Reflectance=%s' % float("{0:.4f}".format(Reflectance))

    print 'Transmittance =%s' % float("{0:.4f}".format(Transmittance))
    
    print 'T color:', mystruct.T_color
    print 'R color:', mystruct.R_color

    show()   

    
        
        
    
        
        
    
    
