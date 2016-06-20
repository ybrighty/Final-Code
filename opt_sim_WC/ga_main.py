'''Import one of the modules from the project folder
(we have ga_3layers, ga_5layers, ga_7layers and ga_7layers).

'''
from ga_4layers_bright import *

def optimize_layer(elem, layer_index):
    '''
This function basically add or subtract 0.5nm from a given layer
until it finds optimum thickness. For every iteration the function
to calculate priority is called and the new priority value gets compared
with the previous one inorder to decide either to keep or stop iterating.
The function returns the new set of thickness values when the iteration is complete.
'''

    if (elem[layer_index]+0.25)<=myConst.max_DLC_thickness:
        priority_value=elem[0]
        elem=elem[:layer_index]+(elem[layer_index]+0.25, ) + elem[(layer_index+1):]
        new_priority = calculate_priority(elem)
    else:
        return elem
    if new_priority>priority_value:
        while True:
            if (elem[layer_index]+0.25)<=myConst.max_DLC_thickness:
                priority_value=new_priority
                elem=elem[:layer_index]+(elem[layer_index]+0.25, ) + elem[(layer_index+1):]
                new_priority=calculate_priority(elem)
            else:
                return elem 
                
            if new_priority<priority_value:
                elem=elem[:layer_index]+(elem[layer_index]-0.25, ) + elem[(layer_index+1):]
                elem=(priority_value, )+elem[1:]
                return elem      
    else:
        elem=elem[:layer_index]+(elem[layer_index]-0.25, ) + elem[(layer_index+1):]
        new_priority = calculate_priority(elem)
        while True:
            
            if (elem[layer_index]-0.25)>10:
                priority_value=new_priority
                elem=elem[:layer_index]+(elem[layer_index]-0.25, ) + elem[(layer_index+1):]
                new_priority = calculate_priority(elem)
            else:
                return elem
            if new_priority<priority_value:
                elem=elem[:layer_index]+(elem[layer_index]+0.25, ) + elem[(layer_index+1):]
                elem=(priority_value, )+elem[1:]
                return elem


                
        

def find_maximum(opt_val):
    '''
This function calls the above function to iterate through each layer of a
given structure until the optimum thicknesses in the structure is obtained.
It starts with the argument (opt_val), then attempts to find a better solution by incrementally
changing the thicknesses. If the change produces a better solution, an incremental change is
made to the new solution, repeating until no further improvements can be found.

'''
    while True:
        old_priority=opt_val[0]
        counter=1
        calc_len=len(opt_val)
        while counter<calc_len:
            opt_val=optimize_layer(opt_val, counter)
            counter+=1
        new_priority=opt_val[0]
        if new_priority==old_priority:
            return opt_val
        elif new_priority<old_priority:
            print ('Eror: new_priority is greater than old_priority')

def adjust_format(arg):
    '''
This function basically adjust long decimal place floating point number to one decimal place.
'''
    counter=0
    while counter < len(arg) - 2:
        arg=arg[:counter]+(float("{0:.1f}".format(arg[counter])), )+arg[(counter+1):]
        counter+=1
    return arg
   


numberOfGenerations=40



#Basically the probaility distribution gets created here.

xk = np.arange(80)
total=sum(xk)
pk = [i/(1.*total) for i in xk]
probability = stats.rv_discrete(name='probability', values=(xk, pk))


create_inital_population()

survivalAges=30 #survival age equals to 30 means 30x4=120 offsprings get created in each generation
counter=0
while counter<numberOfGenerations:
    counter_two=0
    
    if elements.qsize()>80:
        print ("Eror:we have a population greater than 80 and which is %s" %elements.qsize())


    while counter_two<survivalAges:
        counter_two+=1
        
#selection of parents for crossover and mutation

        parent_one_=probability.rvs(size=1)
        parent_two_=probability.rvs(size=1)
        
        while parent_one_==parent_two_:
            parent_two_=probability.rvs(size=1)
            
        
        crossover_and_mutation(parent_one_, parent_two_)

    counter+=1
    print ('Generation %s' % counter)
    change_generation()
    

#Remove 70 elements while keeping the top 10 optimized values
temp_var=0
while temp_var<(numberOfLoser+10):
    elements.get()
    temp_var+=1

final_top_ten = elements.queue[:]
final_top_ten.sort()

for result in final_top_ten:
    print (result)
best_struct = plot_structure(final_top_ten[-1])

def output_TR(fname, check=True):
    global best_struct
    #rev_layers = best_struct._layers_list[:]
    #rev_layers.reverse()
    #rev = ML(rev_layers, min_wl=230, max_wl=2476, n0=1.5, ns=1.0)
    #rev.T, rev.R = correct_TR(rev, incidence="bot")
    rev = corrrect_TR(best_struct, incidence="bot")
    if check:
        TR(rev, min_wl=200, max_wl=2500, show_solar=True, legend=True)
        plot.plt.grid('on')
        show()
    os.chdir("..")
    np.savetxt(fname, np.vstack([
        best_struct.wl, best_struct.T, best_struct.R, rev.R]).T,
               fmt=['%d', '%.4f', '%.4f', '%.4f'],
               header="Wavelength (nm) T fR bR")
    os.chdir("opt_sim")

output_TR('result.txt')

##final_top_ten=PriorityQueue()
##index=0
##
##print 'Reached final stage of optimization and the top ten optimized results will be shown below:'
##
##while index<elements.qsize():
##    opt_val=find_maximum(elements.queue[index])
##    final_top_ten.put(opt_val)
##    opt_val =(float("{0:.4f}".format(opt_val[0])), )+adjust_format(opt_val[1:])
##    print opt_val
##    index+=1
##
##
##print 'The final optimized structure:'
##
##Best_result=max(final_top_ten.queue[:])
##Best_result=(float("{0:.4f}".format(Best_result[0])), )+adjust_format(Best_result[1:])
##
##plot_structure(Best_result[1:])

