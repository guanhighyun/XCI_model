# The code is adapted from https://github.com/elstonlab/mating-dynamics/tree/master by Amy A Pomeroy
import numpy as np
from scipy.integrate import odeint
from deap import base, creator, tools, algorithms
import os
import sys
import pickle
import time as timeski
import math as math
import warnings
warnings.filterwarnings('ignore')


# Conversion matrix
def make_conversion_matrix():
    # want easily savable matrix to hold this info
    # interp boolean, interp range (min,max), power boolean, power number (y)
    len_ind = 15
    arr_IandP = np.zeros((5,len_ind))
    # Set all interp booleans to 1 - everything is going to be interpreted
    arr_IandP[0,:] = 1
    # Set all power booleans to 1 - everything is in the form of powers
    arr_IandP[3,:] = 1
    # Set all power numbers to 10 - everything has a base of 10
    arr_IandP[4,:] = 10

    # Rate constants to be optimized. Each has a range from minimum to maximum.
    # Units see Table 1 in the manuscript.
    # We assume a well-stirred, constant nuclear compartment,
    # so we write units of species in quantities (molecules) instead of concentrations.
    a_act_min = -3 # activator synthesis rate from single X chromosome
    a_act_max =  2
    
    d_act_min = -3 # degradation rate of free activator
    d_act_max = 1

    K_n_min = 0 # quantity of bound SPEN at which activator synthesis rate is half max.
    K_n_max = 4
    
    n_min = 0 # Hill coefficient for SPEN supressing activator synthesis rate
    n_max = 1

    a_x_min = -3 # xist synthesis rate from single X chromosome
    a_x_max = 1
    
    d_x_min = -2.5 # degradation rate of Xist
    d_x_max = -1.5

    m_min = 0  # Hill coefficient for bound SPEN reducing dissociation rate Xist
    m_max = 1
    
    K_S_min = 0 # quantity of SPEN at which Xist dissociation is half max
    K_S_max = 4

    K_a_min = 0 # quantity of activator at which Xist transcription is half max
    K_a_max = 4
    
    k1_min = -4 # rate constant for Xist binding to DNA
    k1_max = 1
    
    k2_min = -3 # dissociation rate for Xist
    k2_max = 1

    k4_min = -1 # dissoication rate for bound SPEN
    k4_max = 1.6
    
    k3_min = -4 # association rate for SPEN
    k3_max = 1

    sT_min = 1.5 # total SPEN quantity
    sT_max = 4.5

    XbsT_min = 2 # total Xist binding sites on one single chromosome
    XbsT_max = 2.3
    
    minimums = [a_act_min, d_act_min, K_n_min, n_min, a_x_min, d_x_min, m_min, K_S_min, K_a_min, k1_min, k2_min, k4_min, k3_min, XbsT_min, sT_min]
    maximums = [a_act_max, d_act_max, K_n_max, n_max, a_x_max, d_x_max, m_max, K_S_max, K_a_max, k1_max, k2_max, k4_max, k3_max, XbsT_max, sT_max]

    for i in range(len(minimums)):
        arr_IandP[1,i] = minimums[i] 
        arr_IandP[2,i] = maximums[i] 
    return arr_IandP


def convert_individual(ea_individual, conversion_matrix):
    # use conversion matrix to convert interp and exponentiate individual:
    # conversion matrix has 5 rows of arrs of size len(individual):
    # conversion_matrix[0] = interp_bool
    # conversion_matrix[1] = interp_range_min
    # conversion_matrix[2] = interp_range_max
    # conversion_matrix[3] = power_bool
    # conversion_matrix[4] = base_val
    
    # copy and get len of individual
    arr_params_conv = np.zeros(15)
    len_ind = len(ea_individual)
    
    # Interp:
    for idx in np.nonzero(conversion_matrix[0])[0]:
        ea_val = ea_individual[idx]
        r_min = conversion_matrix[1][idx]
        r_max = conversion_matrix[2][idx]
        arr_params_conv[idx] = np.interp(ea_val, (0,1), (r_min, r_max))
    
    # Exponentiate:
    for idx in np.nonzero(conversion_matrix[3])[0]:
        ea_val = arr_params_conv[idx]
        base_val = conversion_matrix[4][idx]
        arr_params_conv[idx] = np.power(base_val, ea_val)
    
    return arr_params_conv

# Differential equations for the XY case (a single X chromosome) in males
def DE1(y,t,arr_parameters_IP):
    # Species: activator (act), free Xist (x1f), bound Xist (x1b), and bound SPEN (s1b)
    act, x1f, x1b, s1b = y
    # Model parameters. Please refer to the notations in the function "make_conversion_matrix"
    a_act, d_act, K_n, n, a_x, d_x, m, K_S, K_a, k1, k2, k4, k3, XbsT, sT = arr_parameters_IP
    # Number of SPEN that binds to one Xist. 
    N_S = round(sT/XbsT) # We let each chromosome be able to recruit and bind to all SPEN.

    dy = [
        a_act/(1 +(s1b/K_n)**n) - d_act*act,
        a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)**m),
        k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)**m),
        k3*(sT - s1b)*(N_S*x1b - s1b)  - k4*s1b
        ]
    return dy


# Differential equations for the XX case in females
def DE2(y,t,arr_parameters_IP):
    act, x1f, x1b, s1b, x2f, x2b, s2b = y
    a_act, d_act, K_n, n, a_x, d_x, m, K_S, K_a, k1, k2, k4, k3, XbsT, sT = arr_parameters_IP
    N_S = round(sT/XbsT)
    dy = [
        a_act/(1 +(s1b/K_n)**n) + a_act/(1 +(s2b/K_n)**n) - d_act*act,
        a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)**m),
        k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)**m),
        k3*(sT - s1b - s2b)*(N_S*x1b - s1b)  - k4*s1b,
        a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)**m),
        k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)**m),
        k3*(sT - s1b - s2b)*(N_S*x2b - s2b)  - k4*s2b,
        ]
    return dy

# Differential equations for the XXX case to test n-1 rule.
def DE3(y,t,arr_parameters_IP):
    act, x1f, x1b, s1b, x2f, x2b, s2b, x3f, x3b, s3b = y
    a_act, d_act, K_n, n, a_x, d_x, m, K_S, K_a, k1, k2, k4, k3, XbsT, sT = arr_parameters_IP
    N_S = round(sT/XbsT)
    dy = [
        a_act/(1 +(s1b/K_n)**n) + a_act/(1 +(s2b/K_n)**n) + a_act/(1 +(s3b/K_n)**n) - d_act*act,
        a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)**m),
        k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)**m),
        k3*(sT - s1b - s2b - s3b)*(N_S*x1b - s1b)  - k4*s1b,
        a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)**m),
        k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)**m),
        k3*(sT - s1b - s2b - s3b)*(N_S*x2b - s2b)  - k4*s2b,
        a_x*act/(K_a + act) - d_x*x3f - k1*(XbsT - x3b)*x3f + k2*x3b/(1+(s3b/K_S)**m),
        k1*(XbsT - x3b)*x3f - k2*x3b/(1+(s3b/K_S)**m),
        k3*(sT - s1b - s2b - s3b)*(N_S*x3b - s3b)  - k4*s3b
        ]
    return dy

# Calculate fitness function for XXX case
def results_of_3X(odes, SS_penalty, time_step, XbsT, sT):
    # extract solutions at 6000 min
    sol_6000s = odes[int(6000/time_step)-1,:]
    x1b_6000s = sol_6000s[2] # bound Xist on X1 at 6000 min
    x2b_6000s = sol_6000s[5] # bound Xist on X2 at 6000 min
    x3b_6000s = sol_6000s[8] # bound Xist on X3 at 6000 min

    s1b_6000s = sol_6000s[3] # bound SPEN on X1 at 6000 min
    s2b_6000s = sol_6000s[6] # bound SPEN on X2 at 6000 min
    s3b_6000s = sol_6000s[9] # bound SPEN on X3 at 6000 min

    # extract solution at the final time point
    SS = odes[-1,:]
    
    x1b_SS = SS[2] # bound Xist on X1 at the final time point
    x2b_SS = SS[5] # bound Xist on X2 at the final time point
    x3b_SS = SS[8] # bound Xist on X3 at the final time point

    s1b_SS = SS[3] # bound SPEN on X1 at the final time point
    s2b_SS = SS[6] # bound SPEN on X2 at the final time point
    s3b_SS = SS[9] # bound SPEN on X3 at the final time point

    # fitness function
    fitness = ((np.abs(x1b_SS-x1b_6000s) + np.abs(x2b_SS-x2b_6000s) + np.abs(x3b_SS-x3b_6000s) 
                + np.abs(s1b_SS-s1b_6000s) + np.abs(s2b_SS-s2b_6000s) + np.abs(s3b_SS-s3b_6000s))*SS_penalty 
               + np.abs(XbsT/(1+x2b_SS-x1b_SS)) + 0.1*np.abs(s1b_SS + s2b_SS + s3b_SS)) # penalize on bound SPEN because
                                                                                        # we do not want SPEN to be depleted
                                                                                        # in the activator inhibition model
    return fitness

# Calculate fitness function for XX case
def results_of_2X(odes, SS_penalty, time_step, XbsT, sT):
    # Solution at 6000 min
    sol_6000s = odes[int(6000/time_step)-1,:]
    x1b_6000s = sol_6000s[2] # bound Xist on X1 at 6000 min 
    x2b_6000s = sol_6000s[5] # bound Xist on X2 at 6000 min 

    s1b_6000s = sol_6000s[3] # bound SPEN on X1 at 6000 min 
    s2b_6000s = sol_6000s[6] # bound SPEN on X2 at 6000 min 

    # Solution at the final time point
    SS = odes[-1,:]
    
    x1b_SS = SS[2] # bound Xist on X1 at the final time point
    x2b_SS = SS[5] # bound Xist on X2 at the final time point

    s1b_SS = SS[3] # bound SPEN on X1 at the final time point
    s2b_SS = SS[6] # bound SPEN on X2 at the final time point

    fitness = ((np.abs(x1b_SS-x1b_6000s) + np.abs(x2b_SS-x2b_6000s) + np.abs(s1b_SS-s1b_6000s) + np.abs(s2b_SS-s2b_6000s))*SS_penalty
               + np.abs(XbsT/(1+x2b_SS-x1b_SS)) + 0.1*np.abs(s1b_SS + s2b_SS))
    return fitness

# Calculate fitness function for XY case
def results_of_1X(odes, SS_penalty, time_step, XbsT):
    # Solution at 6000 min
    sol_6000s = odes[int(6000/time_step)-1,:]
    x1b_6000s = sol_6000s[2] # bound Xist at 6000 min
    s1b_6000s = sol_6000s[3] # bound SPEN at 6000 min

    # Solution at the final time point
    SS = odes[-1,:]
    x1b_SS = SS[2] # bound Xist at the final time point
    s1b_SS = SS[3] # bound SPEN at the final time point

    fitness = (np.abs(x1b_SS-x1b_6000s) + np.abs(s1b_SS-s1b_6000s))*SS_penalty + np.abs(x1b_SS)
    return fitness

# Simulate all XY, XX and XXX cases and sum up all the fitness functions
def scorefxn1(arr_parameters, final_time, time_step, SS_penalty):
    # Load parameters sampled from uniform distributions
    arr_conversion_matrix = make_conversion_matrix()
    arr_parameters_IP = convert_individual(arr_parameters, arr_conversion_matrix)
    a_act, d_act, K_n, n, a_x, d_x, m, K_S, K_a, k1, k2, k4, k3, XbsT, sT = arr_parameters_IP  

    # set simulation time
    n_time_point = int(final_time/time_step)
    time = np.linspace(0,final_time,n_time_point)

    ################ XY case #################

    # Initial condition. 
    IC = [0,  #act;
          0,  #x1f;
          XbsT,  #x1b;
          0,  #s1b;
         ]

    # Simulate ODE
    odes = odeint(DE1, IC, time, args=(arr_parameters_IP,) , rtol = 1e-8, atol=1e-8)

    # Calculate fitness function value
    fitness_1X = results_of_1X(odes, SS_penalty, time_step, XbsT)

    ################ XX case #################
    
    # First initial condition (very small amount of initial bound Xist). 
    IC = [0,  #act;
          0,  #x1f;
          1,  #x1b;
          0,  #s1b;
          0,  #x2f;
          2, #x2b
          0, #s2b
         ]
    
    odes = odeint(DE2, IC, time, args=(arr_parameters_IP,) , rtol = 1e-8, atol=1e-8)
    fitness_2X_1 = results_of_2X(odes, SS_penalty, time_step, XbsT, sT)

    # Second initial condition (very large amount of initial bound Xist). 
    IC = [0,  #act;
          0,  #x1f;
          XbsT-1,  #x1b;
          0,  #s1b;
          0,  #x2f;
          XbsT, #x2b
          0, #s2b
         ]
    
    odes = odeint(DE2, IC, time, args=(arr_parameters_IP,) , rtol = 1e-8, atol=1e-8)
    fitness_2X_2 = results_of_2X(odes, SS_penalty, time_step, XbsT, sT)

    ################ XXX case #################
    
    # First initial condition. 
    IC = [0,  #act;
          0,  #x1f;
          1,  #x1b;
          0,  #s1b;
          0,  #x2f;
          2, #x2b
          0, #s2b
          0, #x3f
          3, #x3b
          0 #s3b
         ]
    
    odes = odeint(DE3, IC, time, args=(arr_parameters_IP,) , rtol = 1e-8, atol=1e-8)
    fitness_3X_1 = results_of_3X(odes, SS_penalty, time_step, XbsT, sT)

    # Second initial condition.
    IC = [0,  #act;
          0,  #x1f;
          1,  #x1b;
          0,  #s1b;
          0,  #x2f;
          2, #x2b
          0, #s2b
          0, #x3f
          XbsT, #x3b
          0 #s3b
         ]
    
    odes = odeint(DE3, IC, time, args=(arr_parameters_IP,) , rtol = 1e-8, atol=1e-8)
    fitness_3X_2 = results_of_3X(odes, SS_penalty, time_step, XbsT, sT)

    # Third initial condition. 
    IC = [0,  #act;
          0,  #x1f;
          XbsT-2,  #x1b;
          0,  #s1b;
          0,  #x2f;
          XbsT-1, #x2b
          0, #s2b
          0, #x3f
          XbsT, #x3b
          0 #s3b
         ]
    
    odes = odeint(DE3, IC, time, args=(arr_parameters_IP,) , rtol = 1e-8, atol=1e-8)
    fitness_3X_3 = results_of_3X(odes, SS_penalty, time_step, XbsT, sT)

    # Sum up all fitness function
    fitness = fitness_1X + fitness_2X_1 + fitness_2X_2 + fitness_3X_1 + fitness_3X_2 + fitness_3X_3

    return fitness

# just a helper function that pulls all of scorefxn1 dependencies together
# note the (), <--using single optimization in DEAP for now
def scorefxn_helper(individual):
    SS_penalty_weight = 100000 # penalty for not reaching steady states within the final time point
    time_step = 0.01 # time step
    final_time = 7000 # final time point of simulations
    return scorefxn1(individual, final_time, time_step, SS_penalty_weight),

# file name of output file
def get_filename(i):
    dir_to_use = os.getcwd()
    filename_base = dir_to_use + '/' + str(i)
    return filename_base + '.pickled'

# EA hyperparameters
number_of_runs = 1
number_of_generations = 100
number_of_individuals = 500
mutation_rate = 0.1
crossover_rate = 0.5
arr_conversion_matrix = make_conversion_matrix()
filename_number = sys.argv[1]

for i in range(number_of_runs):
    ###################################################################
    #EVOLUTIONARY ALGORITHM
    ###################################################################
    #TYPE
    #Create minimizing fitness class w/ single objective:
    creator.create('FitnessMin', base.Fitness, weights=(-1.0,))
    #Create individual class:
    creator.create('Individual', list, fitness=creator.FitnessMin)

    #TOOLBOX
    toolbox = base.Toolbox()
    #Register function to create a number in the interval [1-100?]:
    #toolbox.register('init_params', )
    #Register function to use initRepeat to fill individual w/ n calls to rand_num:
    toolbox.register('individual', tools.initRepeat, creator.Individual, 
                     np.random.random, n=48)
    #Register function to use initRepeat to fill population with individuals:
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)

    #GENETIC OPERATORS:
    # Register evaluate fxn = evaluation function, individual to evaluate given later
    toolbox.register('evaluate', scorefxn_helper)
    # Register mate fxn = two points crossover function 
    toolbox.register('mate', tools.cxTwoPoint)
    # Register mutate by swapping two points of the individual:
    toolbox.register('mutate', tools.mutPolynomialBounded, 
                     eta=0.1, low=0.0, up=1.0, indpb=0.2)
    # Register select = size of tournament set to 3
    toolbox.register('select', tools.selTournament, tournsize=3)

    #EVOLUTION!
    pop = toolbox.population(n=number_of_individuals)
    hof = tools.HallOfFame(1)

    stats = tools.Statistics(key = lambda ind: [ind.fitness.values, ind])
    stats.register('all', np.copy)

    # using built in eaSimple algo
    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=crossover_rate, 
                                       mutpb=mutation_rate, 
                                       ngen=number_of_generations, 
                                       stats=stats, halloffame=hof, 
                                       verbose=False)

    #MAKE LISTS
    # Find best scores and individuals in population 
    arr_best_score = []
    arr_best_ind = []
    for a in range(len(logbook)):
        scores = []
        for b in range(len(logbook[a]['all'])):
            scores.append(logbook[a]['all'][b][0][0])
        #print(a, np.nanmin(scores), np.nanargmin(scores))
        arr_best_score.append(np.nanmin(scores))
        #logbook is of type 'deap.creator.Individual' and must be loaded later
        #don't want to have to load it to view data everytime, thus numpy
        ind_np = np.asarray(logbook[a]['all'][np.nanargmin(scores)][1])
        ind_np_conv = convert_individual(ind_np, arr_conversion_matrix)
        arr_best_ind.append(ind_np_conv)
        #arr_best_ind.append(np.asarray(logbook[a]['all'][np.nanargmin(scores)][1]))
        arr_to_pickle = [arr_best_score, arr_best_ind]


    filename = get_filename(filename_number)
    pickle.dump(arr_to_pickle, open(filename,'wb'))
