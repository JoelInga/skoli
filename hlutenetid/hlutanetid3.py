import matplotlib.pyplot as plt
import numpy as np
from numba import njit,jit
from scipy.optimize import brentq
import time

@njit
def number_of_isolated_sensors(lamb, number_of_simulations):
    nr = 0.0 #initialize nr
    for i in range(0,number_of_simulations):
        #lamb = -np.log(0.01)/(np.pi*3.8266853543893022**2) # the rate
        pi = np.pi # pi = 3.14...
        r = np.sqrt(10000.0/np.pi) # the radius of the circle (for the region)
        mean = lamb * pi * r ** 2 # the mean of the Poisson random variable n
        n = np.random.poisson(mean) # the Poisson random variable (i.e., the number of points inside the region)
       
        u_1 = np.random.uniform(0.0, 1.0, n) # generate n uniform (and randomly) distributed points 
        radii = np.zeros(n) # the radial coordinate of the points
        for i in range(n):
            radii[i] = r * (np.sqrt(u_1[i]))
        u_2 = np.random.uniform(0.0, 1.0, n) # generate another n uniformly distributed points 
        angle = np.zeros(n) # the angular coordinate of the points
        for i in range(n):
            angle[i] = 2 * pi * u_2[i]

        x = np.zeros(n) #convert the polar coordinates to cartesian
        y = np.zeros(n)
        for i in range(n):
            x[i] = radii[i] * np.cos(angle[i])
            y[i] = radii[i] * np.sin(angle[i])
            
        shortest_distance = np.zeros(n) #fill an array with zeros (fixed size)

        for i in range(n):
            distance = np.zeros(n)#fill an array with zeros (fixed size)
            for j in range(n):
                if i !=j:
                    distance[j] = np.sqrt((x[j]-x[i])**2.0+(y[j]-y[i])**2.0) #calculate distance between nodes
                else:
                    distance[j] = np.inf #infinite distance for node-same node
                    
            shortest_distance[i] = (np.min(distance)) #find the minimum distance to the next node
        tmp = 0.0
        for i in range(n):
            if shortest_distance[i] > 3.8266853543893022: #if the minimum distance exceeds the coverage radius, it is isolated
                tmp +=1.0
        nr += tmp/n #percentage of isolated sensors for the simulation
        
    nr/=number_of_simulations #average number of isolated sensors over all simulations
                
    return nr*100.0 #percentage of isolated sensors

timeit=False # for performance evaluation
if timeit:
    start = time.time()
    nr = number_of_isolated_sensors(0.1, 1000)
    end = time.time()
    print(end-start)
    
lamb = 0.107
nr_sim = 100

def obj(lamb): #objective function
    return number_of_isolated_sensors(lamb,100)-1.0

def solve_for_lamb():
    lamb = brentq(lambda lamb:obj(lamb),0.01,1.0) #Solver
    return lamb

solve = False #to calculate 
if solve:
    lamb = solve_for_lamb()
    print('ideal density: ', lamb)

print('The percentage of isolated sensors for lamba = {} is {}% after {} simulations'.format(lamb,number_of_isolated_sensors(lamb, nr_sim),nr_sim))

            
plot = True
if plot:
    pi = np.pi # pi = 3.14...
    r = np.sqrt(10000.0/np.pi) # the radius of the circle (for the region)
    mean = lamb * pi * r ** 2 # the mean of the Poisson random variable n
    n = np.random.poisson(mean) # the Poisson random variable (i.e., the number of points inside the region)
    
    u_1 = np.random.uniform(0.0, 1.0, n) # generate n uniform (and randomly) distributed points 
    radii = np.zeros(n) # the radial coordinate of the points
    for i in range(n):
        radii[i] = r * (np.sqrt(u_1[i]))
    u_2 = np.random.uniform(0.0, 1.0, n) # generate another n uniformly distributed points 
    angle = np.zeros(n) # the angular coordinate of the points
    for i in range(n):
        angle[i] = 2 * pi * u_2[i]

    x = np.zeros(n) #convert the polar coordinates to cartesian
    y = np.zeros(n)
    for i in range(n):
        x[i] = radii[i] * np.cos(angle[i])
        y[i] = radii[i] * np.sin(angle[i])
        
        shortest_distance = np.zeros(n) #fill an array with zeros (fixed size)

    for i in range(n):
        distance = np.zeros(n)#fill an array with zeros (fixed size)
        for j in range(n):
            if i !=j:
                distance[j] = np.sqrt((x[j]-x[i])**2.0+(y[j]-y[i])**2.0) #calculate distance between nodes
            else:
                distance[j] = np.inf #infinite distance for node-same node
                
        shortest_distance[i] = (np.min(distance)) #find the minimum distance to the next node
    tmp = 0.0
    fig = plt.gcf()
    ax = fig.gca()
    plt.plot([], [], ' ', label="Lambda = {}".format(lamb))
    plt.plot([], [], ' ', label="Nr of nodes = {}".format(n))
    plt.xlim(-r-5, r+5)
    plt.ylim(-r-5, r+5)
    circ = plt.Circle((0, 0), radius=r, color='g', linewidth=2, fill=False)
    plt.plot(x, y, 'bo')
    x_iso = []
    y_iso = []
    for i in range(n):
        if shortest_distance[i] > 3.8266853543893022: #if the minimum distance exceeds the coverage radius, it is isolated
            x_iso.append(x[i])
            y_iso.append(y[i])#find the isolated sensors
    plt.plot(x_iso, y_iso, 'ro', label = 'isolated nodes')     #mark the isolated sensors
    plt.legend()
    ax.add_artist(circ)
    plt.show()
    
    
    
