from scipy.optimize import fsolve
import numpy as np
n = 0.015
b = 8
vdot = 8.5
s0 = 0.01
yn = fsolve(lambda yn: (n*vdot*(b+2*yn)**(2/3)/(b**(5/3)*np.sqrt(s0)))**(3/5)-yn, [0.9])[0]
print(yn)