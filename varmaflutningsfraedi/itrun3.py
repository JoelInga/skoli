from scipy.optimize import fsolve
import numpy as np


q = 100
k = 205
T_0 = 80
T_inf = 12
t = 0.003
L = 0.03
L_c = L + t/2
B  = 0.5
A =  B*t
P = 2*B+2*t
m = lambda h: np.sqrt((h*P*L**2)/(k*A)) 

#reikna h
h = fsolve(lambda h: np.sqrt(h*P*k*A)*(T_0-T_inf)*np.tanh(m(h))-q, [40])[0]
print(h)
#reikna Re_L
k_loft = 0.0245
Nu_L = h*B/k_loft
Pr = 0.71
c_f = lambda Re_L: 0.523*(np.log(0.06*Re_L))**(-2)-1520/Re_L
Re_L = fsolve(lambda Re_L: 2*Nu_L/(c_f(Re_L)*Pr**(1/3))-Re_L,[10**6])[0]
print(Re_L)