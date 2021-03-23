from scipy.optimize import fsolve,least_squares,minimize,Bounds,root_scalar,brentq
import numpy as np
from CoolProp.CoolProp import PropsSI


T_w = 273.15+2
T_vatn = 273.15+80
d_i = 0.15
d_y = 0.19
L = 0.5
k_st = 0.045

R = np.log(d_y/d_i)/(2*k_st*np.pi*L)

q = (T_vatn-T_w)/R
print('varmastreymi ', q)

nu = lambda T_f:PropsSI('V','T',T_f,'P',1e5,'Air')/PropsSI('D','T',T_f,'P',1e5,'Air')
Pr = lambda T_f:PropsSI('Prandtl','T',T_f,'P',1e5,'Air')
k = lambda T_f:PropsSI('L','T',T_f,'P',1e5,'Air')

h = lambda T_inf: q/((T_w-T_inf)*np.pi*d_y*L)
Gr = lambda T_f,T_inf: 9.8*1/(T_f)*(T_w-T_inf)*d_y**3/(nu(T_f)**2)
Nu_d = lambda Gr,T_f: (0.60 + 0.387*(Gr*Pr(T_f)/(1+(0.559/Pr(T_f))**(9/16))**(16/9))**(1/6))**2


def obj(x):
    return h(x) -  Nu_d( Gr((x+T_w)/2,x), (x+T_w)/2 )*k((x+T_w)/2)/d_y
    

res = brentq(lambda x:obj(x),150,275)
print(res)
print('T_inf = ',res-273.15)

q = h(res)*np.pi*d_y*L*(T_w-res)
print(q)
print('h',h(res))