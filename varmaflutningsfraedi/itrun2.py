
from scipy.optimize import fsolve
import math
T_inf1 = 95
T_inf2 = 18
h_1 = 1200
a = 0.01
A_p = (a/2)**2*math.pi
k = 16
h_2 = 45
d = 0.003
L = 0.01
dx = 0.003
P = math.pi*d
A = math.pi*(d/2)**2
A_b = A_p - A
L_c = L+d/4
m = math.sqrt(P*h_2/(k*A))
M = math.sqrt(h_2*P*k*A)

T_b = fsolve(lambda T_b: k*A_p*(40/49*T_b+855/49-T_b)/dx-(h_2*A_b*(T_b-T_inf2)+(T_b-T_inf2)*M*math.tanh(m*L_c)), [80])[0]

print(T_b)
print(m)
print(M*math.tanh(m*L_c))
print((h_2*A_b*(T_b-T_inf2)+(T_b-T_inf2)*M*math.tanh(m*L_c)))