from scipy.optimize import fsolve,least_squares,minimize,Bounds,root_scalar,brentq
import numpy as np

TL_node_buoy = lambda d: 20*np.log10(d) + 0.0006*d + 10
P_t_node_buoy = lambda d: 10**((20-(173.76-TL_node_buoy(d)-70))/10)
t_tx_node_buoy = 0.8
E_tx_node_buoy = lambda d: P_t_node_buoy(d)*t_tx_node_buoy
t_total_node_buoy = lambda d: t_tx_node_buoy + d/1500
E_tx_buoy_buoy = 20.096e-3
t_tx_buoy_buoy = 0.82e-3

def obj(d):
    return E_tx_node_buoy(d)*2+E_tx_buoy_buoy-0.8776 
res = brentq(lambda d:obj(d),1,100000)
print('Energy savings stop at d = ',res,'m')
def obj_time(d):
    return t_total_node_buoy(d)*2+t_tx_buoy_buoy - 4.8
res = brentq(lambda d:obj_time(d),1,100000)
print('Time savings stop at d = ',res,'m')