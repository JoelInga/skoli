from py2gmsh import (Mesh, Entity, Field)
import axifem
import os
from scipy.optimize import minimize
import numpy as np
h_innan = 1200
T_inf_innan = 95
T_inf_utan = 18
h_utan = 45
k = 16
t_veggur = 0.003
d = 0.003
a = 0.01
L = 0.01
d1 = 0.0005
d2 = 0.0005
V_sk = np.pi*(a/2)**2*t_veggur+np.pi*(d/2)**2*(L)

filename = 'my_mesh'
os.system('gmsh {}.geo -2 -o {}.msh'.format(filename,filename))
x, y, tri, T, V, q = axifem.axiHeatCond('my_mesh.msh', \
            {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})

from matplotlib.pyplot import *

print('Rúmmál: {}'.format(V['ribba']))
print('Varmaflæði: {:g}'.format(q['ytri'][1]))
print('Hámarkshitastig: {:g}'.format(max(T)))
print('Lágmarkshitastig: {:g}'.format(min(T)))
figure(figsize=(16,3))
tricontourf(x,y,tri,T,20)
tricontourf(x,-y,tri,T,20)
colorbar()
axis('equal')
show()


3.063052837250048e-07