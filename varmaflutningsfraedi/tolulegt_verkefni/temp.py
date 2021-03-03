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
d1 = 0.0002
d2 = 0.0005
V_sk = 1.8849555921538741e-06
tolerance = 6e-013 #tolerance for volume


# create Mesh class instance
my_mesh = Mesh()
x = [0.003,0.01,0.01,0.001, 0.0001,0.001]
d = x[0]
a = x[1]
L = x[2]
s1 = x[3]
s2 = x[4]
s3 = x[5]

# create points
p1 = Entity.Point([0., 0., 0.,d1]) #fyrsti punktur neðri vinstri

p2 = Entity.Point([0.,a, 0.,d1])#2. punktur efri vinstri

p3 = Entity.Point([t_veggur, a, 0.,d1])#3. punktur efri hægri


p4 = Entity.Point([t_veggur, a-s1, 0.,d1])#fyrsti ribbup ef



p5 = Entity.Point([t_veggur+L/2,a-s2,0.,d1])# 2 ribbup ef 


p6 = Entity.Point([t_veggur+L,a-s3,0.,d1])#síðast ribbup ef


p7 = Entity.Point([t_veggur+L,s3,0.,d1])#síðast ribbup ne

p8 = Entity.Point([t_veggur+L/2,s2,0,d1])#2. ribbup ne


p9 = Entity.Point([t_veggur, s1, 0.,d1])#1. ribbup ef


p10 = Entity.Point([t_veggur,0.,0.,d1])#síðasti punktur neðri hægri

my_mesh.addEntities([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10])
# create curves
l1 = Entity.Curve([p1, p2]) #innri bein lína upp
l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
l3 = Entity.Curve([p3, p4]) # ytri bein lína upp
l4 = Entity.Curve([p4, p5]) #ribba 1. 2. p e

l5 = Entity.Curve([p5, p6]) #ribba 2. 3. p e
l6 = Entity.Curve([p6, p7]) #ribba endi
l7 = Entity.Curve([p7, p8]) #ribba 3. 2. n
l8 = Entity.Curve([p8, p9]) #ribba 2. 1. n
l9 = Entity.Curve([p9,p10]) #neðri bein lína upp 
l10 = Entity.Curve([p10,p1]) #einangrun neðri

my_mesh.addEntities([l1, l2, l3, l4, l5, l6, l7, l8,l9,l10])


ll1 = Entity.CurveLoop([l1, l2, l3, l4, l5, l6, l7, l8,l9,l10], mesh=my_mesh)



s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




g1 = Entity.PhysicalGroup(name='innri')
g2 = Entity.PhysicalGroup(name='ytri')
g3 = Entity.PhysicalGroup(name='ribba')
g4 = Entity.PhysicalGroup(name='einangrun')
my_mesh.addEntities([g1, g2, g3, g4])
g1.addEntities([l1])
g2.addEntities([l3,l4,l5,l6,l7,l8,l9])
g4.addEntities([l2,l10])
g3.addEntities([s1])


#my_mesh.Options.Mesh.CharacteristicLengthMax = 0.1

# adding Coherence option
my_mesh.Coherence = True
# write the geofile
my_mesh.writeGeo('my_mesh.geo')
os.system('gmsh my_mesh.geo -2 -o my_mesh.msh')
#os.system('gmsh my_mesh.geo')


x, y, tri, T, V, q = axifem.axiHeatCond('my_mesh.msh', \
            {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})

from matplotlib.pyplot import *
print(V)
print('Rúmmál {}'.format(V['ribba']))
print('Varmaflæði: {:g}'.format(q['ytri'][1]))
print('Hámarkshitastig: {:g}'.format(max(T)))
print('Lágmarkshitastig: {:g}'.format(min(T)))
figure(figsize=(16,3))
tricontourf(x,y,tri,T,20)
colorbar()
axis('equal')
show()


