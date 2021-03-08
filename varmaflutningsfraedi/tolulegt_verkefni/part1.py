from py2gmsh import (Mesh, Entity, Field)
import axifem
import os
from scipy.optimize import minimize
import numpy as np
import time
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
V_sk = np.pi*(a/2)**2*t_veggur+np.pi*(d/2)**2*L
print(V_sk)


# create Mesh class instance
my_mesh = Mesh()

# create points
p1 = Entity.Point([0., -a/2, 0.,d1]) #fyrsti punktur neðri vinstri
# add point to mesh
my_mesh.addEntity(p1) 
#create more points
p2 = Entity.Point([0.,a-a/2, 0.,d1])#2. punktur efri vinstri
my_mesh.addEntity(p2)
p3 = Entity.Point([t_veggur, a-a/2, 0.,d1])#3. punktur efri hægri
my_mesh.addEntity(p3)

p4 = Entity.Point([t_veggur, (a-d)/2+d-a/2, 0.,d1])#4. punktur niður frá efri hægri
my_mesh.addEntity(p4)

p5 = Entity.Point([t_veggur+L,(a-d)/2+d-a/2,0.,d1])#5.punktur endi á ribbu efri
my_mesh.addEntity(p5)
p6 = Entity.Point([t_veggur+L,((a-d)/2-a/2)/2,0.,d1])#6. punktur endi á ribbu neðri
my_mesh.addEntity(p6)
p7 = Entity.Point([t_veggur,(a-d)/2-a/2,0.,d1])#7. punktur byrjun á ribbu neðri
my_mesh.addEntity(p7)
p8 = Entity.Point([t_veggur,-a/2,0.,d1])#síðasti punktur neðri hægri
my_mesh.addEntity(p8)
# create curves
l1 = Entity.Curve([p1, p2]) #innri bein lína upp
l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
l3 = Entity.Curve([p3, p4]) # ytri bein lína upp
l4 = Entity.Curve([p4, p5]) #ribba bein lína upp
l5 = Entity.Curve([p5, p6]) #ribba endi
l6 = Entity.Curve([p6, p7]) #ribba bein lína niðri
l7 = Entity.Curve([p7, p8]) # ytri bein lína niðri
l8 = Entity.Curve([p8, p1]) #neðri hlið einangrun

my_mesh.addEntities([l1, l2, l3, l4, l5, l6, l7, l8])


ll1 = Entity.CurveLoop([l1, l2, l3, l4, l5, l6, l7, l8], mesh=my_mesh)



s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




g1 = Entity.PhysicalGroup(name='innri')
g2 = Entity.PhysicalGroup(name='ytri')
g3 = Entity.PhysicalGroup(name='ribba')
g4 = Entity.PhysicalGroup(name='einangrun')
my_mesh.addEntities([g1, g2, g3, g4])
g1.addEntities([l1])
g2.addEntities([l3,l4,l5,l6,l7])
g4.addEntities([l2,l8])
g3.addEntities([s1])
# set max element size
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


