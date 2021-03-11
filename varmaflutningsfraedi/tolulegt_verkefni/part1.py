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
d1 = 0.01

V_sk = np.pi*(a/2)**2*t_veggur+np.pi*(d/2)**2*L
print(V_sk)
res=()

for i in range(5,13):
    d1 = 1.0/(2**(i+1))
    # create Mesh class instance
    my_mesh = Mesh()

    # create points
    p1 = Entity.Point([0.,0., 0.,d1]) #fyrsti punktur neðri vinstri
    # add point to mesh
    my_mesh.addEntity(p1) 
    #create more points
    p2 = Entity.Point([0.,a/2, 0.,d1])#2. punktur efri vinstri
    my_mesh.addEntity(p2)
    p3 = Entity.Point([t_veggur,a/2, 0.,d1])#3. punktur efri hægri
    my_mesh.addEntity(p3)

    p4 = Entity.Point([t_veggur, (a-d)/2+d-a/2, 0.,d1])#4. punktur niður frá efri hægri
    my_mesh.addEntity(p4)

    p5 = Entity.Point([t_veggur+L,(a-d)/2+d-a/2,0.,d1])#5.punktur endi á ribbu efri
    my_mesh.addEntity(p5)
    p6 = Entity.Point([t_veggur+L,0.,0.,d1])#6. punktur endi á ribbu neðri
    my_mesh.addEntity(p6)
    # create curves
    l1 = Entity.Curve([p1, p2]) #innri bein lína upp
    l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
    l3 = Entity.Curve([p3, p4]) # ytri bein lína upp
    l4 = Entity.Curve([p4, p5]) #ribba bein lína upp
    l5 = Entity.Curve([p5, p6]) #ribba endi
    l6 = Entity.Curve([p6, p1]) #ribba bein lína niðri

    my_mesh.addEntities([l1, l2, l3, l4, l5, l6])


    ll1 = Entity.CurveLoop([l1, l2, l3, l4, l5, l6], mesh=my_mesh)



    s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




    g1 = Entity.PhysicalGroup(name='innri')
    g2 = Entity.PhysicalGroup(name='ytri')
    g3 = Entity.PhysicalGroup(name='ribba')
    g4 = Entity.PhysicalGroup(name='einangrun')
    my_mesh.addEntities([g1, g2, g3, g4])
    g1.addEntities([l1])
    g2.addEntities([l3,l4,l5,l6])
    g4.addEntities([l2])
    g3.addEntities([s1])
    # set max element size
    #my_mesh.Options.Mesh.CharacteristicLengthMax = 0.1

    # adding Coherence option
    my_mesh.Coherence = True
    # write the geofile
    my_mesh.writeGeo('my_mesh.geo')
    os.system('gmsh my_mesh.geo -2 -o my_mesh.msh')
    #os.system('gmsh my_mesh.geo')

    nodes, elem, bc, lab = axifem.readMsh4_1('my_mesh.msh')
    
    x, y, tri, T, V, q = axifem.axiHeatCond('my_mesh.msh', \
                {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})

    
    
    d1 = 1.0/(2**(i+1))
    # create Mesh class instance
    my_mesh = Mesh()

    # create points
    p1 = Entity.Point([0.,0., 0.,d1]) #fyrsti punktur neðri vinstri
    # add point to mesh
    my_mesh.addEntity(p1) 
    #create more points
    p2 = Entity.Point([0.,a/2, 0.,d1])#2. punktur efri vinstri
    my_mesh.addEntity(p2)
    p3 = Entity.Point([t_veggur,a/2, 0.,d1])#3. punktur efri hægri
    my_mesh.addEntity(p3)

    p4 = Entity.Point([t_veggur, 0., 0.,d1])#4. punktur niður frá efri hægri
    my_mesh.addEntity(p4)

    # create curves
    l1 = Entity.Curve([p1, p2]) #innri bein lína upp
    l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
    l3 = Entity.Curve([p3, p4]) # ytri bein lína upp
    l4 = Entity.Curve([p4, p1]) #neðri lína

    my_mesh.addEntities([l1, l2, l3, l4])


    ll1 = Entity.CurveLoop([l1, l2, l3, l4], mesh=my_mesh)



    s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




    g1 = Entity.PhysicalGroup(name='innri')
    g2 = Entity.PhysicalGroup(name='ytri')
    g3 = Entity.PhysicalGroup(name='ribba')
    g4 = Entity.PhysicalGroup(name='einangrun')
    my_mesh.addEntities([g1, g2, g3, g4])
    g1.addEntities([l1])
    g2.addEntities([l3,l4])
    g4.addEntities([l2])
    g3.addEntities([s1])
    # set max element size
    #my_mesh.Options.Mesh.CharacteristicLengthMax = 0.1

    # adding Coherence option
    my_mesh.Coherence = True
    # write the geofile
    my_mesh.writeGeo('my_mesh.geo')
    os.system('gmsh my_mesh.geo -2 -o my_mesh.msh')
    #os.system('gmsh my_mesh.geo')
    
    x1, y1, tri1, T1, V1, q1 = axifem.axiHeatCond('my_mesh.msh', \
            {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})
    
    virkni = q['ytri'][1]/q1['ytri'][1]
    
    re = {'d1: ': d1, 'eiginleikar': {'ribbuvirkni':virkni,'fjöldi þríhyrninga': len(elem)+len(bc)+len(lab)+2,'rúmmál': V['ribba'], 'Varmaflæði': q['ytri'][1],'Hámarkshitastig': max(T),'lágmarkshitastig':min(T)}}
    res+=(re,)

print(res)

my_mesh = Mesh()

# create points
p1 = Entity.Point([0.,0., 0.,d1]) #fyrsti punktur neðri vinstri
# add point to mesh
my_mesh.addEntity(p1) 
#create more points
p2 = Entity.Point([0.,a/2, 0.,d1])#2. punktur efri vinstri
my_mesh.addEntity(p2)
p3 = Entity.Point([t_veggur,a/2, 0.,d1])#3. punktur efri hægri
my_mesh.addEntity(p3)

p4 = Entity.Point([t_veggur, 0., 0.,d1])#4. punktur niður frá efri hægri
my_mesh.addEntity(p4)

# create curves
l1 = Entity.Curve([p1, p2]) #innri bein lína upp
l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
l3 = Entity.Curve([p3, p4]) # ytri bein lína upp
l4 = Entity.Curve([p4, p1]) #neðri lína

my_mesh.addEntities([l1, l2, l3, l4])


ll1 = Entity.CurveLoop([l1, l2, l3, l4], mesh=my_mesh)



s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




g1 = Entity.PhysicalGroup(name='innri')
g2 = Entity.PhysicalGroup(name='ytri')
g3 = Entity.PhysicalGroup(name='ribba')
g4 = Entity.PhysicalGroup(name='einangrun')
my_mesh.addEntities([g1, g2, g3, g4])
g1.addEntities([l1])
g2.addEntities([l3,l4])
g4.addEntities([l2])
g3.addEntities([s1])
# set max element size
#my_mesh.Options.Mesh.CharacteristicLengthMax = 0.1

# adding Coherence option
my_mesh.Coherence = True
# write the geofile
my_mesh.writeGeo('my_mesh.geo')
os.system('gmsh my_mesh.geo -2 -o my_mesh.msh')
#os.system('gmsh my_mesh.geo')

x1, y1, tri1, T1, V1, q1 = axifem.axiHeatCond('my_mesh.msh', \
        {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})

virkni = q['ytri'][1]/q1['ytri'][1]

from matplotlib.pyplot import *
print('Ribbuvirkni {}:'.format(virkni))
print('Rúmmál {}'.format(V['ribba']))
print('Varmaflutningur: {:g}'.format(q['ytri'][1]))
print('Hámarkshitastig: {:g}'.format(max(T)))
print('Lágmarkshitastig: {:g}'.format(min(T)))
figure(figsize=(16,3))



tricontourf(x,y,tri,T,20)
tricontourf(x,-y,tri,T,20)
colorbar()
axis('equal')
text(0.005,0.003,'Varmaflutningur: {:g}W'.format(q['ytri'][1]))
text(0.005,0.004,'Ribbuvirkni: {:g}'.format(virkni))
title('Upphafleg ribba')
xlabel('Lengd [m]')
ylabel('Hæð [m]')
text(0.0175,0.006,'Hitastig í °C')
show()


