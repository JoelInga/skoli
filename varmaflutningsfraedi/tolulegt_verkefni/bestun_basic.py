from py2gmsh import (Mesh, Entity, Field)
import axifem
import os
from scipy.optimize import minimize,NonlinearConstraint
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
d1 = 0.0005
d2 = 0.0005
V_sk = np.pi*(a/2)**2*t_veggur+np.pi*(d/2)**2*(L)
print('initial volume: {}'.format(V_sk))

tolerance = 6e-11 #tolerance for volume
# create Mesh class instance
my_mesh = Mesh()
i = 0

def objective(x, sign=-1.0):
    time.sleep(0.1)
    my_mesh = Mesh()
    filename = 'my_mesh'
    d = x[0]
    a = x[1]
    L = x[2]
    # create points
    p1 = Entity.Point([0., 0., 0.,d1]) #fyrsti punktur neðri vinstri
    # add point to mesh
    my_mesh.addEntity(p1) 
    #create more points
    p2 = Entity.Point([0.,a, 0.,d1])#2. punktur efri vinstri
    my_mesh.addEntity(p2)
    p3 = Entity.Point([t_veggur, a, 0.,d1])#3. punktur efri hægri
    my_mesh.addEntity(p3)

    p4 = Entity.Point([t_veggur, (a-d)/2+d, 0.,d1])#4. punktur niður frá efri hægri
    my_mesh.addEntity(p4)

    p5 = Entity.Point([t_veggur+L,(a-d)/2+d,0.,d1])#5.punktur endi á ribbu efri
    my_mesh.addEntity(p5)
    p6 = Entity.Point([t_veggur+L,(a-d)/2,0.,d1])#6. punktur endi á ribbu neðri
    my_mesh.addEntity(p6)
    p7 = Entity.Point([t_veggur,(a-d)/2,0.,d1])#7. punktur byrjun á ribbu neðri
    my_mesh.addEntity(p7)
    p8 = Entity.Point([t_veggur,0.,0.,d1])#síðasti punktur neðri hægri
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
    #os.system('rm .geo')
    my_mesh.writeGeo('{}.geo'.format(filename))
    os.system('gmsh {}.geo -2 -o {}.msh'.format(filename,filename))
    #os.system('gmsh my_mesh.geo')
    try:
        xu, y, tri, T, V, q = axifem.axiHeatCond('{}.msh'.format(filename), \
                    {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})
        print(sign*q['ytri'][1])
    except:
        return 0
    return sign*q['ytri'][1]

def volume(x):
    d = x[0]
    a = x[1]
    L = x[2]
    v=np.pi*(a/2)**2*t_veggur+np.pi*(d/2)**2*(L)
    
    print(v)
    return v


def constraint1(x):
    x_tmp = [*x]
    v = volume(x_tmp)
    return v-V_sk+tolerance

def constraint2(x):
    x_tmp = [*x]
    v = volume(x_tmp)
    return V_sk-v+tolerance

def constraint3(x): #d ekki minna 0.25mm
    d=x[0]
    return d-0.00025

def constraint4(x): #a ekki minna en 2*d
    d = x[0]
    a = x[1]
    return a-2*d

def constraint5(x): #L stærra en 0
    L = x[2]
    return L



nlc = NonlinearConstraint(volume,V_sk-tolerance,V_sk+tolerance)



cons=({'type': 'ineq',
       'fun': constraint1},
      {'type': 'ineq',
       'fun': constraint2},
      {'type': 'ineq',
       'fun': constraint3},
      {'type': 'ineq',
       'fun': constraint4},
      {'type': 'ineq',
       'fun': constraint5})
x0 = [0.003,0.01,0.01]

#vte = volume(x0)
#print(vte)
sol = minimize(objective,x0,method='SLSQP', constraints = cons)
print(sol)


x, y, tri, T, V, q = axifem.axiHeatCond('my_mesh.msh', \
                {'ribba':16}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})

from matplotlib.pyplot import *
print('Rúmmál {}'.format(V['ribba']))
print('Varmaflæði: {:g}'.format(q['ytri'][1]))
print('Hámarkshitastig: {:g}'.format(max(T)))
print('Lágmarkshitastig: {:g}'.format(min(T)))
figure(figsize=(16,3))
tricontourf(x,y,tri,T,20)
colorbar()
axis('equal')
show()


