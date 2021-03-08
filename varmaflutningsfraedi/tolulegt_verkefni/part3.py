from py2gmsh import (Mesh, Entity, Field)
import axifem
import os
from scipy.optimize import minimize,NonlinearConstraint,Bounds,LinearConstraint,BFGS,SR1
import numpy as np
import time
from random import random
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

tolerance = 6e-10 #tolerance for volume
# create Mesh class instance
my_mesh = Mesh()
i = 0




def objective(x, sign=-1.0):
    #x = [d,a,L,n,t,r]
    my_mesh = Mesh()
    filename = 'my_mesh'
    d = x[0]
    a = x[1]
    L = x[2]
    n = int(x[3])
    t = x[4]
    r = x[5]
    dx = x[6]
    
    if volume(x)>10*V_sk:
        d1=0.05
    elif volume(x)>2*V_sk:
        d1=0.005
    else:
        d1=0.0005
            
    # create points
    p1 = Entity.Point([0.,0., 0.,d1]) #fyrsti punktur neðri vinstri

    p2 = Entity.Point([0.,a/2, 0.,d1])#2. punktur efri vinstri

    p3 = Entity.Point([t_veggur, a/2, 0.,d1])#3. punktur efri hægri
    pe=[]
    
    my_mesh.addEntities([p1,p2,p3])
    
    pe.append(Entity.Point([t_veggur,d/2,0.,d1]))#fyrsti punktur á ribbu

    for i in range(0,n):
        delta_x = dx + i*(dx+t) 
        pe.append(Entity.Point([t_veggur+delta_x,d/2,0.,d1]))
        pe.append(Entity.Point([t_veggur+delta_x,r/2,0.,d1]))
        pe.append(Entity.Point([t_veggur+delta_x+t,r/2,0.,d1]))
        pe.append(Entity.Point([t_veggur+delta_x+t,d/2,0.,d1]))
    

        
    my_mesh.addEntities(pe)
  
  
    p10 = Entity.Point([t_veggur+dx*n+t*n,0.,0.,d1])#síðasti punktur neðri hægri

    my_mesh.addEntities([p10])
    # create curves
    l1 = Entity.Curve([p1, p2]) #innri bein lína upp
    l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
    l3 = Entity.Curve([p3, pe[0]]) # ytri bein lína upp
    le = []
    my_mesh.addEntities([l1, l2, l3])
    
    for i in range(0,len(pe)-1):
        le.append(Entity.Curve([pe[i],pe[i+1]]))
        
    le.append(Entity.Curve([pe[-1],p10])) #klára ferð niður


    l9 = Entity.Curve([p10,p1]) #lína niðri

    
    my_mesh.addEntities(le)

    my_mesh.addEntities([l9])

    lines = []
    lines.append(l1)
    lines.append(l2)
    lines.append(l3)
    for line in le: lines.append(line)
    lines.append(l9)
    
    ll1 = Entity.CurveLoop(lines, mesh=my_mesh)



    s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




    g1 = Entity.PhysicalGroup(name='innri')
    g2 = Entity.PhysicalGroup(name='ytri')
    g3 = Entity.PhysicalGroup(name='ribba')
    g4 = Entity.PhysicalGroup(name='einangrun')
    my_mesh.addEntities([g1, g2, g3, g4])
    g1.addEntities([l1])
    g4.addEntities([l2])
    g2.addEntities([l3])
    for line in le: g2.addEntities([line])
    g3.addEntities([s1])
    # set max element size
    #my_mesh.Options.Mesh.CharacteristicLengthMax = 0.1
    # adding Coherence option
    my_mesh.Coherence = True
    # write the geofile
    #os.system('rm .geo')
    try:
        my_mesh.writeGeo('{}.geo'.format(filename))
        os.system('gmsh {}.geo -2 -o {}.msh'.format(filename,filename))
    except:
        return min(-0.5*(dx+t)*n,-0.5*random())
    #os.system('gmsh my_mesh.geo')
    try:
        xu, y, tri, T, V, q = axifem.axiHeatCond('{}.msh'.format(filename), \
                    {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})
        print(sign*q['ytri'][1])
        if sign*q['ytri'][1]!=sign*q['ytri'][1]:
            return -0.1
    except:
        return min(-0.5*(dx+t)*n,-0.5*random())
    
    return sign*q['ytri'][1]



def volume(x):
    
    d = x[0]
    a = x[1]
    L = x[2]
    n = int(x[3])
    t = x[4]
    r = x[5]
    dx = x[6]
    
    v=np.pi*(a/2)**2*t_veggur#base
    
    v+=np.pi*(d/2)**2*dx*n
    v+=np.pi*(r/2)**2*t*n
    
    print(v)
    return v


def volume1(x):
    d = x[0]
    a = x[1]
    L = x[2]
    n = int(x[3])
    t = x[4]
    r = x[5]
    dx = x[6]
    
    v=np.pi*(a/2)**2*t_veggur#base
    
    v+=np.pi*(d/2)**2*dx*n
    v+=np.pi*(r/2)**2*t*n
    print(v)
    return v-V_sk-tolerance

def volume2(x):
    d = x[0]
    a = x[1]
    L = x[2]
    n = int(x[3])
    t = x[4]
    r = x[5]
    dx = x[6]
    
    v=np.pi*(a/2)**2*t_veggur#base
    
    v+=np.pi*(d/2)**2*dx*n
    v+=np.pi*(r/2)**2*t*n
    print(v)
    return V_sk+tolerance-v



def volume_J(x):
    J = []
    a = x[1]
    L = x[2]
    s = x[3:]
    J.append(0)#fyrir d
    x1=np.pi*x[1]/2*t_veggur
    for i in range(0,len(s)-1):
        x1+=x[2]/(len(s)-1)*np.pi/3*( (x[1]/2-x[3+i])+(x[1]/2-x[3+i]/2-x[3+i+1]/2)+(x[1]/2-x[3+i+1]))
    J.append(x1)
    l1 = 0
    for i in range (0,len(s)-1):
        l1+= 1/3*np.pi*( ((x[1]-2*x[3+i])/2)**2 + (x[1]-2*x[3+i])*(x[1]-2*x[3+i+1])/4 + ((x[1]-2*x[3+i+1])/2)**2 )*1/(len(s)-1)
    J.append(l1)
    s1 = 1/3*np.pi*x[2]/(len(s)-1)*( (-x[1]+2*x[3])+(-x[1]/2+x[3+1]) )
    J.append(s1)
    si=[]
    print(len(s))
    if len(s)>2:
        for i in range(0,len(s)-2):#einum minna en síðasti
            si.append(1/3*np.pi*x[2]/(len(s)-1)*((-x[1]+2*x[4+i])+x[3+i]+(-x[1]+2*x[4+i])+(-x[1]/2+x[4+i+1])))  
        for s_j in si:
            J.append(s_j)
    J.append(1/3*np.pi*x[2]/(len(s)-1)*((-x[1]+2*x[-1])+x[-2]))

    return J



#nlc1 = NonlinearConstraint(volume,V_sk-tolerance,V_sk+tolerance, jac=volume_J,hess=BFGS())

bounds = [[0.00025,0.05],[0.00025,0.05],[0.00025,0.02],[0.,10.1],[0.00025,0.01],[0.00025,0.01],[0.00025,0.01]]
#x0 = [d,a,L,n,t,r,dx]
x0 = [0.003,0.01,0.01,2,0.002,0.0031,(0.01-0.002*2)/2]

cons=()
constraint = {'type': 'ineq', 'fun': volume1}
cons+=(constraint,)
constraint = {'type': 'ineq', 'fun': volume2}
cons+=(constraint,)
constraint = {'type': 'ineq', 'fun': lambda x: x[1]-x[5]}
cons+=(constraint,)
#lcst = LinearConstraint(co,lb,ub)
for factor in range(len(bounds)):
    lower, upper = bounds[factor]
    l = {'type': 'ineq',
         'fun': lambda x, lb=lower, i=factor: x[i] - lb}
    u = {'type': 'ineq',
         'fun': lambda x, ub=upper, i=factor: ub - x[i]}
    cons+=(l,)
    cons+=(u,)

sol = minimize(objective,x0,method='COBYLA', constraints = cons)
#sol = minimize(objective,x0,method='SLSQP',options={'gtol': 1e-6, 'disp': True}, constraints = cons)

a = sol['x'][1]
print(a)
#sol = minimize(objective,x0,method='SLSQP',options={'gtol': 1e-6, 'disp': True}, constraints = cons)
print(sol)

print('initial volume: {}'.format(V_sk))
x, y, tri, T, V, q = axifem.axiHeatCond('my_mesh.msh', \
                {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})
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
print('Rúmmál: {}'.format(V['ribba']))
print('Varmaflæði: {:g}'.format(q['ytri'][1]))
print('Hámarkshitastig: {:g}'.format(max(T)))
print('Lágmarkshitastig: {:g}'.format(min(T)))
figure(figsize=(16,3))
tricontourf(x,y,tri,T,20)
tricontourf(x,-y,tri,T,20)
colorbar()
axis('equal')
title('Bestuð ribba með formbreytingum')
text(0.005,0.003,'Varmaflæði: {:g}W'.format(q['ytri'][1]))
text(0.005,0.004,'Ribbuvirkni: {:g}'.format(virkni))
show()



