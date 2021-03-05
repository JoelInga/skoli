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
    s = x[3:]

    ''''
    for i in s:
        if 0.0>a-2*i:
            return 0
    '''     
    
    # create points
    p1 = Entity.Point([0., 0., 0.,d1]) #fyrsti punktur neðri vinstri

    p2 = Entity.Point([0.,a, 0.,d1])#2. punktur efri vinstri

    p3 = Entity.Point([t_veggur, a, 0.,d1])#3. punktur efri hægri
    pe=[]
    pn=[]
    
    my_mesh.addEntities([p1,p2,p3])

    for j in range(0,len(s)):
        pe.append(Entity.Point([t_veggur+L/(len(s)-1)*j,a-s[j],0.,d1]))
        
        pn.append(Entity.Point([t_veggur+L/(len(s)-1)*(len(s)-j-1),s[len(s)-1-j],0.,d1]))
        
    my_mesh.addEntities(pe)
    my_mesh.addEntities(pn)
  
  
    p10 = Entity.Point([t_veggur,0.,0.,d1])#síðasti punktur neðri hægri

    my_mesh.addEntities([p10])
    # create curves
    l1 = Entity.Curve([p1, p2]) #innri bein lína upp
    l2 = Entity.Curve([p2, p3]) # efri hlið einangrun
    l3 = Entity.Curve([p3, pe[0]]) # ytri bein lína upp
    le = []
    ln = []
    my_mesh.addEntities([l1, l2, l3])
    
    #for i in range(0,len(pe)-1):
    
    le.append(Entity.Spline(pe[:]))
        
    ln.append(Entity.Curve([pe[-1],pn[0]]))
    #for i in range(0,len(pn)-1):
        
    ln.append(Entity.Spline(pn[:]))
    
    
    l9 = Entity.Curve([pn[-1],p10])
    l10 = Entity.Curve([p10,p1]) #einangrun neðri
    
    my_mesh.addEntities(le)
    my_mesh.addEntities(ln)
    my_mesh.addEntities([l9])
    my_mesh.addEntities([l10])

    lines = []
    lines.append(l1)
    lines.append(l2)
    lines.append(l3)
    for line in le: lines.append(line)
    for line in ln: lines.append(line)
    lines.append(l9)
    lines.append(l10)
    
    ll1 = Entity.CurveLoop(lines, mesh=my_mesh)



    s1 = Entity.PlaneSurface([ll1], mesh=my_mesh)




    g1 = Entity.PhysicalGroup(name='innri')
    g2 = Entity.PhysicalGroup(name='ytri')
    g3 = Entity.PhysicalGroup(name='ribba')
    g4 = Entity.PhysicalGroup(name='einangrun')
    my_mesh.addEntities([g1, g2, g3, g4])
    g1.addEntities([l1])
    g2.addEntities([l3])
    g2.addEntities(le)
    g2.addEntities(ln)
    g2.addEntities([l9])
    g4.addEntities([l2,l10])
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
    s = x[3:]
    v=np.pi*(a/2)**2*t_veggur#base
    
    for i in range (0,len(s)-1):
        v+= 1/3*np.pi*( ((a-2*s[i])/2)**2 + (a-2*s[i])*(a-2*s[i+1])/4 + ((a-2*s[i+1])/2)**2 )*L/(len(s)-1)

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
    d = x[0]

    return d-0.00025

def constraint4(x):

    L = x[2]
    return L

def constraint5(x):
    a = x[1]
    return a
#byrja constr
def constraint6(x):
    a = x[1]
    s = x[3]
    return a-2*s
def constraint7(x):
    s = x[3]
    return s


def constraint8(x):
    a = x[1]
    s = x[4]
    return a-2*s
def constraint9(x):
    s = x[4]
    return s


def constraint10(x):
    a = x[1]
    s = x[5]
    return a-2*s
def constraint11(x):
    s = x[5]
    return s

def constraint12(x):
    a = x[1]
    s = x[6]
    return a-2*s
def constraint13(x):
    s = x[6]
    return s
def constraint14(x):
    a = x[1]
    s = x[7]
    return a-2*s
def constraint15(x):
    s = x[7]
    return s
'''

s0 = x0[3:]
f = []
for i in range(0,len(s0)):
    num = 3+i
    def c(x,i=i):
        a = x[1]
        s = x[num]

        return a-2*s
    f.append(c)


f_min = []
for i in range(0,len(s0)):
    num = 3+i
    def de(x,i=i):
        return x[num]-0.0001
    f_min.append(de)

print(f)
print(f_min)'''
x0 = [0.003,0.01,0.01,0.001, 0.002,0.003,0.004]
cons=({'type': 'ineq',
       'fun': constraint1},
      {'type': 'ineq',
       'fun': constraint2},
      {'type': 'ineq',
       'fun': constraint3},
      {'type': 'ineq',
       'fun': constraint4},
      {'type': 'ineq',
       'fun': constraint5},
      {'type': 'ineq',
       'fun': constraint6},
      {'type': 'ineq',
       'fun': constraint7},
      {'type': 'ineq',
       'fun': constraint8},
      {'type': 'ineq',
       'fun': constraint9},
      {'type': 'ineq',
       'fun': constraint10},
      {'type': 'ineq',
       'fun': constraint11},
      {'type': 'ineq',
       'fun': constraint12},
      {'type': 'ineq',
       'fun': constraint13})


'''

      {'type': 'ineq',
       'fun': constraint12},
      {'type': 'ineq',
       'fun': constraint13},
      {'type': 'ineq',
       'fun': constraint14},
      {'type': 'ineq',
       'fun': constraint15}
       
       
       
temp = list(cons)
for function in f:
    
    temp.append({'type': 'ineq',
       'fun': function})
               
for function in f_min:
    temp.append({'type': 'ineq',
       'fun': function})
   
cons = tuple(temp)
#vte = volume(x0)
#print(vte)
print(cons)
 '''
sol = minimize(objective,x0,method='SLSQP',options={'gtol': 1e-6, 'disp': True}, constraints = cons)
print(sol)

print('initial volume: {}'.format(V_sk))
x, y, tri, T, V, q = axifem.axiHeatCond('my_mesh.msh', \
                {'ribba':k}, {'ytri':(h_utan,-h_utan*T_inf_utan),'innri':(h_innan,-h_innan*T_inf_innan),'einangrun':(0,0)})

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


