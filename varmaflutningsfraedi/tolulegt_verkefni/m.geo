
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


