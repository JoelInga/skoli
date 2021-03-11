from matplotlib.pyplot import *
import numpy as np

data = ({'d1: ': 0.015625, 'eiginleikar': {'ribbuvirkni': 1.9810039343137793, 'fjöldi þríhyrninga': 16, 'rúmmál': 3.0630528372500486e-07, 'Varmaflæði': 0.5145925537601995, 'Hámarkshitastig': 90.05446841057542, 'lágmarkshitastig': 76.32695241139999}}, 
        {'d1: ': 0.0078125, 'eiginleikar': {'ribbuvirkni': 1.962223002129028, 'fjöldi þríhyrninga': 22, 'rúmmál': 3.063052837250049e-07, 'Varmaflæði': 0.5097139527198153, 'Hámarkshitastig': 90.1302826614785, 'lágmarkshitastig': 75.63496499383355}}, 
        {'d1: ': 0.00390625, 'eiginleikar': {'ribbuvirkni': 1.9445832410369384, 'fjöldi þríhyrninga': 36, 'rúmmál': 3.063052837250049e-07, 'Varmaflæði': 0.5058722425913148, 'Hámarkshitastig': 89.73637293645638, 'lágmarkshitastig': 74.90297871584119}}, 
        {'d1: ': 0.001953125, 'eiginleikar': {'ribbuvirkni': 1.9386782630616903, 'fjöldi þríhyrninga': 56, 'rúmmál': 3.063052837250048e-07, 'Varmaflæði': 0.5043942389517646, 'Hámarkshitastig': 89.7144454492285, 'lágmarkshitastig': 74.70076886755044}}, 
        {'d1: ': 0.0009765625, 'eiginleikar': {'ribbuvirkni': 1.934003707881037, 'fjöldi þríhyrninga': 146, 'rúmmál': 3.0630528372500476e-07, 'Varmaflæði': 0.5032011940754662, 'Hámarkshitastig': 89.74451767333588, 'lágmarkshitastig': 74.43627363670225}}, 
        {'d1: ': 0.00048828125, 'eiginleikar': {'ribbuvirkni': 1.931561910084705, 'fjöldi þríhyrninga': 442, 'rúmmál': 3.063052837250048e-07, 'Varmaflæði': 0.5025688235726488, 'Hámarkshitastig': 89.75628959270291, 'lágmarkshitastig': 74.30140907031269}}, 
        {'d1: ': 0.000244140625, 'eiginleikar': {'ribbuvirkni': 1.9305172653515283, 'fjöldi þríhyrninga': 1444, 'rúmmál': 3.0630528372500465e-07, 'Varmaflæði': 0.5022974907646836, 'Hámarkshitastig': 89.76120169459519, 'lágmarkshitastig': 74.2434763137346}}, 
        {'d1: ': 0.0001220703125, 'eiginleikar': {'ribbuvirkni': 1.9300266211948713, 'fjöldi þríhyrninga': 5142, 'rúmmál': 3.06305283725004e-07, 'Varmaflæði': 0.5021698987248999, 'Hámarkshitastig': 89.7633509517971, 'lágmarkshitastig': 74.21560350398747}}, 
        {'d1: ': 6.103515625e-05, 'eiginleikar': {'ribbuvirkni': 1.9298154213926675, 'fjöldi þríhyrninga': 19566, 'rúmmál': 3.0630528372500523e-07, 'Varmaflæði': 0.5021149562804527, 'Hámarkshitastig': 89.76424087177693, 'lágmarkshitastig': 74.2034697056223}}, 
        {'d1: ': 3.0517578125e-05, 'eiginleikar': {'ribbuvirkni': 1.9297295333957039, 'fjöldi þríhyrninga': 76538, 'rúmmál': 3.0630528372500343e-07, 'Varmaflæði': 0.5020926104028152, 'Hámarkshitastig': 89.76458438039423, 'lágmarkshitastig': 74.19853914101589}})

d1 = []
for item in data: 
    d1.append(item['d1: '])
ribbuvirkni = []
tri = []
q = []
for item in data: 
    ribbuvirkni.append(item['eiginleikar']['ribbuvirkni'])
    tri.append(item['eiginleikar']['fjöldi þríhyrninga'])
    q.append(item['eiginleikar']['Varmaflæði'])
    


figure()
plot(tri,ribbuvirkni,'ro-')
title('Ribbuvirkni sem fall af fjölda þríhyrninga')
per=[]
for i in range(0,len(ribbuvirkni)-1):
    tmp = 1-ribbuvirkni[i+1]/ribbuvirkni[i]
    per.append(tmp*100)
xlabel('Fjöldi þríhyrninga')
ylabel('Ribbuvirkni')
savefig('fj_tri_ribbv.eps', format='eps', dpi=1200)
figure()
log = []
for i in range(1,len(tri)):
    lg = np.log10(tri[i])
    log.append(lg)
plot(log,per,'ro-')
xlabel('log(fjöldi þríhyrninga)')
ylabel('Breyting í ribbuvirkni frá síðasta gildi í %')
title('Breyting í ribbuvirkni (%) sem fall af log(fjölda þríhyrninga)')
savefig('fj_tri_ribbv_bre.eps', format='eps', dpi=1200)


figure()
plot(tri,q,'bo-')
title('Varmaflutningur sem fall af fjölda þríhyrninga')
per=[]
for i in range(0,len(q)-1):
    tmp = 1-q[i+1]/q[i]
    per.append(tmp*100)
xlabel('Fjöldi þríhyrninga')
ylabel('Varmaflutningur')
savefig('fj_tri_varm.eps', format='eps', dpi=1200)
figure()
log = []
for i in range(1,len(tri)):
    lg = np.log10(tri[i])
    log.append(lg)
plot(log,per,'bo-')
title('Breyting í varmaflutning (%) sem fall af log(fjölda þríhyrninga)')
xlabel('log(fjöldi þríhyrninga)')
ylabel('Breyting í varmaflutning frá síðasta gildi í %')
savefig('fj_tri_varm_bre.eps', format='eps', dpi=1200)

show()


