import sys
import os
import math
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import Sasha_functions as SF
from tasz import tfora

today = date.today()
print today

z_int_tolerance = 0.5
fhires_tolerance = 0.94
npart_tolerance = 5e4

doplots = 'y'

if (len(sys.argv)<2):
	print 'syntax blah.py haloN'
	sys.exit()

haloN = int(sys.argv[1])

H = SF.read_halo_history(haloN)

Zs = H['redshift']
ID = H['ID']
x = H['x']
y = H['y']
z = H['z']
vx = H['vx']
vy = H['vy']
vz = H['vz']
M = H['M']
Mstar = H['Mstar']
Mgas = H['Mgas']
host = H['host']
fhires = H['fhires']
ngas = H['ngas']
nstar = H['nstar']
npart = H['npart']
Rvir = H['Rvir']
nDM = npart - nstar - ngas
nDM *= fhires

print 'all IDs' 
print ID
print 'all host'
print host




print 'at first, this halo is the following: '
print ID[0], Zs[0], nDM[0], nstar[0], ngas[0], host[0], fhires[0]
print 'at the end, this halo becomes: '
print ID[-1], Zs[-1], nDM[-1], nstar[-1], ngas[-1], host[-1], fhires[-1]

count = 0

type = -1 
while (count < len(ID)):
	if (host[count] != -1):
		print 'this halo stopped being a main halo at ',count, Zs[count],' N = ',SF.N_for_z(Zs[count])
		override = raw_input('do you agree?')
		if (override != 'n'):
			type = 0
			print 'terminating'
			break
		
	count +=1 

if (host[-1] != -1):
	type = 0

if (count > 0):
	count -= 1

good = True






if (nDM[count] < npart_tolerance):
	print 'does not meet npart tolerance'
	print 'doing backup check '
	if (SF.sasha_max(nDM[0:count]) > npart_tolerance):
		newind =  np.where(nDM == SF.sasha_max(nDM[0:count]))[0][0]
		print 'the maximum was reached at count ', newind
		print 'increase newind? '
		override = raw_input('y/n ')
		if (override == 'y'):
			while (newind < count):
				if (nDM[newind]  > npart_tolerance):
					print 'newind now ',newind, nDM[newind]
					newind+=1 
				else: 
					break
			newind -= 1
		
		print ' do you want to change count from ',count,' to ',newind
		override = raw_input('y/n ')
		if (override == 'y'):
			count = newind
		else:
			good = False
	else:
		print 'still does not meet npart tolerance'
		good = False

if (count == (len(ID) -1)):
	print 'this halo was independent for the entire interval'

Zinterval = Zs[0] - Zs[count]

if (Zinterval < z_int_tolerance):
	good = False
	print 'does not meet zint tolerance'
	
if (fhires[count] < fhires_tolerance):
	good = False 
	print 'does not meet fhires tolerance'



print 'at last available output , z=',Zs[count], 'N = ',SF.N_for_z(Zs[count])

print 'ID, z, nDM, nstar, ngas, host, fhires'
print ID[count], Zs[count], nDM[count], nstar[count], ngas[count], host[count], fhires[count]

print 'it was around for redshift interval, ', Zinterval

if (good):
	print 'this halo is worth doing! '
	if (count != (len(ID)-1)):
		print 'but please stop at z=',Zs[count], ' N = ',SF.N_for_z(Zs[count])
	print 'updating good_halos.txt'
	if (os.path.isfile('good_halos.txt')):
		f = np.loadtxt('good_halos.txt', ndmin=2)
		haloIDs = f[:,3]
		if (ID[count] in haloIDs):
			print 'this already has an entry in good_halos.txt'
			print 'keep going? '
			override = raw_input('y/n ')
			if (override == 'n'):
				print 'ok not including this branch '
				good = False
	if (good):
		g = open('good_halos.txt', 'a')
		line = [haloN,  fhires[count], nDM[count], ID[count], Zs[0], Zs[count], SF.N_for_z(Zs[0]), SF.N_for_z(Zs[count]), type, nstar[count]]
		lala = SF.line_to_string(line)
		g.write(lala)
		g.close()
else:
	print 'do not do this halo!' 

print 'updating tried_halos.txt'
g = open('tried_halos.txt', 'a')
line = [haloN,  fhires[count], nDM[count], ID[count],  Zs[0], Zs[count], SF.N_for_z(Zs[0]), SF.N_for_z(Zs[count]), good, type, nstar[count]]
lala = SF.line_to_string(line)
g.write(lala)
g.close()

if (doplots != 'y'):
	print 'skipping plots '
	sys.exit()



####################################################################
cuts = Zs >= Zs[count]
As = 1.0 / (1.0 + Zs)

haloN = str(haloN)

fbar = (Mgas + Mstar)/M
fgas = Mgas/M
fstar = Mstar/M

print 'doing plots for the valid interval'
fig1 = plt.figure(figsize=(10, 6))
plt.xlabel('z')
plt.ylabel('fraction of total mass ')
plt.plot(Zs[cuts], fbar[cuts], color='k', lw=3)
plt.plot(Zs[cuts], fgas[cuts], color='b', ls='--',lw=3)
plt.plot(Zs[cuts], fstar[cuts], color='r', ls='--', lw=3)
plt.gca().invert_xaxis()

plt.legend(['baryons','gas','stars'], loc = 'best', prop={'size':10})
plt.savefig('./' + str(today) + 'baryonfractionN'+haloN+'.pdf')
plt.clf()


fig2  = plt.figure(figsize=(10, 6))
plt.xlabel('z')
plt.ylabel(r'$R_{vir}$ (comoving kpc/h)')
plt.plot(Zs[cuts], Rvir[cuts], color='k', lw=2)
plt.gca().invert_xaxis()

plt.savefig('./' + str(today) + 'NRvir_evo'+haloN+'.pdf')
plt.clf()




fig3  = plt.figure(figsize=(10, 6))
plt.xlabel('z')
plt.ylabel(r'$M_{vir}$ ($M_{\odot}$)')
plt.yscale('log')
plt.plot(Zs[cuts], M[cuts], color='k', lw=2)
plt.gca().invert_xaxis()

plt.savefig('./' + str(today) + 'NTotalMass_evo'+haloN+'.pdf')
plt.clf()

fig4  = plt.figure(figsize=(10, 6))
plt.xlabel('z')
plt.ylabel(r'$M_*$ ($M_{\odot}$)')
plt.yscale('log')
if (np.sum(Mstar[cuts]) < 0.1):
	Mstar[cuts]+=1
plt.plot(Zs[cuts], Mstar[cuts], color='k', lw=2)
plt.gca().invert_xaxis()

plt.savefig('./' + str(today) + 'StellarMass_evo'+haloN+'.pdf')
plt.clf()

fig5  = plt.figure(figsize=(10, 6))
plt.xlabel('z')
plt.ylabel(r'$M_{gas}$ ($M_{\odot}$)')
plt.yscale('log')
if (np.sum(Mgas[cuts]) < 0.1):
	Mgas[cuts]+=1
plt.plot(Zs[cuts], Mgas[cuts], color='k', lw=2)
plt.gca().invert_xaxis()
plt.savefig('./' + str(today) + 'GasMass_evo'+haloN+'.pdf')
plt.clf()



dx = [0]
dy = [0]
dz = [0]

countb = 1

while (countb < len(x)):
	#print 'lalala ',countb
	dx.append(x[countb] - x[countb-1])
	dy.append(y[countb] - y[countb-1])
	dz.append(z[countb] - z[countb-1])
	countb+=1
	
print 'lala ', len(dx), len(x)

countbb = 1
ddx = [0]
ddy = [0]
ddz = [0]
while (countbb < len(x)):
	ddx.append(dx[countbb] - dx[countbb-1])
	ddy.append(dy[countbb] - dy[countbb-1])
	ddz.append(dz[countbb] - dz[countbb-1])
	countbb+=1


dvx = [0]
dvy = [0]
dvz = [0]
countb=1
while (countb < len(x)):
	dvx.append(vx[countb] - vx[countb-1])
	dvy.append(vy[countb] - vy[countb-1])
	dvz.append(vz[countb] - vz[countb-1])
	countb+=1

dx = np.array(dx)
dy = np.array(dy)
dz = np.array(dz)

ddx = np.array(ddx)
ddy = np.array(ddy)
ddz = np.array(ddz)

dvx = np.array(dvx)
dvy = np.array(dvy)
dvz = np.array(dvz)


fig6 = plt.figure(figsize=(10, 6))
plt.xlabel('redshift')
plt.ylabel('dx, dy, dz (kpc/h)')
plt.plot(Zs[cuts], dx[cuts], color='k', lw=2)
plt.plot(Zs[cuts], dy[cuts], color='b', lw=2)
plt.plot(Zs[cuts], dz[cuts], color='r',  lw=2)
plt.gca().invert_xaxis()
plt.legend(['dx','dy','dz'], loc = 'best', prop={'size':10})
plt.savefig('./' + str(today) + 'xyzN'+haloN+'.pdf')
plt.clf()



problemcut = (np.fabs(dx[cuts])>5 ) + (np.fabs(dy[cuts])>5) + (np.fabs(dz[cuts])>5) 



print 'problematic? snapshots/redshifts positions'
count = 0
while (count < len(Zs[cuts][problemcut])):
	print Zs[cuts][problemcut][count], SF.N_for_z(Zs[cuts][problemcut][count]), dx[cuts][problemcut][count], dy[cuts][problemcut][count], dz[cuts][problemcut][count]
	count += 1

problemcut2 =  (np.fabs(ddx[cuts]) + np.fabs(ddy[cuts]) + np.fabs(ddz[cuts]))>5  
print 'real problematic snapshots d(dx)'
count = 0
while (count < len(Zs[cuts][problemcut2])):
	print Zs[cuts][problemcut2][count], SF.N_for_z(Zs[cuts][problemcut2][count]), dx[cuts][problemcut2][count], dy[cuts][problemcut2][count], dz[cuts][problemcut2][count]
	count += 1

fig7 = plt.figure(figsize=(10, 6))
plt.xlabel('redshift')
plt.ylabel('dvx, dvy, dvz (km/s)')
plt.plot(Zs[cuts], dvx[cuts], color='k', lw=2)
plt.plot(Zs[cuts], dvy[cuts], color='b', lw=2)
plt.plot(Zs[cuts], dvz[cuts], color='r',  lw=2)
plt.gca().invert_xaxis()
plt.legend(['dvx','dvy','dvz'], loc = 'best', prop={'size':10})
plt.savefig('./' + str(today) + 'vxvyvzN'+haloN+'.pdf')

problemcut = (np.fabs(dvx[cuts])>1.5 ) + (np.fabs(dvy[cuts])>1.5) + (np.fabs(dvz[cuts])>1.5) 

print 'problematic? snapshots/redshifts velocity'

while (count < len(Zs[cuts][problemcut])):
	print Zs[cuts][problemcut][count], SF.N_for_z(Zs[cuts][problemcut][count]), dvx[cuts][problemcut][count], dvy[cuts][problemcut][count], dvz[cuts][problemcut][count]
	count += 1
