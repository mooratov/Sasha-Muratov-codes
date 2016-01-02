#re-centers halos according to kernel density estimator routine. Useful for halos with significant baryonic and DM components. Uses amiga halo finder input

import matplotlib
matplotlib.use('Agg')
import sys
import os
import math
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import Sasha_functions as SF
import graphics_library as GL
from tasz import tfora
from readsnap import readsnap
use_fixed_halos = 0


today = date.today()
today = '2015-05-15'
print today

fake = False
multifile = True
adaptive_Rstar = True
Inner_SF_thresh = 0.2
adaptivejump = True
jumpfac = 0.02
doplots = True
DMcheck = True

def calc_stellar_dens(Stars, haloX, haloY, haloZ, theRvir, rvirFac, boxsize):
	dists = SF.calcdist2(haloX, haloY, haloZ, Stars['p'][:,0], Stars['p'][:,1], Stars['p'][:,2], boxsize)
	cut = dists < theRvir * rvirFac
	StellarMass = np.sum(Stars['m'][cut])
	StellarDens = StellarMass / ((4.0/3.0) * 3.14159 * (theRvir*  rvirFac)**3.0)
	return StellarDens

if (len(sys.argv)<2):
	print 'syntax blah.py haloN [Zmin]'
	sys.exit()

haloN = int(sys.argv[1])
halo_number_string = str(haloN)
if (haloN < 10):
	halo_number_string='0'+halo_number_string

H = SF.read_halo_history(haloN, rod=use_fixed_halos)

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
Vsig = H['Vsig']

if (len(sys.argv)>2):
	Zlimit = float(sys.argv[2])
	print 'using Z limit', Zlimit
	cut = Zs > Zlimit
	ID = ID[cut]
	x = x[cut]
	y = y[cut]
	z = z[cut]
	vx = vx[cut]
	vy = vy[cut]
	vz = vz[cut]
	M = M[cut]
	Mstar = Mstar[cut]
	Mgas = Mgas[cut]
	host = host[cut]
	fhires = fhires[cut]
	ngas = ngas[cut]
	nstar = nstar[cut]
	npart = npart[cut]
	Rvir = Rvir[cut]
	Vsig = Vsig[cut]
	Zs = Zs[cut]

Ns = []
for theZ in Zs:
   N = SF.N_for_z(theZ)
   Ns.append(N)



count = 0
jump_use = []

iSinX = []
iSinY = []
iSinZ = []

ShalX = []
ShalY = []
ShalZ = []
TheLen = len(Zs)

DinX = []
DinY = []
DinZ = []

iDhalX = []
iDhalY = []
iDhalZ = []




while (count < TheLen and not fake):
	theN = Ns[count]
	print 'reading particle type 4 for ',theN
	if (theN < 10):
		Nsnapstring = '00'+str(theN)
	elif (theN < 100):
		Nsnapstring = '0'+str(theN)
	else:
		Nsnapstring = str(theN)
	the_snapdir = './hdf5/'
	if (multifile):
		the_snapdir = './snapdir_'+Nsnapstring+'/'
	the_prefix ='snapshot'
	the_suffix = '.hdf5'
	S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
	thetime = S['header'][2]
	redshift = S['header'][3]
	boxsize = S['header'][9]
	omega_matter = S['header'][10]
	omega_L = S['header'][11]
	h = S['header'][12]
	AHF_haloX = x[count]
	AHF_haloY = y[count] 
	AHF_haloZ = z[count]
	AHF_Rvir = Rvir[count]
	AHF_a = 1.0 / (1.0 + Zs[count])
	
	
	if (adaptive_Rstar):
		Rstars = Inner_SF_thresh*Rvir[count]
		RstarsPhys = Rstars * AHF_a / h
	Sdists = SF.calcdist2(AHF_haloX, AHF_haloY, AHF_haloZ, S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize)
	Sinner = Sdists<Rstars
	S_COM = SF.gaussian_KDE_center(S['p'][:,0][Sinner], S['p'][:,1][Sinner], S['p'][:,2][Sinner], downsample=True)	
	Shalo = Sdists < AHF_Rvir
	S_COM_grand = SF.gaussian_KDE_center(S['p'][:,0][Shalo], S['p'][:,1][Shalo], S['p'][:,2][Shalo], downsample=True)	
	iSinX.append(S_COM[0])
	iSinY.append(S_COM[1])
	iSinZ.append(S_COM[2])
	
	ShalX.append(S_COM_grand[0])
	ShalY.append(S_COM_grand[1])
	ShalZ.append(S_COM_grand[2])	

	D = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix)
	Ddists = SF.calcdist2(AHF_haloX, AHF_haloY, AHF_haloZ, D['p'][:,0], D['p'][:,1], D['p'][:,2], boxsize)
	Dinner = Ddists < Rstars	
	Dhalo = Ddists < AHF_Rvir
	D_COM_inner = 	SF.gaussian_KDE_center(D['p'][:,0][Dinner], D['p'][:,1][Dinner], D['p'][:,2][Dinner], downsample=True)
	D_COM_grand = 	SF.gaussian_KDE_center(D['p'][:,0][Dhalo], D['p'][:,1][Dhalo], D['p'][:,2][Dhalo], downsample=True)

	iDhalX.append( D_COM_grand[0])
	iDhalY.append( D_COM_grand[1])
	iDhalZ.append( D_COM_grand[2])

	DinX.append(D_COM_inner[0])
	DinY.append(D_COM_inner[1])
	DinZ.append(D_COM_inner[2])
	count += 1
	del(S)
	del(D)


if (fake):
	iSinX = x
	iSinY = y
	iSinZ = z
	ShalX = x
	ShalY = y
	ShalZ = z
	iDhalX = x
	iDhalY = y 
	iDhalZ = z
	DinX = x
	DinY = y
	DinZ = z

dx = [0]
dy = [0]
dz = [0]

countb = 1

DhalX = iDhalX
DhalY = iDhalY
DhalZ = iDhalZ

SinX = iSinX
SinY = iSinY
SinZ = iSinZ

print 'searching for jumps with dx' 
while (countb < len(x)):
	#print 'lalala ',countb
	dx.append(x[countb] - x[countb-1])
	dy.append(y[countb] - y[countb-1])
	dz.append(z[countb] - z[countb-1])
	countb+=1

countb = 1
ddx = [0]
ddy = [0]
ddz = [0]
print 'searching for jumps with ddx' 
while (countb < len(x)):
	ddx.append(dx[countb] - dx[countb-1])
	ddy.append(dy[countb] - dy[countb-1])
	ddz.append(dz[countb] - dz[countb-1])
	countb +=1
	
	
Sindx = [0]
Sindy = [0]
Sindz = [0]

countb = 1

print 'searching for jumps with Sin dx' 
while (countb < len(x)):
	#print 'lalala ',countb
	Sindx.append(SinX[countb] - SinX[countb-1])
	Sindy.append(SinY[countb] - SinY[countb-1])
	Sindz.append(SinZ[countb] - SinZ[countb-1])
	countb+=1

countb = 1
Sinddx = [0]
Sinddy = [0]
Sinddz = [0]

print 'searching for jumps with Sin ddx' 
while (countb < len(x)):
	Sinddx.append(Sindx[countb] - Sindx[countb-1])
	Sinddy.append(Sindy[countb] - Sindy[countb-1])
	Sinddz.append(Sindz[countb] - Sindz[countb-1])
	countb +=1

OSinX = np.copy(SinX)
OSinY = np.copy(SinY)
OSinZ = np.copy(SinZ)

OSindx = np.copy(Sindx)
OSindy = np.copy(Sindy)
OSindz = np.copy(Sindz)

OSinddx = np.copy(Sinddx)
OSinddy = np.copy(Sinddy)
OSinddz = np.copy(Sinddz)

finalx = [ShalX[0]]
finaly = [ShalY[0]]
finalz = [ShalZ[0]]

finaldx = [Sindx[0]]
finaldy = [Sindy[0]]
finaldz = [Sindz[0]]

finalddx = [0]
finalddy = [0]
finalddz = [0]

countb=1
while (countb < len(SinX)):
	if (adaptivejump):
		jumpthresh =Rvir[count-1]*jumpfac
	ddsum = np.fabs(Sinddx[countb]) + np.fabs(Sinddy[countb]) + np.fabs(Sinddz[countb])
	if (ddsum>jumpthresh):
		print countb, 'uh oh there is a jump with',jumpthresh, 'dd sum',ddsum
		print 'lets see if we can fix it'	
		trialSx = ShalX[countb]
		trialSy = ShalY[countb]
		trialSz = ShalZ[countb]
		
		trial_dSx = trialSx - finalx[countb-1]
		trial_dSy = trialSy - finaly[countb-1]
		trial_dSz = trialSz - finalz[countb-1]
		
		trial_ddSx = trial_dSx -  finaldx[countb-1]
		trial_ddSy = trial_dSy -  finaldy[countb-1]
		trial_ddSz = trial_dSz -  finaldz[countb-1]
		
		trialddsum = np.fabs(trial_ddSx) + np.fabs(trial_ddSy) + np.fabs(trial_ddSz)
		
		trialSx2 = x[countb]
		trialSy2 = y[countb]
		trialSz2 = z[countb]
		
		trial_dSx2 = trialSx2 - finalx[countb-1]
		trial_dSy2 = trialSy2 - finaly[countb-1]
		trial_dSz2 = trialSz2 - finalz[countb-1]
		
		trial_ddSx2 = trial_dSx2 -  finaldx[countb-1]
		trial_ddSy2 = trial_dSy2 -  finaldy[countb-1]
		trial_ddSz2 = trial_dSz2 -  finaldz[countb-1]				
		trialddsum2 = np.fabs(trial_ddSx2) + np.fabs(trial_ddSy2) + np.fabs(trial_ddSz2)

		trialSx3 = DhalX[countb]
		trialSy3 = DhalY[countb]
		trialSz3 = DhalZ[countb]

		trial_dSx3 = trialSx3 - finalx[countb-1]
		trial_dSy3 = trialSy3 - finaly[countb-1]
		trial_dSz3 = trialSz3 - finalz[countb-1]

		trial_ddSx3 = trial_dSx3 -  finaldx[countb-1]
		trial_ddSy3 = trial_dSy3 -  finaldy[countb-1]
		trial_ddSz3 = trial_dSz3 -  finaldz[countb-1]				
		trialddsum3 = np.fabs(trial_ddSx3) + np.fabs(trial_ddSy3) + np.fabs(trial_ddSz3)

		trialSx4 = DinX[countb]
		trialSy4 = DinY[countb]
		trialSz4 = DinZ[countb]

		trial_dSx4 = trialSx4 - finalx[countb-1]
		trial_dSy4 = trialSy4 - finaly[countb-1]
		trial_dSz4 = trialSz4 - finalz[countb-1]

		trial_ddSx4 = trial_dSx4 -  finaldx[countb-1]
		trial_ddSy4 = trial_dSy4 -  finaldy[countb-1]
		trial_ddSz4 = trial_dSz4 -  finaldz[countb-1]				
		trialddsum4 = np.fabs(trial_ddSx4) + np.fabs(trial_ddSy4) + np.fabs(trial_ddSz4)

		print 'with all halo stars, i get :',trialddsum,' is this an improvement from ',ddsum, ' ? we will also try original coordinates ',trialddsum2, 'and inner stars? ',trialddsum3, ' how about inner DM?',trialddsum4
		winner = min([ddsum, trialddsum, trialddsum2, trialddsum3, trialddsum4])
		
		if (winner == trialddsum):
			print 'i guess total stellar center is improvement '
			finalx.append(trialSx)
			finaly.append(trialSy)
			finalz.append(trialSz)
			
			finaldx.append(trial_dSx)
			finaldy.append(trial_dSy)
			finaldz.append(trial_dSz)
			
			finalddx.append(trial_ddSx)
			finalddy.append(trial_ddSy)
			finalddz.append(trial_ddSz)
			
			if (countb+1 < len(SinX)):
				print 'adjusting next steps' 
				Sindx[countb+1] = SinX[countb+1] - finalx[countb]
				Sinddx[countb+1] = Sindx[countb+1] - finaldx[countb]
				Sindy[countb+1] = SinY[countb+1] - finaly[countb]
				Sinddy[countb+1] = Sindy[countb+1] - finaldy[countb]
				Sindz[countb+1] = SinZ[countb+1] - finalz[countb]
				Sinddz[countb+1] = Sindz[countb+1] - finaldz[countb]
		elif (winner == trialddsum2):
			print 'with original, i get :',trialddsum2,' & this is an improvement from ',ddsum
			finalx.append(trialSx2)
			finaly.append(trialSy2)
			finalz.append(trialSz2)
			
			finaldx.append(trial_dSx2)
			finaldy.append(trial_dSy2)
			finaldz.append(trial_dSz2)
			
			finalddx.append(trial_ddSx2)
			finalddy.append(trial_ddSy2)
			finalddz.append(trial_ddSz2)
			
			if (countb+1 < len(SinX)):
				print 'adjusting next steps' 
				Sindx[countb+1] = SinX[countb+1] - finalx[countb]
				Sinddx[countb+1] = Sindx[countb+1] - finaldx[countb]
				Sindy[countb+1] = SinY[countb+1] - finaly[countb]
				Sinddy[countb+1] = Sindy[countb+1] - finaldy[countb]
				Sindz[countb+1] = SinZ[countb+1] - finalz[countb]
				Sinddz[countb+1] = Sindz[countb+1] - finaldz[countb]	
		elif (winner == trialddsum3):
			print 'with inner stars, i get :',trialddsum3,' & this is an improvement from ',ddsum
			finalx.append(trialSx3)
			finaly.append(trialSy3)
			finalz.append(trialSz3)
	
			finaldx.append(trial_dSx3)
			finaldy.append(trial_dSy3)
			finaldz.append(trial_dSz3)
	
			finalddx.append(trial_ddSx3)
			finalddy.append(trial_ddSy3)
			finalddz.append(trial_ddSz3)
			if (countb+1 < len(SinX)):
				print 'adjusting next steps' 
				Sindx[countb+1] = SinX[countb+1] - finalx[countb]
				Sinddx[countb+1] = Sindx[countb+1] - finaldx[countb]
				Sindy[countb+1] = SinY[countb+1] - finaly[countb]
				Sinddy[countb+1] = Sindy[countb+1] - finaldy[countb]
				Sindz[countb+1] = SinZ[countb+1] - finalz[countb]
				Sinddz[countb+1] = Sindz[countb+1] - finaldz[countb]
		elif (winner == trialddsum4):
			print 'with inner DM, i get :',trialddsum4,' & this is an improvement from ',ddsum
			finalx.append(trialSx4)
			finaly.append(trialSy4)
			finalz.append(trialSz4)
	
			finaldx.append(trial_dSx4)
			finaldy.append(trial_dSy4)
			finaldz.append(trial_dSz4)
	
			finalddx.append(trial_ddSx4)
			finalddy.append(trial_ddSy4)
			finalddz.append(trial_ddSz4)
	
			if (countb+1 < len(SinX)):
				print 'adjusting next steps' 
				Sindx[countb+1] = SinX[countb+1] - finalx[countb]
				Sinddx[countb+1] = Sindx[countb+1] - finaldx[countb]
				Sindy[countb+1] = SinY[countb+1] - finaly[countb]
				Sinddy[countb+1] = Sindy[countb+1] - finaldy[countb]
				Sindz[countb+1] = SinZ[countb+1] - finalz[countb]
				Sinddz[countb+1] = Sindz[countb+1] - finaldz[countb]								
		else:
			print 'guess theres nothing i can do'
			finalx.append(SinX[countb])
			finaly.append(SinY[countb])
			finalz.append(SinZ[countb])
			finaldx.append(Sindx[countb])
			finaldy.append(Sindy[countb])
			finaldz.append(Sindz[countb])
			finalddx.append(Sinddx[countb])
			finalddy.append(Sinddy[countb])
			finalddz.append(Sinddz[countb])
	else:
		finalx.append(SinX[countb])
		finaly.append(SinY[countb])
		finalz.append(SinZ[countb])
		finaldx.append(Sindx[countb])
		finaldy.append(Sindy[countb])
		finaldz.append(Sindz[countb])
		finalddx.append(Sinddx[countb])
		finalddy.append(Sinddy[countb])
		finalddz.append(Sinddz[countb])
	countb += 1
		

extra_artifact = ''
if (use_fixed_halos == 1):
	extra_artifact ='_rod'
if (use_fixed_halos == 2):
	extra_artifact = '_wod'
finname = 'halos/halo_000'+halo_number_string+extra_artifact+'.dat'
f = open(finname, 'r')
foutname = 'halos/halo_000'+halo_number_string+'_ffzod.dat'
print 'foutname ',foutname

strang1 = f.readline()
strang2 = f.readline()	

dars = f.readlines()
f.close()
gork = open(foutname, 'w')
gork.write(strang1)
gork.write(strang2)
count = 0
for line in dars:
	xsd = line.split()
	xsd[6] = finalx[count]
	xsd[7] = finaly[count]
	xsd[8] = finalz[count]
	myline = SF.line_to_string(xsd, newline='y', tabs='y')
	#print myline
	gork.write(myline)
	count +=1
	if (count >= len(finalx)):
		break
gork.close()

OSinddsum = np.fabs(OSinddx) + np.fabs(OSinddy) + np.fabs(OSinddz)
finalddsum = np.fabs(finalddx) + np.fabs(finalddy) + np.fabs(finalddz)
initddsum = np.fabs(ddx) + np.fabs(ddy) + np.fabs(ddz)

fig3 = plt.figure(figsize=(10, 6))
plt.xlabel('z')
plt.ylabel('ddx, ddy, ddz (kpc/h)')
plt.plot(Zs, initddsum, color='r', ls ='-', lw=1)
plt.plot(Zs, OSinddsum, color='b', ls ='-', lw=1)
plt.plot(Zs, finalddsum, color='k', ls='-')
plt.gca().invert_xaxis()
plt.legend(['initial','OSin','final'], loc = 'best', prop={'size':10})
plt.savefig('./' + str(today) + 'thenewddxyzN'+str(haloN)+'.pdf')
plt.clf()

if doplots:
	count = 0
	while (count < TheLen and not fake):
		theN = Ns[count]
		print 'reading particle type 4 for ',theN
		if (theN < 10):
			Nsnapstring = '00'+str(theN)
		elif (theN < 100):
			Nsnapstring = '0'+str(theN)
		else:
			Nsnapstring = str(theN)
		the_snapdir = './hdf5/'
		if (multifile):
			the_snapdir = './snapdir_'+Nsnapstring+'/'
		the_prefix ='snapshot'
		the_suffix = '.hdf5'
		foutname1a = str(today)+'image_stars'+Nsnapstring+'N'+str(haloN)+'.pdf'
		foutname1b = str(today)+'image_DM'+Nsnapstring+'N'+str(haloN)+'.pdf'

		S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
		Rstars = 15.0
		GL.pluscuts_2d_desnity_hist(S, Rstars, finalx[count], finaly[count], finalz[count],  foutname1a, dolabel=True, extrapoints = [[finalx[count], finaly[count], finalz[count]], [x[count], y[count], z[count]], [ShalX[count], ShalY[count], ShalZ[count]]],  thevmax=1e2)
		D = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix)
		GL.pluscuts_2d_desnity_hist(D, Rstars, finalx[count], finaly[count], finalz[count],  foutname1b, dolabel=True, extrapoints = [[finalx[count], finaly[count], finalz[count]], [x[count], y[count], z[count]], [ShalX[count], ShalY[count], ShalZ[count]]],  thevmax=1e2)
		count += 1
		del(S)
		del(D)
