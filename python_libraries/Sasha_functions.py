#General Functions for analysis

import numpy as np
from tasz import tfora
import sys
from scipy.interpolate import interp1d
import scipy.stats as stats

def convertTemp(Tb, Neb, rho, h):
	#converts temperature and density to physical units
	H_MASSFRAC = 0.76
	PROTONMASS  = 1.6726e-24
	GAMMA = 5.0/3.0
	BOLTZMANN  = 1.3806e-16
	
	UnitLength_in_cm=3.085678e21 / h
	UnitMass_in_g=1.989e43 / h
	UnitVelocity_in_cm_per_s=1.e5
	UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
	UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2)

	MeanWeights = 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb) * PROTONMASS
	MeanWeights_con = 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb) 

	converted_rho = rho * UnitMass_in_g * 6.02e23 / MeanWeights_con
	converted_rho /= (UnitLength_in_cm * UnitLength_in_cm *UnitLength_in_cm)
	
	TrueTemp = Tb *  MeanWeights/BOLTZMANN * (GAMMA-1.) * UnitEnergy_in_cgs/ UnitMass_in_g
	return (TrueTemp, converted_rho)

def cut_stellar_age(a, N, age, omeganot, littleh):
	#returns the time since the previous snapshot, as well as an array that gives the actual age of all stellar particles in Myr (as opposed to formation time, as given by code)
	finname = 'output_times.txt'
	f = open(finname)
	SdarsA = np.loadtxt(f)
	prev_epoch = SdarsA[max((N-1), 0)]
	SFR_time_range = (tfora(a, omeganot, littleh) - tfora(prev_epoch, omeganot, littleh)) * 1e9
	print 'using exact SFR', a, prev_epoch, (SFR_time_range)/1e7, ' * 1e7 yrs'
	f.close()
	#print age
	time_in_Gyr = tfora(age, omeganot, littleh)
	age_in_Gyr = tfora(a, omeganot, littleh) - time_in_Gyr
	print ' i interpolated your ages'
	#print age_in_Gyr
	#print ' the times of formation are: ',time_in_Gyr
	return (SFR_time_range, age_in_Gyr)


def read_halo_catalog(Nsnapstring, redshiftstring, usembp='n'):
	finname = 'Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstring+'.AHF_halos' #or whatever prefix
	f = open(finname)
	C = np.loadtxt(f)
	halID = C[:,0]
	x = C[:,5]
	y = C[:,6]
	z = C[:,7]
	vx = C[:,8]
	vy = C[:,9]
	vz = C[:,10]
	M = C[:,3]
	Vmax = C[:,16]
	Vsig = C[:,18]
	Rvir = C[:,11]
	Mgas = C[:,53]
	Mstar = C[:,73]
	#if (use_mbp=='y'): not even going to bother writing this for now - it isnt very accurate for halo center in Amiga. May want to use COM instead though
		
	return {'ID':halID, 'x':x, 'y':y, 'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'M':M, 'Vmax':Vmax, 'Vsig':Vsig, 'Rvir':Rvir, 'Mgas':Mgas, 'Mstar':Mstar}

def read_halo_history( halocount, rod=0):
	halo_number_string = str(halocount)
	if (halocount < 10):
		halo_number_string='0'+halo_number_string
	if (rod==1):
		halo_in_name = 'halos/halo_000' + halo_number_string + '_rod.dat'
	elif (rod==2):
		halo_in_name = 'halos/halo_000' + halo_number_string + '_wod.dat'
	elif (rod==3):
		halo_in_name = 'halos/halo_000' + halo_number_string + '_sod.dat'
	elif (rod==4):
		halo_in_name = 'halos/halo_000' + halo_number_string + '_zod.dat'
 	elif (rod==5):
 		halo_in_name = 'halos/halo_000' + halo_number_string + '_ffzod.dat'
	else:
		halo_in_name = 'halos/halo_000' + halo_number_string + '.dat'
	hfile = open(halo_in_name)
	H = np.loadtxt(hfile)
	halozs = H[:,0]
	haloIDs = H[:,1]
	host = H[:,2]
	x = H[:,6]
	y = H[:,7]
	z = H[:,8]
	vx = H[:,9]
	vy = H[:,10]
	vz = H[:,11]
	M = H[:,4]
	Vmax = H[:,17]
	Vsig = H[:,19]
	Rvir = H[:,12]
	Mgas = H[:,54]
	Mstar = H[:,74]
	
	npart = H[:,5]
	ngas= H[:,53]
	nstar = H[:,73]
	fhires = H[:,38]

	return {'redshift':halozs, 'ID':haloIDs, 'x':x, 'y':y, 'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'M':M, 'Vmax':Vmax, 'Vsig':Vsig, 'Rvir':Rvir, 'Mgas':Mgas, 'Mstar':Mstar, 'ngas':ngas, 'nstar':nstar, 'fhires':fhires, 'host':host, 'npart':npart}

def find_halo_now(halocount, a, therod=0):
	theredshift = 1.0 / a - 1.0
	#if (theredshift >= 4.501):
	#	print 'this code only works for z >= 4.5'
	#	print 'for now i am adjusting your redshift to 4.5.'
	#	theredshift = 4.50
	#	a = 1.0 / (1.0 + theredshift)
	H = read_halo_history(int(halocount), rod=therod)
	halozs = H['redshift']
	haloIDs = H['ID']
	epoch_trace_count = 0
	epoch = 0 
	the_haloN = -9999
	for haloredshift in halozs:
		if ( abs(haloredshift - theredshift)<0.0001):
			the_haloN = int(haloIDs[epoch_trace_count])
			epoch = epoch_trace_count
		epoch_trace_count += 1
	if (the_haloN == -9999):
		print 'oops there wasnt a match in the assembly history. Quitting'
		sys.exit()
	else:
		print 'found it! ',the_haloN
		HaloLine = [halozs[epoch], haloIDs[epoch], H['x'][epoch], H['y'][epoch], H['z'][epoch], H['vx'][epoch], H['vy'][epoch], H['vz'][epoch], H['M'][epoch], H['Vmax'][epoch], H['Vsig'][epoch], H['Rvir'][epoch], H['Mgas'][epoch], H['Mstar'][epoch]]
	return HaloLine
	
def shell_check(dists, Rend, R_bins, a, h, use_physical_bins='n', physical_extent=50):

	TheCuts = []
	bincount = 0
	ShellList = np.zeros(len(dists)) - 1 
	while (bincount < R_bins):
		DistMin = Rend * bincount/R_bins
		DistMax = Rend * (bincount+1)/R_bins
	
		#by default, the virial radius is cut into 10 bins. Here i'm also allowing for a fixed physical scale.
		if (use_physical_bins == 'y'):
			print 'using physical bins ',physical_extent, 'kpc/h'
			#gotta divide by a to go to comoving coordinates
			Distmin = bincount/(Rbins) * physical_extent / a
			Distmax = (bincount+1)/(Rbins) * physical_extent / a 

		ShellThickness = (DistMax - DistMin)*a/h

		Cut1 = dists > DistMin
		Cut2 = dists < DistMax
		Cuts = Cut1 * Cut2
		ShellList[Cuts] = bincount
		TheCuts.append(Cuts)
		bincount+=1
	TheCuts = np.array(TheCuts)
	return (TheCuts, ShellList, ShellThickness)

def shell_check_custom(dists, Rend, R_bin_array_min, R_bin_array_max, a, h):

	TheCuts = []
	CrossCuts = []
	ShellThickness_ar = []
	bincount = 0
	#R_bin_array should be something like [0.0, 0.1, 0.2, 0.3, etc... 1.0... 2.0?]  length = Numbins+1
	
	while (bincount < len(R_bin_array_min)):
		DistMin = Rend * R_bin_array_min[bincount]
		DistMax = Rend * R_bin_array_max[bincount]
		
		CrossMin = (DistMin + DistMax)/2.0
	
		#by default, the virial radius is cut into 10 bins. Here i'm also allowing for a fixed physical scale.
		
		ShellThickness_ar.append((DistMax - DistMin)*a/h)

		Cut1 = dists > DistMin
		Cut2 = dists < DistMax
		Cuts = Cut1 * Cut2		
		TheCuts.append(Cuts)
		
		#changing this to be a cut for stuff that's crossed
		CCut = dists > CrossMin
		CrossCuts.append(CCut)
		
		bincount+=1
	TheCuts = np.array(TheCuts)
	ShellThickness_ar = np.array(ShellThickness_ar)
	CrossCuts = np.array(CrossCuts)
	return (TheCuts, CrossCuts, ShellThickness_ar)

def calcdist2(x1, y1, z1, posx, posy, posz, boxsize):
#calculates distances of an array of points (particles) from a single given point (halo center)
#for your reference: x1, y1, z1 are floats for halo center or hwatever, posx, posy posz are arrays for positions of particles or whatever
	 xdist = abs(posx - x1)
	 ydist = abs(posy - y1)
	 zdist = abs(posz - z1)
	 
	 #adjust for periodic boundary conditions
	 x_straggler = xdist/boxsize > 0.5
	 if (len(xdist[x_straggler])) > 0:
	 	#print 'adjusting x stragglers ', len(xdist[x_straggler])
	 	xdist[x_straggler] = boxsize - xdist[x_straggler]
	 y_straggler = ydist/boxsize > 0.5
	 if (len(ydist[y_straggler])) > 0:
	 	ydist[y_straggler] = boxsize - ydist[y_straggler]	 	
	 	#print 'adjusting y stragglers ', len(ydist[y_straggler])
	 z_straggler = zdist/boxsize > 0.5
	 if (len(zdist[z_straggler])) > 0:
	 	zdist[z_straggler] = boxsize - zdist[z_straggler]
	 	#print 'adjusting z stragglers ', len(zdist[z_straggler])

	 
	 dists = pow(((xdist*xdist + ydist*ydist + zdist*zdist)), 0.5)
 	 return dists
 	 
def line_to_string(line, newline = 'y', tabs = 'n'):
	linestring = ''
	for q in line:
		if (type(q) != str):
			linestring += "{0:.6g}".format(q)  
		else:
			linestring += q 
		if (tabs == 'y'):
			linestring += '\t'
		else:
			linestring += '  '
	if (newline == 'y'):
		linestring+=' \n'
	return linestring
	
def array_to_string(line):
	a = []
	for q in line:
		a.append(str(q))
	return a
	
	
def sasha_max(a):
	if (len(a) > 0):
		return np.max(a)
	else:
		print 'warning there was an empty array max'
		return -99999999

def sasha_min(a):
	if (len(a) > 0):
		return np.min(a)
	else:
		print 'warning there was an empty array min'
		return 99999999

def check_Rvir_growth(halocount, a, Rvir, Vsig, Mass, therod=0):
	#check to make sure that the snapshot you are using is not affected by a glitch in the halo catalog that causes a temporary dip in Rvir, Mvir, Vsig 
	H = read_halo_history(int(halocount), rod=therod)
	history_Rvirs = H['Rvir']
	history_M = H['M']
	history_Vsig = H['Vsig']
	history_Zs = H['redshift']
	history_As = 1.0 / (1.0 + history_Zs)
	history_Rvirs_phys = history_Rvirs * history_As
	RvirPhys = Rvir * a
	cut = history_As < a
	historical_max_phys = sasha_max(history_Rvirs_phys[cut])  
	#historical_max_cmv = sasha_max(history_Rvirs[cut])
	if (RvirPhys >= historical_max_phys):
		print 'Rvir is already the biggest'
		return (Rvir, Vsig, Mass)
	maxindex = np.where(history_Rvirs_phys[cut]==historical_max_phys)[0][0]
	#maxindex_cmv = np.where(history_Rvirs[cut]==historical_max_cmv)[0][0]	
	print 'adjusting back to z = ',history_Zs[cut][maxindex]
	print 'Rvir is now ',history_Rvirs_phys[cut][maxindex]/a
	print 'instead of meager ',Rvir
	return ((history_Rvirs_phys[cut][maxindex]/a), history_Vsig[cut][maxindex], history_M[cut][maxindex])

def check_Rvir_growth_extra(halocount, a, Rvir, Vsig, Mass, Mstar, Mgas, therod=0):
	#check to make sure that the snapshot you are using is not affected by a glitch in the halo catalog that causes a temporary dip in Rvir, Mvir, Vsig 
	H = read_halo_history(int(halocount), rod=therod)
	history_Rvirs = H['Rvir']
	history_M = H['M']
	history_Vsig = H['Vsig']
	history_Mstar = H['Mstar']
	history_Mgas = H['Mgas']
	theZ = 1.0 / (1.0 + theZ)

	history_Zs = H['redshift']
	history_As = 1.0 / (1.0 + history_Zs)
	history_Rvirs_phys = history_Rvirs * history_As
	RvirPhys = Rvir * a
	cut = history_As < a
	historical_max_phys = sasha_max(history_Rvirs_phys[cut])  
	#historical_max_cmv = sasha_max(history_Rvirs[cut])
	if (RvirPhys >= historical_max_phys):
		print 'Rvir is already the biggest'
		return (Rvir, Vsig, Mass, Mstar, Mgas)
	maxindex = np.where(history_Rvirs_phys[cut]==historical_max_phys)[0][0]
	#maxindex_cmv = np.where(history_Rvirs[cut]==historical_max_cmv)[0][0]	
	print 'adjusting back to z = ',history_Zs[cut][maxindex]
	deltaZ = history_Zs[cut][maxindex] - theZ 
	print 'Rvir is now ',history_Rvirs_phys[cut][maxindex]/a
	print 'instead of meager ',Rvir
	
	the_new_Rvir = history_Rvirs_phys[cut][maxindex]/a
	the_new_Vsig = history_Vsig[cut][maxindex]
	the_new_M = Mass
	the_new_Mstar = Mstar
	the_new_Mgas = Mgas

	if (Mass < 1.01 * history_M[cut][maxindex]):
		the_new_M = history_M[cut][maxindex]
		print 'also adjusting Mass ', 		the_new_M/Mass
		if (Mstar < 1.5 * history_M[cut][maxindex] and 1.5*RvirPhys < historical_max_phys ):
			the_new_Mstar = history_M[cut][maxindex]
			print 'also adjusting Mstar ', 		the_new_Mstar/Mstar

			if (Mgas < 2.0 * history_Mgas[cut][maxindex]):
				the_new_Mgas = history_Mgas[cut][maxindex]
				print 'full monty, also adjusting Mgas ', 		the_new_Mgas/Mgas

	
	return (the_new_Rvir , the_new_Vsig, the_new_M, the_new_Mstar, the_new_Mgas)

def N_for_z(z):
	finname = 'output_times.txt'
	l =  np.loadtxt(finname)
	the_zs = 1.0 / l  - 1.0
	count = 0 
	for haloredshift in the_zs:
		if ( abs(haloredshift - z)<0.0001):
			#print 'z ',z,' corresponds to N ',count
			return count	
		count +=1
	return -1

def z_for_N(N):
	finname = 'output_times.txt'
	l =  np.loadtxt(finname)
	the_zs = 1.0 / l  - 1.0
	if (N < len(the_zs)):
		#print 'N ',N,' corresponds to z ',the_zs[N]
		return the_zs[N]

def find_COM(x, y, z, mass):
	#assumes you don't have to deal with periodic boundary...
	totalmass = np.sum(mass)
	totalXmass = np.sum(x*mass)
	totalYmass = np.sum(y*mass)
	totalZmass = np.sum(z*mass)
	
	COMx = totalXmass / totalmass
	COMy = totalYmass / totalmass
	COMz = totalZmass / totalmass
		
	return [COMx, COMy, COMz]

def find_rhalf(dists, mass):
	totalmass = np.sum(mass)
	distorder = np.argsort(dists)
	halfmass = totalmass / 2.0
	count = 0
	thecumsum = 0
	halfMR = 0
	for pmass in mass[distorder]:
		thecumsum+=pmass
		if (thecumsum > halfmass):
			halfMR = dists[distorder][count]
			break
		count += 1
	return [halfMR, halfmass]

def gaussian_KDE_center(x, y, z, downsample=False, downsample_lim=1e4):

	xfornow = x
	yfornow = y
	zfornow = z
	if (len(x)<1):
		return [0,0,0]
	if (downsample):
		if (len(x) > downsample_lim):		
			print 'down sampling from ',len(x),' to ',downsample_lim
			f = np.random.choice(np.arange(len(x)), int(downsample_lim), replace=False)
			xfornow = x[f]
			yfornow = y[f]
			zfornow = z[f]
	xyz = np.vstack([xfornow,yfornow,zfornow])
	kde = stats.gaussian_kde(xyz)
	density = kde(xyz)
	f = xyz.T[np.argmax(density)]
	return f
	
def parametric_bin_center(M, x, y, z, mins, maxes, nbins):
	#print 'a'
	dx = (maxes[0] - mins[0])/float(nbins)
	dy = (maxes[1] - mins[1])/float(nbins)
	dz = (maxes[2] - mins[2])/float(nbins)
	
	x0 = mins[0]
	y0 = mins[1]
	z0 = mins[2]
	
	countx = 0
	county = 0
	countz = 0

	maxpos = [x0,y0,z0]
	maxM = 0
	#print 'a'
	while (countx < nbins):
		thex = x0 + float(countx)*dx
		xnext = thex + dx
		xcut = (x<xnext)*(x>thex)
		county = 0 
		while (county < nbins):
			they = y0 + float(county)*dy
			ynext = they + dy
			ycut = (y<ynext)*(y>they)
			countz = 0
			while (countz < nbins):
				thez = z0 + float(countz)*dz
				znext = thez + dz
				zcut = (z<znext)*(z>thez)
				
				M_in_cut = len(M[xcut*ycut*zcut])
				if (M_in_cut > maxM):
					maxM = M_in_cut			
					maxpos = [thex+(dx/2.0), they+(dy/2.0), thez + (dz/2.0)]
				#print countz, maxpos
				countz+=1
			#print county, maxpos
			county+=1
		countx+=1
		#print countx, thex, maxpos
	return maxpos
	
def read_profile(Nsnapstring, redshiftstring, haloN):
	finname = 'Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstring+'.AHF_profiles' #or whatever prefix
	f = open(finname)
	dars = np.loadtxt(f)
	dtable = []
	curtable = []
	curtable2 =  []
	count =0
	prev = -1 
	for line in dars:
		r = line[0]
		if (count > 0 and r<0 and prev>0):
			dtable.append(curtable2)
			curtable = []
			curtable2 = []
		if (r > 0):
			if (len(curtable2) > 1):
				curtable2 = np.vstack([curtable2,line])
			else: 
				curtable2 = np.array(line)
			curtable.append(line)
		prev = r
		count +=1 
	dtable.append(curtable)
	dtable = np.array(dtable)
	return dtable[haloN]		 
	

