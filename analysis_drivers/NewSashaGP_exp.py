#Analysis routine to read FIRE simulations, calculate outflow rates and other properties of halos in concentric shells

import sys
import os
from datetime import date
import numpy as np
import scipy.stats as stats
from readsnap import readsnap
import Sasha_functions as SF
import particle_tracking_functions as PTF
import graphics_library as GL
from tasz import hubble_param

focusshell = -2.0

today = date.today()
today = '2015-09-05'
final_Rvirfac = 0.02
print today
useprev = False
rhalf_analysis = True
recenter_to_KDE_stars = False
recenter_to_KDE_dark = False
extra_rhalf_analysis = True
use_fixed_halos = 4
#1 for rod, 2 for wod, #3 for sod, #4 for zod

fbar = 0.044 #from paramvalue

multifile = 'n'
do_graphics = 'n'
ContamCheck = True


halo_to_do = [0]
pc = 3.08567758e18
kpc = pc * 1e3
year = 3.15569e7
custombins = True
use_no_Pep = True
Rvirfac = 2.0 #maximum Rvir for which we gather gas particles


use_tempcut = False
use_metalcut = False
use_rhocut = False


use_halo_tree = 'y'
use_physical_bins = 'n'
fixed_OutflowCut = -9999
fixed_inflow_cut = -9999

write_outflow_particles = 'y'
Rvir_must_grow = 'y'

rho_cut_lower = 1e0 #separates 'galactic' gas. 0 for using all gas
rho_cut_upper = 1e99 #separates 'galactic' gas. 0 for using all gas
Temp_cut_lower = 1e0 #seperates 'cool' gas 
Temp_cut_upper = 1e99 #seperates 'cool' gas 
metal_cut_value = 1e-99 #separates 'metal rich' from 'metal poor'
metal_cut_value = 1e99 #separates 'metal rich' from 'metal poor'

Inner_SF_thresh = 0.2 #distinguish between "inner" and "outer" SF

R_bins = 20 # this is the number of bins you'll get from 0 to Rvir. IN reality the outflows are computed starting with bin focus_shell. Default focus_shell=2 (0.2 - 0.3 Rvir)

general_dL = 0.1 #if custom bins are used

extra_physical_bins = True
PhysicalBin_values = [10.0, 25.0, 50.0]
PhysicalBin_width_values = [5.0, 5.0, 5.0]

#histogram properties
vmax_lim = 2000
vmin_lim = -2000
num_bins = 200


if (len(sys.argv) < 2):
	print 'syntax: blah.py Nsnapshot'
	sys.exit()

if (int(sys.argv[1]) < 10):
	Nsnapstring = '00'+str(sys.argv[1])
elif (int(sys.argv[1]) < 100):
	Nsnapstring = '0'+str(sys.argv[1])
else:
	Nsnapstring = str(sys.argv[1])
Nsnap = int(sys.argv[1])

the_snapdir = './'
the_prefix ='snap_convertedshot'
the_suffix = ''

#for multiple snapshots

#for HDF5
the_snapdir = './hdf5/'
if (multifile == 'y'):
	the_snapdir = './snapdir_'+Nsnapstring+'/'
the_prefix ='snapshot'
the_suffix = '.hdf5'



print 'reading particle type 0 (gas?)' 
G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)

print 'reading particle type 1 (fine DM?)'
DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix)

print 'reading particle type 2'
P2 = readsnap(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix)

print 'reading particle type 3'
P3 = readsnap(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix)

print 'reading particle type 4'
S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)

print 'reading particle type 5 (black holes? I hope not)'
P5 = readsnap(the_snapdir, Nsnapstring, 5, snapshot_name=the_prefix, extension=the_suffix,skip_bh=1)

#also lets read in previous snapshot
if (useprev):
	Previous_pos_table, prev_a = PTF.simple_previous_snap_read((Nsnap-1), use_multi=multifile)



#newheader = [npart, massarr, time, redshift, flag_sfr, flag_feedbacktp, npartTotal, flag_cooling, numfiles, boxsize, omega_matter, omega_lambda, hubble, flag_stellarage, flag_metals]

#ensure cosmological parameters are here

thetime = G['header'][2]
redshift = G['header'][3]
boxsize = G['header'][9]
omega_matter = G['header'][10]
omega_L = G['header'][11]
h = G['header'][12]

print [thetime, redshift, boxsize, omega_matter, omega_L, h]



NP0, NP1, NP2, NP3, NP4, NP5 = G['header'][6]
M0, M1, M2, M3, M4, M5 = G['header'][1]

if (M1 > 0):
	print 'adjusting P1 mass'
	DM['m'] = np.zeros(NP1)
	DM['m'][0:NP1] = M1
	#print DM['m']

if (M2 > 0):
	print 'adjusting P2 mass'
	P2['m'] = np.zeros(NP2)
	P2['m'][0:NP2] = M2

if (M3 > 0):
	print 'adjusting P3 mass'
	P3['m'] = np.zeros(NP3)
	P3['m'][0:NP3] = M3

if (M5 > 0):
	print 'adjusting P5 mass'
	P5['m'] = np.zeros(NP5)
	P5['m'][0:NP3] = M5
	
a = float(thetime)
print 'using a = ',a
theredshift = (1.0/a - 1)
redshiftstring = "{0:.3f}".format(theredshift)
print 'z = ',"{0:.3f}".format(theredshift)
if (useprev):
	prev_redshift = (1.0/prev_a - 1)
	prev_redshiftstring = "{0:.3f}".format(prev_redshift)

UnitLength_in_cm=3.085678e21 / h
UnitMass_in_g=1.989e43 / h
UnitVelocity_in_cm_per_s=1.e5
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2)

print 'using the following cosmos ',boxsize, omega_matter, omega_L, h

totalGmass = np.sum(G['m'])
totalDmass = np.sum(DM['m'])
totalSmass = np.sum(S['m'])
if ('m' in P2):
	total2mass = np.sum(P2['m'])
else:
	total2mass = 0
if ('m' in P3):
	total3mass = np.sum(P3['m'])
else:
	total3mass = 0
if ('m' in P5):
	total5mass = np.sum(P5['m'])
else:
	total5mass = 0
totalMass = totalGmass + totalDmass + total2mass +totalSmass + total3mass + total5mass

#take care of density cut
gadget_volume = boxsize * boxsize * boxsize
gas_avg_density = totalMass * (fbar/omega_matter) /gadget_volume

rhocut = (G['rho'] > rho_cut_lower*gas_avg_density) * (G['rho'] < rho_cut_upper*gas_avg_density)

#take care of temperature cut
PhysTemp, PhysRho = SF.convertTemp(G['u'], G['ne'], G['rho'], h)
tempcut = (PhysTemp > Temp_cut_lower)  * (PhysTemp < Temp_cut_upper)



SFR_time_range, age_in_Gyr = SF.cut_stellar_age(a, Nsnap, S['age'], omega_matter, h)
stellar_age_cut =  age_in_Gyr < SFR_time_range/1e9

Metallicity = G['z'][:,0]  #Phil says that the "total metallicity" is more accurate than the sum of individual yields, since the two are discrepant.
Stellar_met = S['z'][:,0]

metalcut = (Metallicity > metal_cut_lower) * (Metallicity < metal_cut_upper)


C = SF.read_halo_catalog(Nsnapstring, redshiftstring)

Shell_Outflow_Line_ar = []
Shell_Inflow_Line_ar = []

Shell_Outflow_LineB_ar = []
Shell_Inflow_LineB_ar = []

Shell_Outflow_Line0_ar = []
Shell_Inflow_Line0_ar = []

Shell_Outflow_LineB0_ar = []
Shell_Inflow_LineB0_ar = []	

Shell_Outflow_LineC_ar = []
Shell_Inflow_LineC_ar = []

Shell_Outflow_LineC0_ar = []
Shell_Inflow_LineC0_ar = []	


rhalf_analysis_ar = []
extra_rhalf_analysis_ar = []


workcount = 0 
Hubble = hubble_param(a, omega_matter, h)

while (workcount < len(halo_to_do)):
	if (use_halo_tree == 'y'):
		halostats = SF.find_halo_now(halo_to_do[workcount], a, therod=use_fixed_halos)
		haloN = halostats[1]
	else:
		haloN = halocount
	
	Rvir = C['Rvir'][haloN]
	Vsig = C['Vsig'][haloN]
	haloX = C['x'][haloN]
	haloY = C['y'][haloN]
	haloZ = C['z'][haloN]
	haloVX = C['vx'][haloN]
	haloVY = C['vy'][haloN]
	haloVZ = C['vz'][haloN]
	M = C['M'][haloN]
	Mstar = C['Mstar'][haloN]
	Mgas = C['Mgas'][haloN]
	Vmax = C['Vmax'][haloN]
	
	if (use_no_Pep and use_halo_tree == 'y'):
		Rvir = halostats[11]
		Vsig = halostats[10]
		haloX = halostats[2]
		haloY = halostats[3]
		haloZ = halostats[4]
		haloVX =halostats[5]
		haloVY = halostats[6]
		haloVZ = halostats[7]
		M = halostats[8]
		Mstar = halostats[13]
		Mgas = halostats[12]
		Vmax = halostats[9]
	
	if (Rvir_must_grow == 'y'):
		print 'instead of ',Rvir,Vsig,M
		Rvir, Vsig, M = SF.check_Rvir_growth(halo_to_do[workcount], a, Rvir, Vsig, M, therod=use_fixed_halos)
		
	print 'using Rvir ', Rvir
	print 'Halo Vmax, Vsig, Vsig/sqrt(3) (1D)',Vmax, Vsig, (Vsig/np.sqrt(3))
	
	Rmax = Rvir*Rvirfac
	
	HaloCat = SF.read_halo_history(halo_to_do[workcount],  rod=use_fixed_halos)
	Final_Rvir = HaloCat['Rvir'][-1]
	
	
	if (extra_physical_bins and custombins):
		biggest_physical_distance = max(PhysicalBin_values) + max(PhysicalBin_width_values)
		biggest_comoving_distance = biggest_physical_distance * (h/a)
		Rmax = max([Rmax, biggest_comoving_distance])
	print 'using Rmax ',Rmax
	


	if (recenter_to_KDE_stars):
		print 'recentering to stars'
		Sdists = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize) 
		Sinner = Sdists < Inner_SF_thresh*Rvir
		S_COM = SF.gaussian_KDE_center(S['p'][:,0][Sinner], S['p'][:,1][Sinner], S['p'][:,2][Sinner], downsample=True)
		haloX = S_COM[0]
		haloY = S_COM[1]
		haloZ = S_COM[2]
	if (recenter_to_KDE_dark):
		print 'recentering to dark'
		DMdists =  SF.calcdist2(haloX, haloY, haloZ, DM['p'][:,0], DM['p'][:,1], DM['p'][:,2], boxsize) 
		Dinner = DMdists <  Inner_SF_thresh*Rvir
		D_COM = SF.gaussian_KDE_center(DM['p'][:,0][Dinner], DM['p'][:,1][Dinner], DM['p'][:,2][Dinner], downsample=True)
		haloX = D_COM[0]
		haloY = D_COM[1]
		haloZ = D_COM[2]

	
	#make cuts to rule out all particles more than Rvir away in each direction.  
	GcloseHx = abs(G['p'][:,0] - haloX) < Rmax
	GcloseHy = abs(G['p'][:,1] - haloY) < Rmax
	GcloseHz = abs(G['p'][:,2] - haloZ) < Rmax
	#combine all cuts
	Gclose = GcloseHx * GcloseHy * GcloseHz  

	#various gas properties cuts
	if (use_rhocut):
		Gclose *=  rhocut 
	if (use_tempcut):
		Gclose *= tempcut
	if (use_metalcut):
		Gclose *=  metalcut		
	
	#now get an actual distance estimate for all remaining gas particles
	dists = SF.calcdist2(haloX, haloY, haloZ, G['p'][:,0][Gclose], G['p'][:,1][Gclose], G['p'][:,2][Gclose], boxsize)
	
	#make another distance cut to get only stuff within a sphere of halo center. And cut everything from the "ISM" within 0.2 Rvir
	Ghalo = dists < Rmax
	#GhaloB = dists > (float(0.0) / float(R_bins)) * Rmax
	#Ghalo *= GhaloB

	#dm particles within the virial radius
	DMcloseHx = abs(DM['p'][:,0] - haloX) < Rmax
	DMcloseHy = abs(DM['p'][:,1] - haloY) < Rmax
	DMcloseHz = abs(DM['p'][:,2] - haloZ) < Rmax
	DMclose = DMcloseHx * DMcloseHy * DMcloseHz
	DMdists = SF.calcdist2(haloX, haloY, haloZ, DM['p'][:,0][DMclose], DM['p'][:,1][DMclose], DM['p'][:,2][DMclose], boxsize)


	if (ContamCheck):
		print 'checking for contamination'
		p2_mass_in_halo = 0
		p3_mass_in_halo = 0
		p5_mass_in_halo = 0
		if (total2mass > 0):
			print 'there are type 2 particles '
			P2closeHx = abs(P2['p'][:,0] - haloX) < Rmax
			P2closeHy = abs(P2['p'][:,1] - haloY) < Rmax
			P2closeHz = abs(P2['p'][:,2] - haloZ) < Rmax
			P2close = P2closeHx * P2closeHy * P2closeHz
			P2dists = SF.calcdist2(haloX, haloY, haloZ, P2['p'][:,0][P2close], P2['p'][:,1][P2close], P2['p'][:,2][P2close], boxsize)
			p2_mass_in_halo = np.sum(P2['m'][P2close])
		if (total3mass > 0):
			print 'there are type 3 particles '
			P3closeHx = abs(P3['p'][:,0] - haloX) < Rmax
			P3closeHy = abs(P3['p'][:,1] - haloY) < Rmax
			P3closeHz = abs(P3['p'][:,2] - haloZ) < Rmax
			P3close = P3closeHx * P3closeHy * P3closeHz
			P3dists = SF.calcdist2(haloX, haloY, haloZ, P3['p'][:,0][P3close], P3['p'][:,1][P3close], P3['p'][:,2][P3close], boxsize)	
			p3_mass_in_halo = np.sum(P3['m'][P3close])
		if (total5mass > 0):
			print 'there are type 5 particles '
			P5closeHx = abs(P5['p'][:,0] - haloX) < Rmax
			P5closeHy = abs(P5['p'][:,1] - haloY) < Rmax
			P5closeHz = abs(P5['p'][:,2] - haloZ) < Rmax
			P5close = P5closeHx * P5closeHy * P5closeHz
			P5dists = SF.calcdist2(haloX, haloY, haloZ, P5['p'][:,0][P5close], P5['p'][:,1][P5close], P5['p'][:,2][P5close], boxsize)	
			p5_mass_in_halo = np.sum(P5['m'][P5close])


	#stars particles within the virial radius
	ScloseHx = abs(S['p'][:,0] - haloX) < Rmax
	ScloseHy = abs(S['p'][:,1] - haloY) < Rmax
	ScloseHz = abs(S['p'][:,2] - haloZ) < Rmax
	
	#mostly i'm interested only in young stars to compute SFR, but also considering all stars with Sclose_comp
	Sclose_young = ScloseHx * ScloseHy * ScloseHz * stellar_age_cut
	Sclose_all = ScloseHx * ScloseHy * ScloseHz
	Sdists_young = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0][Sclose_young], S['p'][:,1][Sclose_young], S['p'][:,2][Sclose_young], boxsize) 
	Sdists_all = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0][Sclose_all], S['p'][:,1][Sclose_all], S['p'][:,2][Sclose_all], boxsize)
	#differentiate between stars formed in "inner" and "outer" parts of the halo right now Inner_SF_thresh has always been 0.2
	Sinner_young = Sdists_young < Inner_SF_thresh*Rvir	 
	Sinner_all = Sdists_all < Inner_SF_thresh*Rvir
	S_OuterA = Sdists_young < 1.0*Rvir
	S_OuterB = Sdists_young > (Inner_SF_thresh*Rvir)
	Souter_young = S_OuterA * S_OuterB
		
	HaloSFR = np.sum(S['m'][Sclose_young][Sinner_young])
	HaloSFR = ((HaloSFR * 1e10 / h) / SFR_time_range)
	HaloSFR_Outer = np.sum(S['m'][Sclose_young][Souter_young])
	HaloSFR_Outer = ((HaloSFR_Outer * 1e10 / h) / SFR_time_range)


	#ditching velocity diagnostics for now.
	#calculate relative position and relative velocity
	#the sqrt(a) is for gadget's weird velocity units
	#the extra factor is to adjust for hubble flow

	GprelX = G['p'][:,0][Gclose][Ghalo] - haloX
	GprelY = G['p'][:,1][Gclose][Ghalo] - haloY
	GprelZ = G['p'][:,2][Gclose][Ghalo] - haloZ	
	GvrelX = G['v'][:,0][Gclose][Ghalo] * np.sqrt(a) - haloVX + GprelX*a*Hubble/(h*1000)
	GvrelY = G['v'][:,1][Gclose][Ghalo] * np.sqrt(a) - haloVY + GprelY*a*Hubble/(h*1000)
	GvrelZ = G['v'][:,2][Gclose][Ghalo] * np.sqrt(a) - haloVZ + GprelZ*a*Hubble/(h*1000)

	
	#calculate some properties of gas
	PhysRho_dense= PhysRho[Gclose][Ghalo] > 1.0 #per cc
	Gmass_galactic_dense = np.sum(G['m'][Gclose][Ghalo][PhysRho_dense]) * 1e10 / h
	
	if (rhalf_analysis):	
		#all dense gas within 2 Rvir right?
		dense_gas_COM = SF.find_COM(G['p'][:,0][Gclose][Ghalo][PhysRho_dense], G['p'][:,1][Gclose][Ghalo][PhysRho_dense], G['p'][:,2][Gclose][Ghalo][PhysRho_dense], G['m'][Gclose][Ghalo][PhysRho_dense])
		dists_dense_COM = SF.calcdist2(dense_gas_COM[0], dense_gas_COM[1], dense_gas_COM[2], G['p'][:,0][Gclose][Ghalo][PhysRho_dense], G['p'][:,1][Gclose][Ghalo][PhysRho_dense], G['p'][:,2][Gclose][Ghalo][PhysRho_dense], boxsize)
		[dense_gas_rhalf, dense_gas_mhalf] = SF.find_rhalf(dists_dense_COM, G['m'][Gclose][Ghalo][PhysRho_dense])
		#alternate method using center of halo
		[dense_gas_rhalf_alt, dense_gas_mhalf_alt] = SF.find_rhalf(dists[Ghalo][PhysRho_dense], G['m'][Gclose][Ghalo][PhysRho_dense])
		
		#extra close distances
		G_extra_close = dists[Ghalo] < Inner_SF_thresh*Rvir
		G_extra_close *= PhysRho_dense
		Xdense_gas_COM = SF.find_COM(G['p'][:,0][Gclose][Ghalo][G_extra_close], G['p'][:,1][Gclose][Ghalo][G_extra_close], G['p'][:,2][Gclose][Ghalo][G_extra_close], G['m'][Gclose][Ghalo][G_extra_close])
		Xdists_dense_COM = SF.calcdist2(Xdense_gas_COM[0], Xdense_gas_COM[1], Xdense_gas_COM[2], G['p'][:,0][Gclose][Ghalo][G_extra_close], G['p'][:,1][Gclose][Ghalo][G_extra_close], G['p'][:,2][Gclose][Ghalo][G_extra_close], boxsize)
		[Xdense_gas_rhalf, Xdense_gas_mhalf] = SF.find_rhalf(Xdists_dense_COM, G['m'][Gclose][Ghalo][G_extra_close])
		#alternate method using center of halo
		[Xdense_gas_rhalf_alt, Xdense_gas_mhalf_alt] = SF.find_rhalf(dists[Ghalo][G_extra_close], G['m'][Gclose][Ghalo][G_extra_close])

		#new stars 
		Sinner_youngCOM = SF.find_COM(S['p'][:,0][Sclose_young][Sinner_young], S['p'][:,1][Sclose_young][Sinner_young], S['p'][:,2][Sclose_young][Sinner_young], S['m'][Sclose_young][Sinner_young])
		dists_Sinner_youngCOM = SF.calcdist2(Sinner_youngCOM[0], Sinner_youngCOM[1], Sinner_youngCOM[2], S['p'][:,0][Sclose_young][Sinner_young], S['p'][:,1][Sclose_young][Sinner_young], S['p'][:,2][Sclose_young][Sinner_young], boxsize)
		[Sinner_rhalf, Sinner_mhalf] = SF.find_rhalf(dists_Sinner_youngCOM, S['m'][Sclose_young][Sinner_young])

		#alternative method using center of halo 
		[Sinner_rhalf_alt, Sinner_mhalf_alt] = SF.find_rhalf(Sdists_young[Sinner_young], S['m'][Sclose_young][Sinner_young])
		
		#all stars
		Sinner_allCOM = SF.find_COM(S['p'][:,0][Sclose_all][Sinner_all], S['p'][:,1][Sclose_all][Sinner_all], S['p'][:,2][Sclose_all][Sinner_all], S['m'][Sclose_all][Sinner_all])
		dists_Sinner_allCOM = SF.calcdist2(Sinner_allCOM[0], Sinner_allCOM[1], Sinner_allCOM[2], S['p'][:,0][Sclose_all][Sinner_all], S['p'][:,1][Sclose_all][Sinner_all], S['p'][:,2][Sclose_all][Sinner_all], boxsize)
		[Sinner_rhalf_all, Sinner_mhalf_all] = SF.find_rhalf(dists_Sinner_allCOM, S['m'][Sclose_all][Sinner_all])

		#alternative method using center of halo
		[Sinner_rhalf_all_alt, Sinner_mhalf_all_alt] = SF.find_rhalf(Sdists_all[Sinner_all], S['m'][Sclose_all][Sinner_all])
		
		rhalf_line = [Nsnap, theredshift, SFR_time_range,  -1, haloX, haloY, haloZ, 0, dense_gas_COM[0], dense_gas_COM[1], dense_gas_COM[2], dense_gas_rhalf, dense_gas_rhalf_alt, dense_gas_mhalf, dense_gas_mhalf_alt, 2,  Xdense_gas_COM[0], Xdense_gas_COM[1], Xdense_gas_COM[2],Xdense_gas_rhalf, Xdense_gas_rhalf_alt, Xdense_gas_mhalf, Xdense_gas_mhalf_alt, 3, Sinner_youngCOM[0], Sinner_youngCOM[1], Sinner_youngCOM[2], Sinner_rhalf, Sinner_rhalf_alt, Sinner_mhalf, Sinner_mhalf_alt, 4, Sinner_allCOM[0], Sinner_allCOM[1], Sinner_allCOM[2], Sinner_rhalf_all, Sinner_rhalf_all_alt, Sinner_mhalf_all, Sinner_mhalf_all_alt]
		
		linestring = ''
		for q in rhalf_line:
			linestring += "{0:.6g}".format(q) + '  '
		linestring+=' \n'
		
		rhalf_analysis_ar.append(linestring)

		if (extra_rhalf_analysis):
			dists_Sall_youngCOM = SF.calcdist2(Sinner_youngCOM[0], Sinner_youngCOM[1], Sinner_youngCOM[2], S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize)
			dists_DMall_youngCOM = SF.calcdist2(Sinner_youngCOM[0], Sinner_youngCOM[1], Sinner_youngCOM[2], DM['p'][:,0], DM['p'][:,1], DM['p'][:,2], boxsize)
			dists_Gall_youngCOM = SF.calcdist2(Sinner_youngCOM[0], Sinner_youngCOM[1], Sinner_youngCOM[2], G['p'][:,0], G['p'][:,1], G['p'][:,2], boxsize)
			
			DM_within_rhalf_cut = dists_DMall_youngCOM < Sinner_rhalf
			G_within_rhalf_cut = dists_Gall_youngCOM < Sinner_rhalf
			S_within_rhalf_cut = dists_Sall_youngCOM < Sinner_rhalf	
			TotalM_within_rhalf = np.sum(DM['m'][DM_within_rhalf_cut]) + np.sum(G['m'][G_within_rhalf_cut]) + np.sum(S['m'][S_within_rhalf_cut])
			
			DM_within_2rhalf_cut = dists_DMall_youngCOM < 2*Sinner_rhalf
			G_within_2rhalf_cut = dists_Gall_youngCOM < 2*Sinner_rhalf
			S_within_2rhalf_cut = dists_Sall_youngCOM < 2*Sinner_rhalf
			TotalM_within_2rhalf = np.sum(DM['m'][DM_within_2rhalf_cut]) + np.sum(G['m'][G_within_2rhalf_cut]) + np.sum(S['m'][S_within_2rhalf_cut])
			
			DM_within_1kpc_cut = dists_DMall_youngCOM < 1.0 *h / a
			G_within_1kpc_cut = dists_Gall_youngCOM < 1.0 *h / a
			S_within_1kpc_cut = dists_Sall_youngCOM < 1.0 *h / a
			TotalM_within_1kpc = np.sum(DM['m'][DM_within_1kpc_cut]) + np.sum(G['m'][G_within_1kpc_cut]) + np.sum(S['m'][S_within_1kpc_cut])
			dist_of_youngCOM = SF.calcdist2(haloX, haloY, haloZ, [Sinner_youngCOM[0]],  [Sinner_youngCOM[1]],[Sinner_youngCOM[2]], boxsize)[0]
			DM_within_youngCOM_cut =  DMdists < dist_of_youngCOM
			G_within_youngCOM_cut =   dists < dist_of_youngCOM
			S_within_youngCOM_cut = Sdists_all < dist_of_youngCOM
			
			M_within_young_COM = np.sum(DM['m'][DM_within_youngCOM_cut]) + np.sum(G['m'][G_within_youngCOM_cut]) + np.sum(S['m'][S_within_youngCOM_cut])
			
			extra_rhalf_line = [Nsnap, theredshift, Sinner_rhalf, dist_of_youngCOM, TotalM_within_rhalf, TotalM_within_2rhalf, TotalM_within_1kpc, M_within_young_COM]
			
			linestring = ''
			for q in extra_rhalf_line:
				linestring += "{0:.6g}".format(q) + '  '
			linestring+=' \n'
		
			extra_rhalf_analysis_ar.append(linestring)			
	
	#rescale positions so that magnitude = 1 (unit vectors)
	GprelX_RS = GprelX / dists[Ghalo]
	GprelY_RS = GprelY / dists[Ghalo]
	GprelZ_RS = GprelZ / dists[Ghalo]
	
	#find dot product of position unit vector and relative velocity
	v_rad = GprelX_RS*GvrelX  + GprelY_RS*GvrelY + GprelZ_RS*GvrelZ


	#find particles that will meet outflow criteria
	if (fixed_OutflowCut > -1 ):
		OutflowCut = v_rad > fixed_OutflowCut
		justfast = abs(v_rad) > fixed_OutflowCut  
	else:
		OutflowCut = (v_rad > Vsig/np.sqrt(3)) 
		justfast = abs(v_rad) >  Vsig/np.sqrt(3)  
	if (fixed_inflow_cut > -1 ):
		InflowCut = (v_rad < -1*fixed_inflow_cut) 
	else:
		InflowCut = (v_rad < -1*Vsig/np.sqrt(3)) 

	OutflowCut0 =  v_rad > 0
	InflowCut0 = v_rad < 0

	if (custombins):
		R_bin_array_min = np.array([0.0, 0.45, 0.95, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9])
		
		R_bin_array_max = R_bin_array_min+0.1
		
		#R_bin_array_max[1]  = 0.02
		R_bin_array_max[0]  = (Final_Rvir * final_Rvirfac)/Rvir
		
		if (extra_physical_bins):
			PhysicalBin_values = np.array(PhysicalBin_values)
			PhysicalBin_width_values = np.array(PhysicalBin_width_values)				
			PhysicalBin_meds = PhysicalBin_values * h / (a*Rvir)
			PhysicalBin_widths = PhysicalBin_width_values * h / (a*Rvir)
			PhysicalBin_maxesds = PhysicalBin_meds + PhysicalBin_widths
			PhysicalBin_mins = PhysicalBin_meds - PhysicalBin_widths
			R_bin_array_min = np.append(PhysicalBin_mins, R_bin_array_min)
			R_bin_array_max = np.append(PhysicalBin_maxesds, R_bin_array_max)
			less_than_zero = R_bin_array_min<0.0
			R_bin_array_min[less_than_zero] = 0.0
			
		R_bins = len(R_bin_array_min)		
		R_bin_array_med = (R_bin_array_min + R_bin_array_max)/2.0
		print 'doing custom bins ', R_bin_array_min
		
		
		ShellCuts, CrossCuts, ShellThickness_ar = SF.shell_check_custom(dists[Ghalo], Rvir, R_bin_array_min , R_bin_array_max, a, h)
		print 'widths are this:',ShellThickness_ar

	else:
		ShellCuts, ShellList, ShellThickness = SF.shell_check(dists[Ghalo], (Rmax), R_bins, a, h)

	#find which new outflows originated in ISM
	if (useprev):
		print 'doing work for ',haloN
		match_previous_particles = np.in1d( Previous_pos_table[:,0], G['id'][Gclose][Ghalo])
		
		prevID = Previous_pos_table[:,0][match_previous_particles]
		prevx = Previous_pos_table[:,1][match_previous_particles]
		prevy = Previous_pos_table[:,2][match_previous_particles]
		prevz = Previous_pos_table[:,3][match_previous_particles]
		pH = SF.find_halo_now(halo_to_do[workcount], prev_a, therod=use_fixed_halos)
		
		prev_sort_order = np.argsort(prevID)
		
		prevhaloX = pH[2]
		prevhaloY = pH[3]
		prevhaloZ = pH[4]
		
		previous_dists = SF.calcdist2(prevhaloX, prevhaloY, prevhaloZ, prevx, prevy, prevz, boxsize)
		
		
		bincount = 0
		
		noswept_shell_cut = []
		noswept_shell_in = []
		#noswept_shell_cut.append([])
		if (custombins):
			while (bincount < R_bins):
				AllOuter = CrossCuts[bincount]
				AllInner = np.invert(CrossCuts[bincount])
				previous_close = (previous_dists[prev_sort_order] <  R_bin_array_med[bincount]*Rvir)
				previous_far = (previous_dists[prev_sort_order] >R_bin_array_med[bincount]*Rvir)
				sortorder = np.argsort(G['id'][Gclose][Ghalo])
				AllOuter[sortorder] *= previous_close
				AllInner[sortorder] *= previous_far
				noswept_shell_cut.append(AllOuter)
				noswept_shell_in.append(AllInner)
				bincount+=1
		else:
			while (bincount < R_bins):
				AllOuter = np.copy(ShellCuts[bincount])
				AllInner = np.copy(ShellCuts[bincount])
				subbincount = bincount
				while (subbincount < R_bins):
					AllOuter += ShellCuts[subbincount]
					subbincount+=1
				subbincount = bincount
				while (subbincount >= 0):
					AllInner += ShellCuts[subbincount]
					subbincount -= 1
				previous_close = (previous_dists[prev_sort_order] <  ((float(bincount) / float(R_bins)) * Rmax))
				previous_far = (previous_dists[prev_sort_order] >  ((float(bincount+1) / float(R_bins)) * Rmax))
				#previous_far = np.invert(previous_close)
				#has to be 
			#need to have compatible sortorder
				sortorder = np.argsort(G['id'][Gclose][Ghalo])
			
				#print 'out of ', len(dists[Ghalo][OutflowCut*BinCut]),' outflows'
				#BinCuts *= previous_close  - no this isnt right because i want to allow particles that are in other shells - further shells to be counted as part of the outflow rate  - BUT I do need a check to make sure the outflow is NOW at least past the shell 
				AllOuter[sortorder] *= previous_close
				AllInner[sortorder] *= previous_far
				noswept_shell_cut.append(AllOuter)
				noswept_shell_in.append(AllInner)
				bincount+=1
		
	


	 ######## OUTPUT 2: main calculation: outflows per shell
	Shell_Outflow_Line = []
	Shell_Inflow_Line = []
	
	Shell_Outflow_LineB = []
	Shell_Inflow_LineB = []
	
	Shell_Outflow_Line0 = []
	Shell_Inflow_Line0 = []
	
	Shell_Outflow_LineB0 = []
	Shell_Inflow_LineB0 = []	

	Shell_Outflow_LineC = []
	Shell_Inflow_LineC = []
	Shell_Outflow_LineC0 = []
	Shell_Inflow_LineC0 = []	
	


	
	
	if (ContamCheck):
		ContamLine = []
	
	SumOfOutflows = 0
	SumOfInflows = 0

	
	bincount = 0
	while (bincount < R_bins):
		if (custombins):
			ShellThickness = ShellThickness_ar[bincount]
		BinCuts = ShellCuts[bincount]
		if (useprev):
			ThreshCut = noswept_shell_cut[bincount] 
			InThreshCut = noswept_shell_in[bincount]
		else:
			ThreshCut = BinCuts
			InThreshCut = BinCuts
		#print 'a ',BinCuts[0:50]
		#print 'b' ,ThreshCut[0:50]
		NParticles_in_shell = len(dists[Ghalo][BinCuts])
		TotalShellMass = np.sum(G['m'][Gclose][Ghalo][BinCuts])
		TotalShellMetal = np.sum(G['m'][Gclose][Ghalo][BinCuts]*Metallicity[Gclose][Ghalo][BinCuts])
		NParticles_in_outflow = len(dists[Ghalo][OutflowCut * BinCuts])
		NParticles_in_outflow0 = len(dists[Ghalo][OutflowCut0 * BinCuts])
		NParticles_in_inflow = len(dists[Ghalo][InflowCut * BinCuts])
	#rate of mass transfer through the shell in Msun/yr? -totally obsolete atm
		NParticles_in_inflow0 = len(dists[Ghalo][InflowCut0 * BinCuts])
		ratetotal = ((G['m'][Gclose][Ghalo][BinCuts]*1e10/h) * (v_rad[BinCuts]*1e5*year/kpc)) / ShellThickness
		#rate of mass flow that exceeds the velocity cut
		outflowtotal = ((G['m'][Gclose][Ghalo][OutflowCut * BinCuts]*1e10/h) * (v_rad[OutflowCut*BinCuts]*1e5*year/kpc)) / ShellThickness
		outflowtotal0 = ((G['m'][Gclose][Ghalo][OutflowCut0 * BinCuts]*1e10/h) * (v_rad[OutflowCut0*BinCuts]*1e5*year/kpc)) / ShellThickness
		inflowtotal = ((G['m'][Gclose][Ghalo][InflowCut * BinCuts]*1e10/h) * (v_rad[InflowCut*BinCuts]*1e5*year/kpc)) / ShellThickness
		inflowtotal0 = ((G['m'][Gclose][Ghalo][InflowCut0 * BinCuts]*1e10/h) * (v_rad[InflowCut0*BinCuts]*1e5*year/kpc)) / ShellThickness
		metalinflowtotal = ((G['m'][Gclose][Ghalo][InflowCut * BinCuts]*1e10/h) *  Metallicity[Gclose][Ghalo][InflowCut * BinCuts]* (v_rad[InflowCut*BinCuts]*1e5*year/kpc)) / ShellThickness
		metalinflowtotal0 = ((G['m'][Gclose][Ghalo][InflowCut0 * BinCuts]*1e10/h) *  Metallicity[Gclose][Ghalo][InflowCut0 * BinCuts]* (v_rad[InflowCut0*BinCuts]*1e5*year/kpc)) / ShellThickness

		Complete_Outflow_Total = ((G['m'][Gclose][Ghalo][OutflowCut * ThreshCut]*1e10/h) * (v_rad[OutflowCut*ThreshCut]*1e5*year/kpc)) / ShellThickness
		Complete_Outflow_Npart = len((G['m'][Gclose][Ghalo][OutflowCut * ThreshCut]))		
		Complete_Outflow_Mass = np.sum(G['m'][Gclose][Ghalo][OutflowCut * ThreshCut])
		Complete_Outflow_Metal = np.sum(G['m'][Gclose][Ghalo][OutflowCut * ThreshCut] * Metallicity[Gclose][Ghalo][OutflowCut * ThreshCut])
		metaloutflowtotal = ((G['m'][Gclose][Ghalo][OutflowCut * BinCuts]*1e10/h) *  Metallicity[Gclose][Ghalo][OutflowCut * BinCuts] * (v_rad[OutflowCut*BinCuts]*1e5*year/kpc)) / ShellThickness

		Complete_Outflow_Total0 = ((G['m'][Gclose][Ghalo][OutflowCut0 * ThreshCut]*1e10/h) * (v_rad[OutflowCut0*ThreshCut]*1e5*year/kpc)) / ShellThickness
		Complete_Outflow_Npart0 = len((G['m'][Gclose][Ghalo][OutflowCut0 * ThreshCut]))		
		Complete_Outflow_Mass0 = np.sum(G['m'][Gclose][Ghalo][OutflowCut0 * ThreshCut])
		Complete_Outflow_Metal0 = np.sum(G['m'][Gclose][Ghalo][OutflowCut0 * ThreshCut] * Metallicity[Gclose][Ghalo][OutflowCut0 * ThreshCut])
		metaloutflowtotal0 = ((G['m'][Gclose][Ghalo][OutflowCut0 * BinCuts]*1e10/h) *  Metallicity[Gclose][Ghalo][OutflowCut0 * BinCuts] * (v_rad[OutflowCut0*BinCuts]*1e5*year/kpc)) / ShellThickness

		
		
		
		Complete_Inflow_Total = ((G['m'][Gclose][Ghalo][InflowCut * InThreshCut]*1e10/h) * (v_rad[InflowCut*InThreshCut]*1e5*year/kpc)) / ShellThickness
		Complete_Inflow_Npart = len((G['m'][Gclose][Ghalo][InflowCut * InThreshCut]))		
		Complete_Inflow_Mass = np.sum(G['m'][Gclose][Ghalo][InflowCut * InThreshCut])		
		Complete_Inflow_Metal = np.sum(G['m'][Gclose][Ghalo][InflowCut * InThreshCut] * Metallicity[Gclose][Ghalo][InflowCut * InThreshCut])	

		Complete_Inflow_Total0 = ((G['m'][Gclose][Ghalo][InflowCut0 * InThreshCut]*1e10/h) * (v_rad[InflowCut0*InThreshCut]*1e5*year/kpc)) / ShellThickness
		Complete_Inflow_Npart0 = len((G['m'][Gclose][Ghalo][InflowCut0 * InThreshCut]))		
		Complete_Inflow_Mass0 = np.sum(G['m'][Gclose][Ghalo][InflowCut0 * InThreshCut])		
		Complete_Inflow_Metal0 = np.sum(G['m'][Gclose][Ghalo][InflowCut0 * InThreshCut] * Metallicity[Gclose][Ghalo][InflowCut0 * InThreshCut])	


		Vmag = G['v'][:,0]*G['v'][:,0] + G['v'][:,1]*G['v'][:,1] + G['v'][:,2]*G['v'][:,2]
		Vmag = np.sqrt(Vmag)
		Vmag *= np.sqrt(a)

		dOKinE = (G['m'][Gclose][Ghalo][OutflowCut * BinCuts]*1e10/h) * Vmag[Gclose][Ghalo][OutflowCut * BinCuts] * Vmag[Gclose][Ghalo][OutflowCut * BinCuts] * 0.5
		dOKinE0 = (G['m'][Gclose][Ghalo][OutflowCut0 * BinCuts]*1e10/h) * Vmag[Gclose][Ghalo][OutflowCut0 * BinCuts] * Vmag[Gclose][Ghalo][OutflowCut0 * BinCuts] * 0.5
		dIKinE = (G['m'][Gclose][Ghalo][InflowCut * BinCuts]*1e10/h) * Vmag[Gclose][Ghalo][InflowCut * BinCuts] * Vmag[Gclose][Ghalo][InflowCut * BinCuts] * 0.5		
		dIKinE0 = (G['m'][Gclose][Ghalo][InflowCut0 * BinCuts]*1e10/h) * Vmag[Gclose][Ghalo][InflowCut0 * BinCuts] * Vmag[Gclose][Ghalo][InflowCut0 * BinCuts] * 0.5		
		dKinETot = (G['m'][Gclose][Ghalo][BinCuts]*1e10/h) * Vmag[Gclose][Ghalo][BinCuts] * Vmag[Gclose][Ghalo][BinCuts] * 0.5		
		dOthreshKinE = (G['m'][Gclose][Ghalo][OutflowCut * ThreshCut]*1e10/h) * Vmag[Gclose][Ghalo][OutflowCut * ThreshCut] * Vmag[Gclose][Ghalo][OutflowCut * ThreshCut] * 0.5
		dOthreshKinE0 = (G['m'][Gclose][Ghalo][OutflowCut0 * ThreshCut]*1e10/h) * Vmag[Gclose][Ghalo][OutflowCut0 * ThreshCut] * Vmag[Gclose][Ghalo][OutflowCut0 * ThreshCut] * 0.5
		dIthreshKinE = (G['m'][Gclose][Ghalo][InflowCut * InThreshCut]*1e10/h) * Vmag[Gclose][Ghalo][InflowCut * InThreshCut] * Vmag[Gclose][Ghalo][InflowCut * InThreshCut] * 0.5
		dIthreshKinE0 = (G['m'][Gclose][Ghalo][InflowCut0 * InThreshCut]*1e10/h) * Vmag[Gclose][Ghalo][InflowCut0 * InThreshCut] * Vmag[Gclose][Ghalo][InflowCut0 * InThreshCut] * 0.5
		InflowTemp0 = PhysTemp[Gclose][Ghalo][InflowCut0 * BinCuts]
		InflowTemp = PhysTemp[Gclose][Ghalo][InflowCut * BinCuts]
		OutflowTemp0 = PhysTemp[Gclose][Ghalo][OutflowCut0 * BinCuts]
		OutflowTemp = PhysTemp[Gclose][Ghalo][OutflowCut * BinCuts]
		
		# should be multiplied *2e33*1e5*1e5 for cgs units
		# i guess i will do that, and also divide by 1e51. 
		# hence, 2e-8 
		#units: 1e51 erg		

		OKinE = np.sum(dOKinE*2e-8)
		OKinE0 = np.sum(dOKinE0*2e-8)
		IKinE = np.sum(dIKinE*2e-8)
		IKinE0 = np.sum(dIKinE0*2e-8)
		KinETot = np.sum(dKinETot*2e-8)
		OthreshKinE = np.sum(dOthreshKinE*2e-8)
		OthreshKinE0 = np.sum(dOthreshKinE0*2e-8)
		IthreshKinE = np.sum(dIthreshKinE*2e-8)
		IthreshKinE0 = np.sum(dIthreshKinE0*2e-8)
		
		AvgInflowTemp = 10**np.mean(np.log10(InflowTemp))
		AvgOutflowTemp = 10**np.mean(np.log10(OutflowTemp))
		AvgInflowTemp0 = 10**np.mean(np.log10(InflowTemp0))
		AvgOutflowTemp0 =  10**np.mean(np.log10(OutflowTemp0))
		

		
		OIntE = np.mean(PhysTemp[Gclose][Ghalo][OutflowCut * BinCuts])
		OIntE0 = np.mean(PhysTemp[Gclose][Ghalo][OutflowCut0 * BinCuts]) 
		IIntE = np.mean(PhysTemp[Gclose][Ghalo][InflowCut * BinCuts])
		IIntE0 = np.mean(PhysTemp[Gclose][Ghalo][InflowCut0 * BinCuts])
		OthreshIntE = np.mean(PhysTemp[Gclose][Ghalo][OutflowCut * ThreshCut])
		OthreshIntE0 = np.mean(PhysTemp[Gclose][Ghalo][OutflowCut0 * ThreshCut])
		IthreshIntE = np.mean(PhysTemp[Gclose][Ghalo][InflowCut * InThreshCut])
		IthreshIntE0 = np.mean(PhysTemp[Gclose][Ghalo][InflowCut0 * InThreshCut])
		IntETot =  np.mean(PhysTemp[Gclose][Ghalo][BinCuts])
		
		
		
		#add up all particles
		sumoutflow = np.sum(outflowtotal)
		suminflow = np.sum(inflowtotal)
		sumoutflow0 = np.sum(outflowtotal0)
		suminflow0 = np.sum(inflowtotal0)
		
		sumComplete = np.sum(Complete_Outflow_Total)
		#add up a grand total for all shells
		#this is a useless quantity now. 

		SumOfOutflows +=sumoutflow
		SumOfInflows +=suminflow

		summetaloutflow = np.sum(metaloutflowtotal)
		summetalinflow = np.sum(metalinflowtotal)
		summetaloutflow0 = np.sum(metaloutflowtotal0)
		summetalinflow0 = np.sum(metalinflowtotal0)
		
		
		#tack on the relevant values to the line
		Shell_Outflow_Line.extend([bincount, sumoutflow, NParticles_in_shell, NParticles_in_outflow, Complete_Outflow_Npart, Complete_Outflow_Mass ])
		Shell_Outflow_Line0.extend([bincount, sumoutflow0, NParticles_in_shell, NParticles_in_outflow0, Complete_Outflow_Npart0, Complete_Outflow_Mass0])

	
		Shell_Outflow_LineB.extend([bincount, Complete_Outflow_Metal, summetaloutflow, TotalShellMetal, TotalShellMass])
		Shell_Outflow_LineB0.extend([bincount, Complete_Outflow_Metal0, summetaloutflow0, TotalShellMetal, TotalShellMass])
		
		Shell_Outflow_LineC0.extend([bincount, OKinE0, OIntE0, OthreshKinE0, OthreshIntE0, KinETot, IntETot,AvgOutflowTemp0])
		Shell_Outflow_LineC.extend([bincount, OKinE, OIntE, OthreshKinE, OthreshIntE, KinETot, IntETot, AvgOutflowTemp])

		#used to be : Shell_Outflow_Line.extend([bincount, sumoutflow, NParticles_in_shell, NParticles_in_outflow, NSecondTimers, SecondTimerRate])

		Shell_Inflow_Line.extend([bincount, suminflow, NParticles_in_shell, NParticles_in_inflow, Complete_Inflow_Npart, Complete_Inflow_Mass])
		Shell_Inflow_Line0.extend([bincount, suminflow0, NParticles_in_shell, NParticles_in_inflow0, Complete_Inflow_Npart0, Complete_Inflow_Mass0])


		Shell_Inflow_LineB.extend([bincount, Complete_Inflow_Metal, summetalinflow]) 
		Shell_Inflow_LineB0.extend([bincount, Complete_Inflow_Metal0, summetalinflow0])
		Shell_Inflow_LineC.extend([bincount, IKinE, IIntE, IthreshKinE, IthreshIntE, AvgInflowTemp])
		Shell_Inflow_LineC0.extend([bincount, IKinE0, IIntE0, IthreshKinE0, IthreshIntE0, AvgInflowTemp0])
		
		
		if (ContamCheck and custombins):
			#print 'computing contam'
			P2mass_in_bin = 0
			P2N_in_bin = 0
			P3mass_in_bin = 0
			P3N_in_bin = 0
			P5mass_in_bin = 0
			P5N_in_bin = 0
			DMmass_in_bin = 0
			DMN_in_bin = 0		
			if (p2_mass_in_halo > 0):
				P2bincut = (P2dists/Rvir > R_bin_array_min[bincount]) * (P2dists/Rvir <  R_bin_array_max[bincount])
				P2mass_in_bin = np.sum(P2['m'][P2close][P2bincut])
				P2N_in_bin = len(P2['m'][P2close][P2bincut])
			if (p3_mass_in_halo > 0):
				P3bincut = (P3dists/Rvir > R_bin_array_min[bincount]) * (P3dists/Rvir <  R_bin_array_max[bincount])
				P3mass_in_bin = np.sum(P3['m'][P3close][P3bincut])
				P3N_in_bin = len(P3['m'][P3close][P3bincut])				
			if (p5_mass_in_halo > 0):
				P5bincut = (P5dists/Rvir > R_bin_array_min[bincount]) * (P5dists/Rvir <  R_bin_array_max[bincount])
				P5mass_in_bin = np.sum(P5['m'][P5close][P5bincut])
				P5N_in_bin = len(P5['m'][P5close][P5bincut])							
			DMbincut = (DMdists/Rvir > R_bin_array_min[bincount]) * (DMdists/Rvir <  R_bin_array_max[bincount])
			DMmass_in_bin = np.sum(DM['m'][DMclose][DMbincut])
			DMN_in_bin = len(DM['m'][DMclose][DMbincut])
			
			fhiresbin = DMmass_in_bin / (P5mass_in_bin + P3mass_in_bin + P2mass_in_bin + DMmass_in_bin)   

			Gbincut = (dists/Rvir > R_bin_array_min[bincount]) * (dists/Rvir <  R_bin_array_max[bincount])
			Gmass_in_bin = np.sum(G['m'][Gclose][Gbincut])
			GN_in_bin = len(G['m'][Gclose][Gbincut])

			Sbincut = (Sdists_all/Rvir > R_bin_array_min[bincount]) * (Sdists_all/Rvir <  R_bin_array_max[bincount])
			Smass_in_bin = np.sum(S['m'][Sclose_all][Sbincut])
			SN_in_bin = len(S['m'][Sclose_all][Sbincut])
			Smetal_in_bin = np.sum(S['m'][Sclose_all][Sbincut]*Stellar_met[Sclose_all][Sbincut])

			ContamLine.extend( [bincount, DMN_in_bin, DMmass_in_bin, P2N_in_bin, P2mass_in_bin, P3N_in_bin, P3mass_in_bin, P5N_in_bin, P5mass_in_bin, fhiresbin, Gmass_in_bin,GN_in_bin, Smass_in_bin,SN_in_bin, Smetal_in_bin])

		bincount += 1


	foutname = str(today)+'Outflow'+'N'+str(halo_to_do[workcount])+'.txt'
	#g = open(foutname, 'a')
	
	foutname2 = str(today)+'B_Outflow'+'N'+str(halo_to_do[workcount])+'.txt'
	#g2 = open(foutname2, 'a')
	
	#here's where i put in the quantities i want
	line = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, SFR_time_range]
	lineB = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, SFR_time_range]
	line0 = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, SFR_time_range]
	lineB0 = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, SFR_time_range]
	lineC = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, SFR_time_range]
	lineC0 = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, SFR_time_range]

	#used to be : line = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfOutflows, NumNewParticles, NewParticleFlux,  NewParticleMass, SFR_time_range, NewParticleMassInner, NewParticleMass_plusswept]

	#add in all outflows for each shell
	line.extend(Shell_Outflow_Line)
	lineB.extend(Shell_Outflow_LineB)
	line0.extend(Shell_Outflow_Line0)
	lineB0.extend(Shell_Outflow_LineB0)
	lineC.extend(Shell_Outflow_LineC)
	lineC0.extend(Shell_Outflow_LineC0)
	
	
	#annoying procedure to convert list to string
	linestring = ''
	linestringB = ''
	linestring0 = ''
	linestringB0 = ''
	linestringC = ''
	linestringC0 = ''

	for q in line:
		linestring += "{0:.6g}".format(q) + '  '
	linestring+=' \n'
	
	for q in lineB:
		linestringB += "{0:.6g}".format(q) + '  '
	linestringB+=' \n'
	
	for q in line0:
		linestring0 += "{0:.6g}".format(q) + '  '
	linestring0+=' \n'
	
	for q in lineB0:
		linestringB0 += "{0:.6g}".format(q) + '  '
	linestringB0+=' \n'

	for q in lineC0:
		linestringC0 += "{0:.6g}".format(q) + '  '
	linestringC0+=' \n'

	for q in lineC:
		linestringC += "{0:.6g}".format(q) + '  '
	linestringC+=' \n'


	Shell_Outflow_Line_ar.append(linestring)
	Shell_Outflow_LineB_ar.append(linestringB)
	Shell_Outflow_Line0_ar.append(linestring0)
	Shell_Outflow_LineB0_ar.append(linestringB0)
	Shell_Outflow_LineC_ar.append(linestringC)
	Shell_Outflow_LineC0_ar.append(linestringC0)

	#g.write(linestring)
	#g.close()

	#g2.write(linestringB)
	#g2.close()

	
	linestring = ''
	linestringB = ''
	linestring0 = ''
	linestringB0 = ''
	linestringC = ''
	linestringC0 = ''

	#do the same for inflows
	foutname = str(today)+'Inflow'+'N'+str(halo_to_do[workcount])+'.txt'
	foutname2 = str(today)+'B_Inflow'+'N'+str(halo_to_do[workcount])+'.txt'

	#g = open(foutname, 'a')
	#g2 = open(foutname2, 'a')

	line = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfInflows]
	lineB = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfInflows]
	line0 = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfInflows]
	lineB0 = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfInflows]
	lineC = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfInflows]
	lineC0 = [Nsnap, theredshift, Rvir, HaloSFR, HaloSFR_Outer, Gmass_galactic_dense, M, Mgas, Mstar, totalSmass, SumOfInflows]


	line0.extend(Shell_Inflow_Line0)
	lineB0.extend(Shell_Inflow_LineB0)

	line.extend(Shell_Inflow_Line)
	lineB.extend(Shell_Inflow_LineB)

	lineC.extend(Shell_Inflow_LineC)
	lineC0.extend(Shell_Inflow_LineC0)
	
	
	for q in line:
		linestring += "{0:.6g}".format(q) + '  '
	linestring +=' \n'
	
	for q in lineB:
		linestringB += "{0:.6g}".format(q) + '  '
	linestringB +=' \n'
	
	for q in line0:
		linestring0 += "{0:.6g}".format(q) + '  '
	linestring0 +=' \n'
	
	for q in lineB0:
		linestringB0 += "{0:.6g}".format(q) + '  '
	linestringB0 +=' \n'

	for q in lineC:
		linestringC += "{0:.6g}".format(q) + '  '
	linestringC +=' \n'

	for q in lineC0:
		linestringC0 += "{0:.6g}".format(q) + '  '
	linestringC0 +=' \n'


	Shell_Inflow_Line_ar.append(linestring)
	Shell_Inflow_LineB_ar.append(linestringB)
	Shell_Inflow_Line0_ar.append(linestring0)
	Shell_Inflow_LineB0_ar.append(linestringB0)
	Shell_Inflow_LineC_ar.append(linestringC)
	Shell_Inflow_LineC0_ar.append(linestringC0)
	
	if (ContamCheck):
		print 'writing contam string'
		contamstring = ''
		contamintro = [Nsnap, theredshift, Rvir, M, Mgas, Mstar]
		contamintro.extend(ContamLine)
		foutname = str(today)+'Contam'+'N'+str(halo_to_do[workcount])+'.txt'

		for q in contamintro:
			contamstring += "{0:.6g}".format(q) + '  '
		contamstring+=' \n'
		g3 = open(foutname, 'a')
		g3.write(contamstring)
		g3.close()
		
	
	######OUTPUT 3 velocity histogram & velocity stats 
	# a little redundant from OUTPUT 2, but i think its worth it to keep it separate
	#while (bincount < R_bins):
	#if (len(v_rad[OutflowCut])>1):
	#	vmax_lim = max(200, max(dotprods))
	#	vmin_lim = min(-200, min(dotprods))
	#else:
	bincount = 0
	histresults = []
	masshistresults = []
	velocity_medians = [] 
	velocity_maxima = []
	
	#first do things for all material - no shells
	
	
	ratetotal = ((G['m'][Gclose][Ghalo]*1e10/h) * (v_rad*1e5*year/kpc)) / ShellThickness
	histresults_mass = np.histogram(v_rad, num_bins, range=(vmin_lim, vmax_lim), weights = G['m'][Gclose][Ghalo])
	#print 'hist results length for ',bincount, len(histresults_mass[0])
	maxflowrate = SF.sasha_max(histresults_mass[0])
	bulkindex = np.where(histresults_mass[0] == maxflowrate)[0][0]
	bulkflowvelocity = histresults_mass[1][bulkindex]
	outflow_median = stats.nanmedian(v_rad[OutflowCut])
	inflow_median = stats.nanmedian(v_rad[InflowCut])
	general_median = stats.nanmedian(v_rad)
	maxoutflowvelocity = SF.sasha_max(v_rad[OutflowCut]) 
	maxinflowvelocity = SF.sasha_min(v_rad[InflowCut])
	if (outflow_median > 0):
		ind90 = 0.9 * len(v_rad[OutflowCut])
		velocity90 = np.sort(v_rad[OutflowCut])[ind90]
	else: 
		velocity90 = 0.0
	
	
	#print 'medians for -1',outflow_median, general_median, inflow_median, maxoutflowvelocity
	velocity_table = []
	velocity_stats_line = [-1, 	general_median, bulkflowvelocity, outflow_median, inflow_median, maxoutflowvelocity, maxinflowvelocity, velocity90]
	#print 'starting, appending line of length ',len(velocity_stats_line)
	velocity_table.append(velocity_stats_line)
	
	
	while (bincount < R_bins):
	
		if (custombins):
			ShellThickness = ShellThickness_ar[bincount]
		#first get the histogram[s]
		BinCuts = ShellCuts[bincount]
		if (useprev):
			ThreshCut = noswept_shell_cut[bincount] 
			InThreshCut = noswept_shell_in[bincount]
		else:
			ThreshCut = BinCuts
			InThreshCut = BinCuts
		ratetotal = ((G['m'][Gclose][Ghalo][ThreshCut]*1e10/h) * (v_rad[ThreshCut]*1e5*year/kpc)) / ShellThickness
		shell_histresults = np.histogram(v_rad[ThreshCut], num_bins, range=(vmin_lim, vmax_lim), weights = ratetotal)
		shell_histresults_mass = np.histogram(v_rad[ThreshCut], num_bins, range=(vmin_lim, vmax_lim), weights = G['m'][Gclose][Ghalo][ThreshCut])
		#print 'hist results length for',bincount, len(shell_histresults[0])
		histresults.append(shell_histresults[0])
		#print 'hist results length for',bincount, len(shell_histresults_mass[0])

		masshistresults.append(shell_histresults_mass[0])
		
		#now nitty gritty		
		maxflowrate = SF.sasha_max(shell_histresults_mass[0])
		bulkindex = np.where(shell_histresults_mass[0] == maxflowrate)[0][0]
		bulkflowvelocity = shell_histresults_mass[1][bulkindex]
		outflow_median = stats.nanmedian(v_rad[ThreshCut*OutflowCut])
		maxoutflowvelocity = SF.sasha_max(v_rad[OutflowCut*ThreshCut])
		inflow_median = stats.nanmedian(v_rad[InThreshCut*InflowCut])
		general_median = stats.nanmedian(v_rad[BinCuts])
		#print 'medians for',bincount, inflow_median, general_median, outflow_median, maxoutflowvelocity
		maxinflowvelocity = SF.sasha_min(v_rad[InflowCut*InThreshCut])
		if (outflow_median > 0):
			ind90 = 0.9 * len(v_rad[OutflowCut*ThreshCut])
			velocity90 = np.sort(v_rad[OutflowCut*ThreshCut])[ind90]
		else: 
			velocity90 = 0.0
		velocity_stats_line = [bincount, 	general_median, bulkflowvelocity, outflow_median, inflow_median, maxoutflowvelocity, maxinflowvelocity, velocity90]
		#print 'appending line of length ',len(velocity_stats_line)
		velocity_table.append(velocity_stats_line)
		
		bincount+=1



	histlim1 = shell_histresults[1][:-1]
	histlim2 = shell_histresults[1][1:]
	
	#print 'ok hist results ',len(histlim1)
	
	histresults.append(histlim1)
	histresults.append(histlim2)
	masshistresults.append(histlim1)
	masshistresults.append(histlim2)

	histresultsT = (np.array(histresults)).T
	masshistresultsT = (np.array(masshistresults)).T
	
	vhist_header = [(Vsig/np.sqrt(3.0)), a, Nsnap, num_bins, vmin_lim, vmax_lim, M]
	vhist_header_string = SF.line_to_string(vhist_header, newline='n')
	vhist_foutname = str(today)+'histogram'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.txt'
	vhist_foutname2 = str(today)+'histogram_mass'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.txt'
	#print histresultsT.shape
	np.savetxt(vhist_foutname, histresultsT,  fmt='%1.6e', header=vhist_header_string)
	np.savetxt(vhist_foutname2, masshistresultsT,  fmt='%1.6e', header=vhist_header_string)
	velocity_table = np.array(velocity_table)
	vstats_foutname = str(today)+'v_stats'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.txt'
	np.savetxt(vstats_foutname, velocity_table,  fmt='%1.6e', header=vhist_header_string)


	#OUTPUT 4 simple density diagrams
	
	if (do_graphics == 'y'):
		print 'creating images '

		foutname1a = str(today)+'image_gasoutflows'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname1b = str(today)+'image_gasinflows'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname1c = str(today)+'image_gasstagnant'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		
		foutname11a = str(today)+'image_Wgasoutflows'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname11b = str(today)+'image_Wgasinflows'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname11c = str(today)+'image_Wgasstagnant'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		
		
		#foutname7 = str(today)+'phase_diagram_inner'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname7a = str(today)+'phase_diagram_total'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname7b = str(today)+'phase_diagram_outflow'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname7c = str(today)+'phase_diagram_0.25_total'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname7d = str(today)+'phase_diagram_0.25_outflow'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname7e = str(today)+'phase_diagram_0.95_total'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname7f = str(today)+'phase_diagram_0.95_outflow'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'		
		
		
		foutname8a = str(today)+'phase_histogram_total'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname8b = str(today)+'phase_histogram_outflow'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname8c = str(today)+'phase_histogram_0.25_total'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname8d = str(today)+'phase_histogram_0.25_outflow'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname8e = str(today)+'phase_histogram_0.95_total'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutname8f = str(today)+'phase_histogram_0.95_outflow'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'		

		#first doing cuts for the density histograms
		sendCuts = np.copy(Gclose) 
		sendCuts[sendCuts] *= Ghalo
		
		sendCutsB = np.copy(sendCuts)
		sendCutsC = np.copy(sendCuts)
		
		sendCuts[sendCuts] *= OutflowCut
		sendCutsB[sendCutsB]*= InflowCut
		sendCutsC[sendCutsC] *= np.invert(OutflowCut) * np.invert(InflowCut)
		

		#cut these figures at rvir
		#GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1a, Cuts=sendCuts, focusshell=focusshell)
		#GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1b, Cuts=sendCutsB, focusshell=focusshell)
		#GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1c, Cuts=sendCutsC, focusshell=focusshell)


		ThePhysWeights = G['m']*1e10/h
		binuse = 300
		
		GunitLength = 2 * (Rvir/h) * a / binuse
		GunitArea = GunitLength * GunitLength
		ThePhysWeights /= GunitArea
		meangas = np.mean(ThePhysWeights)


		myvmin = 1 * meangas
		myvmax = 1e3 * meangas

		myvmin = 2e6 / 4
		myvmax = 2e9 / 4
		
		
		GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname11a, Cuts=sendCuts, focusshell=focusshell, thevmin = myvmin, thevmax = myvmax, PhysWeight=ThePhysWeights, dolabel=True)
		GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname11b, Cuts=sendCutsB, focusshell=focusshell, thevmin = myvmin, thevmax = myvmax, PhysWeight=ThePhysWeights, dolabel=True)
		GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname11c, Cuts=sendCutsC, focusshell=focusshell, thevmin = myvmin, thevmax = myvmax, PhysWeight=ThePhysWeights, dolabel=True)
		
		
		sendCuts = np.copy(Gclose) 
		sendCuts[sendCuts] *= Ghalo
		
		#dont even do the general case
		#phase_diagram(PhysTemp, PhysRho, G['m'], foutname7a, sendCuts)
		
		sendCuts2 = np.copy(sendCuts)
		sendCuts3 = np.copy(sendCuts)
		sendCuts4 = np.copy(sendCuts)
		sendCuts5 = np.copy(sendCuts)
		sendCuts6 = np.copy(sendCuts)
		
		#lets do one for all gas within Rvir, one for outflows within Rvir, one for all gas at 0.25 Rvir, one for outflows at 0.25 Rvir, one for all gas at 0.95 Rvir, one for outflows at 0.95 Rvir
		
		sendCuts[sendCuts] *= ( ShellCuts[0] + ShellCuts[2] + ShellCuts[1] + ShellCuts[4] + ShellCuts[3] + ShellCuts[5] + ShellCuts[6] + ShellCuts[7] + ShellCuts[8] + ShellCuts[9])
		sendCuts2[sendCuts2] *= OutflowCut * (ShellCuts[2] + ShellCuts[1] + ShellCuts[4] + ShellCuts[3] + ShellCuts[5] + ShellCuts[6] + ShellCuts[7] + ShellCuts[8] + ShellCuts[9])
		
		# for plotting outer stuff
		#sendCuts3[sendCuts3] *= np.invert(ShellCuts[0])
		
		sendCuts3[sendCuts3] *= ShellCuts[2]
		
		#for plotting inner outflows
		#sendCuts2[sendCuts2]*=OutflowCut*ShellCuts[0]
		sendCuts4[sendCuts4]*= OutflowCut * ShellCuts[2]
		
		sendCuts5[sendCuts5] *= ShellCuts[9]
		
		#for plotting inner outflows
		#sendCuts2[sendCuts2]*=OutflowCut*ShellCuts[0]
		sendCuts6[sendCuts6]*= OutflowCut * ShellCuts[9]
		
		#GL.phase_diagram(PhysTemp, PhysRho, G['m'], foutname7a, sendCuts, zstring=redshiftstring)
		#GL.phase_diagram(PhysTemp, PhysRho, G['m'], foutname7b, sendCuts2, zstring=redshiftstring)
		#GL.phase_diagram(PhysTemp, PhysRho, G['m'], foutname7c, sendCuts3, zstring=redshiftstring)
		#GL.phase_diagram(PhysTemp, PhysRho, G['m'], foutname7d, sendCuts4, zstring=redshiftstring)
		#GL.phase_diagram(PhysTemp, PhysRho, G['m'], foutname7e, sendCuts5, zstring=redshiftstring)
		#GL.phase_diagram(PhysTemp, PhysRho, G['m'], foutname7f, sendCuts6, zstring=redshiftstring)


		#GL.phase_histogram_delux(PhysTemp, sendCuts, sendCuts2, foutname8a, zstring=redshiftstring)
		#GL.phase_histogram_delux(PhysTemp, sendCuts3, sendCuts4, foutname8c, zstring=redshiftstring)
		#GL.phase_histogram_delux(PhysTemp, sendCuts5, sendCuts6, foutname8e, zstring=redshiftstring)

		
		
	print 'done with halo', halo_to_do[workcount]
	workcount+=1
	
print 'writing outflows '
count = 0
while (count < len(halo_to_do)):

	foutname = str(today)+'Outflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname2 = str(today)+'B_Outflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname3 = str(today)+'C_Outflow'+'N'+str(halo_to_do[count])+'.txt'
	g = open(foutname, 'a')
	g2 = open(foutname2, 'a')
	g3 = open(foutname3, 'a')
	g.write(Shell_Outflow_Line_ar[count])
	g2.write(Shell_Outflow_LineB_ar[count]) 
	g3.write(Shell_Outflow_LineC_ar[count])
	g.close()
	g2.close()
	g3.close()

	foutname = str(today)+'Inflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname2 = str(today)+'B_Inflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname3 = str(today)+'C_Inflow'+'N'+str(halo_to_do[count])+'.txt'
	g = open(foutname, 'a')
	g2 = open(foutname2, 'a')
	g3 = open(foutname3, 'a')
	g.write(Shell_Inflow_Line_ar[count])
	g2.write(Shell_Inflow_LineB_ar[count]) 
	g3.write(Shell_Inflow_LineC_ar[count]) 
	g.close()
	g2.close()
	g3.close()
	
	foutname = str(today)+'0_Inflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname2 = str(today)+'0B_Inflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname3 = str(today)+'0C_Inflow'+'N'+str(halo_to_do[count])+'.txt'
	g = open(foutname, 'a')
	g2 = open(foutname2, 'a')
	g3 = open(foutname3, 'a')
	g.write(Shell_Inflow_Line0_ar[count])
	g2.write(Shell_Inflow_LineB0_ar[count]) 
	g3.write(Shell_Inflow_LineC0_ar[count]) 
	g.close()
	g2.close()
	g3.close()

	foutname = str(today)+'0_Outflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname2 = str(today)+'0B_Outflow'+'N'+str(halo_to_do[count])+'.txt'
	foutname3 = str(today)+'0C_Outflow'+'N'+str(halo_to_do[count])+'.txt'
	g = open(foutname, 'a')
	g2 = open(foutname2, 'a')
	g3 = open(foutname3, 'a')
	g.write(Shell_Outflow_Line0_ar[count])
	g2.write(Shell_Outflow_LineB0_ar[count]) 
	g3.write(Shell_Outflow_LineC0_ar[count]) 
	g.close()
	g2.close()
	g3.close()
	
	if (rhalf_analysis):
		foutname = str(today)+'rhalf'+'N'+str(halo_to_do[count])+'.txt'
		foutname2 = str(today)+'rhalf'+'N'+str(halo_to_do[count])+'.txt'
		g = open(foutname, 'a')
		g.write(rhalf_analysis_ar[count])
		g.close()
		if (extra_rhalf_analysis):
			foutname = str(today)+'extrarhalf'+'N'+str(halo_to_do[count])+'.txt'
			foutname2 = str(today)+'extrarhalf'+'N'+str(halo_to_do[count])+'.txt'
			g = open(foutname, 'a')
			g.write(extra_rhalf_analysis_ar[count])
			g.close()			
	
	count += 1




print 'all done'
sys.exit()

