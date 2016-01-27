#meant to get info about the halo in shells

#follow the same general pattern of NewSashaGP_exp.py
#but don't bother with outflow rate calculations; only get average properties of shells. 

import sys
import numpy as np
import Sasha_functions as SF
import outflow_rates as OR
from tasz import hubble_param
from readsnap import readsnap
from datetime import date
today = str(date.today())


year = 3.15569e7
pc = 3.08567758e18
kpc = pc * 1e3

halo_to_do = 0
multifile = True
use_fixed_halos = 0
readDM = False
use_empirical_yield = True

InnerDivide = 0.1 #In units of Rvir
UsePhysicalDivide = False
PhysicalDivideISM = 10.0  #in units of physical kpc

ColdTempBound = 10**4.0
LowIonUpperBound = 10**4.7
OVIUpperBound = 10**5.3

use_OuterBoundary = False 
OuterBoundary_comoving = 100 

if (len(sys.argv) < 2):
	print 'syntax: blah.py Nsnapshot'
	sys.exit()

#setting name of snapshot
if (int(sys.argv[1]) < 10):
	Nsnapstring = '00'+str(sys.argv[1])
elif (int(sys.argv[1]) < 100):
	Nsnapstring = '0'+str(sys.argv[1])
else:
	Nsnapstring = str(sys.argv[1])
Nsnap = int(sys.argv[1])

the_snapdir = './hdf5/'
if (multifile):
	the_snapdir = './snapdir_'+Nsnapstring+'/'
the_prefix ='snapshot'
the_suffix = '.hdf5'

#read the snapshot
print 'reading particle type 0 (gas?)' 
G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)

print 'reading particle type 4'
S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)

if (readDM):
	print 'reading particle type 5 (black holes? I hope not)'
	P5 = readsnap(the_snapdir, Nsnapstring, 5, snapshot_name=the_prefix, extension=the_suffix,skip_bh=1)

	print 'reading particle type 1 (fine DM?)'
	DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix)

	print 'reading particle type 2'
	P2 = readsnap(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix)

	print 'reading particle type 3'
	P3 = readsnap(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix)


Metallicity = G['z'][:,0]  #Phil says that the "total metallicity" is more accurate than the sum of individual yields, since the two are discrepant.
Stellar_met = S['z'][:,0]
Gas_oxygen = G['z'][:,4]
Stellar_oxygen = S['z'][:,4]
Total_metal_mass =  (np.sum(Stellar_met*S['m']) + np.sum(Metallicity*G['m']))

OxyYield = (np.sum(Stellar_oxygen*S['m']) + np.sum(Gas_oxygen*G['m'])) / Total_metal_mass

EmpiricalYield = Total_metal_mass / ((np.sum(S['m']) + np.sum(G['m'])))

a = G['header'][2]
redshift = G['header'][3]
boxsize = G['header'][9]
omega_matter = G['header'][10]
omega_L = G['header'][11]
h = G['header'][12]
redshiftstring = "{0:.3f}".format(redshift)
Hubble = hubble_param(a, omega_matter, h)
PhysTemp, PhysRho = SF.convertTemp(G['u'], G['ne'], G['rho'], h)

halostats = SF.find_halo_now(halo_to_do, a, therod=use_fixed_halos)

print 'halo stats', halostats
haloN = halostats[1]

#read halo catalog, assign halo properties
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


Rvir, Vsig, M = SF.check_Rvir_growth(halo_to_do, a, Rvir, Vsig, M, therod=use_fixed_halos)

if (use_OuterBoundary):
	Rvir = OuterBoundary_comoving

#setting the locations of each spherical shell, as a fraction of Rvir
R_bin_array_min = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
R_bin_array_max = R_bin_array_min+0.1

Gclose = OR.prepare_X_close(G, Rvir, haloX, haloY, haloZ)	
Gdists = SF.calcdist2(haloX, haloY, haloZ, G['p'][:,0][Gclose], G['p'][:,1][Gclose], G['p'][:,2][Gclose], boxsize)
Ghalo = Gdists < Rvir

Sclose = OR.prepare_X_close(S, Rvir, haloX, haloY, haloZ)	
Sdists = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0][Sclose], S['p'][:,1][Sclose], S['p'][:,2][Sclose], boxsize)
Shalo = Sdists < Rvir

SVclose = Sdists[Shalo] < Rvir*InnerDivide

ShellCuts, CrossCuts, ShellThickness_ar = SF.shell_check_custom(Gdists[Ghalo], Rvir, R_bin_array_min , R_bin_array_max, a, h)

S_ShellCuts, S_CrossCuts, S_ShellThickness_ar = SF.shell_check_custom(Sdists[Shalo], Rvir, R_bin_array_min , R_bin_array_max, a, h)

TempCutHot = PhysTemp[Gclose][Ghalo] > OVIUpperBound
TempCutO6 = (PhysTemp[Gclose][Ghalo] <= OVIUpperBound) * (PhysTemp[Gclose][Ghalo] >= LowIonUpperBound)
TempCutLowIon = (PhysTemp[Gclose][Ghalo] <= LowIonUpperBound) * (PhysTemp[Gclose][Ghalo] >= ColdTempBound)
TempCutCold = (PhysTemp[Gclose][Ghalo] <= ColdTempBound)

bincount = 0
R_bins = len(ShellThickness_ar)
MetalMassHot = []
MetalMassO6 = []
MetalMassLowIon = []
MetalMassCold = []

OxyMassHot = []
OxyMassO6 = []
OxyMassLowIon = []
OxyMassCold = []


#divide things up in bins - this step is done even when divide between ISM and CGM is physical... it's useful to have it this way so that we can do full profile of where hot gas is. 
while (bincount < R_bins):
	BinCuts = ShellCuts[bincount]
	aMetalMassHot = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutHot] * Metallicity[Gclose][Ghalo][BinCuts * TempCutHot]) * 1e10 / h
	aMetalMassO6 = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutO6] * Metallicity[Gclose][Ghalo][BinCuts * TempCutO6]) * 1e10 / h
	aMetalMassLowIon = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutLowIon] * Metallicity[Gclose][Ghalo][BinCuts * TempCutLowIon]) * 1e10 / h
	aMetalMassCold = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutCold] * Metallicity[Gclose][Ghalo][BinCuts * TempCutCold]) * 1e10 / h

	aOxyMassHot = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutHot] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutHot]) * 1e10 / h
	aOxyMassO6 = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutO6] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutO6]) * 1e10 / h
	aOxyMassLowIon = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutLowIon] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutLowIon]) * 1e10 / h
	aOxyMassCold = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutCold] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutCold]) * 1e10 / h


	MetalMassHot.append(aMetalMassHot)
	MetalMassO6.append(aMetalMassO6)
	MetalMassLowIon.append(aMetalMassLowIon)
	MetalMassCold.append(aMetalMassCold)

	OxyMassHot.append(aOxyMassHot)
	OxyMassO6.append(aOxyMassO6)
	OxyMassLowIon.append(aOxyMassLowIon)
	OxyMassCold.append(aOxyMassCold)

	bincount+=1

MetalMassHot = np.array(MetalMassHot)
MetalMassO6 = np.array(MetalMassO6)
MetalMassLowIon = np.array(MetalMassLowIon)
MetalMassCold = np.array(MetalMassCold)

OxyMassHot = np.array(OxyMassHot)
OxyMassO6 = np.array(OxyMassO6)
OxyMassLowIon = np.array(OxyMassLowIon)
OxyMassCold = np.array(OxyMassCold)

S_ShellClose =  S_ShellCuts[0] + S_ShellCuts[1]
Gvclose = ShellCuts[0]

if (UsePhysicalDivide):
	Gvclose = Gdists[Ghalo]*a/h < PhysicalDivideISM
	Gfar = np.invert(Gvclose)
else:
	Gvclose =  ShellCuts[0]

MetalMassExpected = np.sum(S['m'][Sclose][Shalo][S_ShellClose])


S_ShellFar = np.invert(S_ShellClose)
MetalMassExpectedOuter = np.sum(S['m'][Sclose][Shalo][S_ShellFar])

theyield = 0.02
if (use_empirical_yield):
	theyield = EmpiricalYield 
boost_for_oldstars = 1.3

#when using emperical, aging is already kind of accounted for.. 
if (use_empirical_yield):
	MetalMassExpected *= 1e10 * theyield / h
	MetalMassExpectedOuter *= 1e10 * theyield / h
else:
	MetalMassExpected *= 1e10 * theyield * boost_for_oldstars / h
	MetalMassExpectedOuter *= 1e10 * theyield * boost_for_oldstars / h



MetalMassStars = np.sum(S['m'][Sclose][Shalo][S_ShellClose] * Stellar_met[Sclose][Shalo][S_ShellClose]) 
MetalMassStarsOuter = np.sum(S['m'][Sclose][Shalo][S_ShellFar] * Stellar_met[Sclose][Shalo][S_ShellFar]) 


OxyMassStars = np.sum(S['m'][Sclose][Shalo][S_ShellClose] * Stellar_oxygen[Sclose][Shalo][S_ShellClose] )

MetalMassStars *=1e10 / h
MetalMassStarsOuter *=1e10 / h

if (UsePhysicalDivide):
#same as when we were binning but now only separating the particles out
	MetalPartsHot = (G['m'][Gclose][Ghalo][TempCutHot] * Metallicity[Gclose][Ghalo][TempCutHot]) * 1e10 / h
	MetalPartsO6 = (G['m'][Gclose][Ghalo][TempCutO6] * Metallicity[Gclose][Ghalo][TempCutO6]) * 1e10 / h
	MetalPartsLowIon = (G['m'][Gclose][Ghalo][TempCutLowIon] * Metallicity[Gclose][Ghalo][TempCutLowIon]) * 1e10 / h
	MetalPartsCold = (G['m'][Gclose][Ghalo][TempCutCold] * Metallicity[Gclose][Ghalo][TempCutCold]) * 1e10 / h

	OxyPartsHot = (G['m'][Gclose][Ghalo][TempCutHot] * Gas_oxygen[Gclose][Ghalo][ TempCutHot]) * 1e10 / h
	OxyPartsO6 = (G['m'][Gclose][Ghalo][ TempCutO6] * Gas_oxygen[Gclose][Ghalo][ TempCutO6]) * 1e10 / h
	OxyPartsLowIon = (G['m'][Gclose][Ghalo][ TempCutLowIon] * Gas_oxygen[Gclose][Ghalo][ TempCutLowIon]) * 1e10 / h
	OxyPartsCold = (G['m'][Gclose][Ghalo][ TempCutCold] * Gas_oxygen[Gclose][Ghalo][ TempCutCold]) * 1e10 / h


	MetalMassISM = np.sum(MetalPartsHot[Gvclose[TempCutHot]]) + np.sum(MetalPartsO6[Gvclose[TempCutO6]]) + np.sum(MetalPartsLowIon[Gvclose[TempCutLowIon]]) + np.sum(MetalPartsCold[Gvclose[TempCutCold]])
	OxyMetalMassISM = np.sum(OxyPartsHot[Gvclose[TempCutHot]]) + np.sum(OxyPartsO6[Gvclose[TempCutO6]]) + np.sum(OxyPartsLowIon[Gvclose[TempCutLowIon]]) + np.sum(OxyPartsCold[Gvclose[TempCutCold]])
	
	
	G_not_vclose = np.invert(Gvclose)
	
	
	HaloMassHot = np.sum(MetalPartsHot[G_not_vclose[TempCutHot]])
	HaloMassO6 = np.sum(MetalPartsO6[G_not_vclose[TempCutO6]])
	HaloMassLowIon = np.sum(MetalPartsLowIon[G_not_vclose[TempCutLowIon]])
	HaloMassCold = np.sum(MetalPartsCold[G_not_vclose[TempCutCold]])
	
	OxyHaloMassHot = np.sum(OxyPartsHot[G_not_vclose[TempCutHot]])
	OxyHaloMassO6 = np.sum(OxyPartsO6[G_not_vclose[TempCutO6]])
	OxyHaloMassLowIon = np.sum(OxyPartsLowIon[G_not_vclose[TempCutLowIon]])
	OxyHaloMassCold = np.sum(OxyPartsCold[G_not_vclose[TempCutCold]])
	
else:
#when we can just use shells
	MetalMassISM = MetalMassHot[0] + MetalMassO6[0] + MetalMassLowIon[0] + MetalMassCold[0]
	OxyMetalMassISM = OxyMassHot[0] + OxyMassO6[0] + OxyMassLowIon[0] + OxyMassCold[0]

	HaloMassHot = np.sum(MetalMassHot[1:])
	HaloMassO6 = np.sum(MetalMassO6[1:])
	HaloMassLowIon = np.sum(MetalMassLowIon[1:])
	HaloMassCold = np.sum(MetalMassCold[1:])

	OxyHaloMassHot = np.sum(OxyMassHot[1:])
	OxyHaloMassO6 = np.sum(OxyMassO6[1:])
	OxyHaloMassLowIon = np.sum(OxyMassLowIon[1:])
	OxyHaloMassCold = np.sum(OxyMassCold[1:])

#if (UsePhysicalDivide):

MetalOutLine = [Mstar, MetalMassExpected, MetalMassStars, MetalMassISM, HaloMassHot, HaloMassO6, HaloMassLowIon, HaloMassCold, OxyYield, OxyMassStars, OxyMetalMassISM, OxyHaloMassHot, OxyHaloMassO6, OxyHaloMassLowIon, OxyHaloMassCold, theyield]

print 'here is your metals'
print SF.line_to_string(MetalOutLine, newline='n')

####################################
## stellar age metal distribution ##
####################################

amin_lim = 0.0
amax_lim = 1.0
num_bins = 200
histresults_mass = np.histogram(S['age'][Sclose][Shalo][SVclose], num_bins, range=(amin_lim, amax_lim), weights = S['m'][Sclose][Shalo][SVclose])

histresults_metalmass = np.histogram(S['age'][Sclose][Shalo][SVclose], num_bins, range=(amin_lim, amax_lim), weights = S['m'][Sclose][Shalo][SVclose]*Stellar_met[Sclose][Shalo][SVclose])

histlim1 = histresults_metalmass[1][:-1]
histlim2 = histresults_metalmass[1][1:]
cummetals = histresults_metalmass[0]
cummass = histresults_mass[0]

thetable = np.column_stack((histlim1, histlim2, cummetals, cummass))
foutname = today+'stellar_cum_table.txt'
#indMstar = np.cumsum(S['m'][Sclose][Shalo])
#indMstarClose = np.cumsum(S['m'][Sclose][Shalo][SVclose])
indMstar = np.sum(S['m'][Sclose][Shalo])
indMstarClose = np.sum(S['m'][Sclose][Shalo][SVclose])
theheader = [Mstar, indMstar, indMstarClose]
print 'header ',theheader
header_string = SF.line_to_string(theheader, newline='n')
np.savetxt(foutname, thetable, fmt='%1.6e', header=header_string)

