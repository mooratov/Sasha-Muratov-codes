#meant to get info about the halo in shells

#follow the same general pattern of NewSashaGP_exp.py
#but don't bother with outflow rate calculations; only get average properties of shells. 

#only counts outflows

import sys
import numpy as np
import Sasha_functions as SF
import outflow_rates as OR
from tasz import hubble_param
from readsnap import readsnap
from datetime import date
fastwinds = True
today = str(date.today())
#today = '2016-03-10'

do_graphics = False
year = 3.15569e7
pc = 3.08567758e18
kpc = pc * 1e3

halo_to_do = 0
multifile = False
use_fixed_halos = 0
readDM = True
use_empirical_yield = True
plotcontam = False

InnerDivide = 0.1 #In units of Rvir
UsePhysicalDivide = False
PhysicalDivideISM = 10.0  #in units of physical kpc

ColdTempBound = 10**4.0
LowIonUpperBound = 10**4.7
OVIUpperBound = 10**5.3

use_OuterBoundary = False 
OuterBoundary_physical = 300 

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

OxyYield = (np.sum(Stellar_oxygen*S['m']) + np.sum(Gas_oxygen*G['m'])) / np.sum(S['m'])

EmpiricalYield = Total_metal_mass / np.sum(S['m'])

a = G['header'][2]
z = 1.0/a - 1.0
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

print Rvir, Vsig, M, 'before'
Rvir, Vsig, M = SF.check_Rvir_growth(halo_to_do, a, Rvir, Vsig, M, therod=use_fixed_halos)
print Rvir, Vsig, M, 'after'

if (use_OuterBoundary):
	Rvir = OuterBoundary_physical * h / a

#setting the locations of each spherical shell, as a fraction of Rvir
R_bin_array_min = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
R_bin_array_max = R_bin_array_min+0.1

Gclose = OR.prepare_X_close(G, Rvir, haloX, haloY, haloZ)	
Gdists = SF.calcdist2(haloX, haloY, haloZ, G['p'][:,0][Gclose], G['p'][:,1][Gclose], G['p'][:,2][Gclose], boxsize)
Ghalo = Gdists < Rvir



#Hubble = hubble_param(a, omega_matter, h)

GprelX = G['p'][:,0][Gclose][Ghalo] - haloX
GprelY = G['p'][:,1][Gclose][Ghalo] - haloY
GprelZ = G['p'][:,2][Gclose][Ghalo] - haloZ	
GvrelX = G['v'][:,0][Gclose][Ghalo] * np.sqrt(a) - haloVX + GprelX*a*Hubble/(h*1000)
GvrelY = G['v'][:,1][Gclose][Ghalo] * np.sqrt(a) - haloVY + GprelY*a*Hubble/(h*1000)
GvrelZ = G['v'][:,2][Gclose][Ghalo] * np.sqrt(a) - haloVZ +  GprelZ*a*Hubble/(h*1000)
#rescale positions so that magnitude = 1 (unit vectors)
GprelX_RS = GprelX / Gdists[Ghalo]
GprelY_RS = GprelY / Gdists[Ghalo]
GprelZ_RS = GprelZ / Gdists[Ghalo]
v_rad = GprelX_RS*GvrelX  + GprelY_RS*GvrelY + GprelZ_RS*GvrelZ
#print 'vrad ',v_rad
#cut particles that are "outflowing"
if (fastwinds):
	OutflowCut = v_rad > Vsig/np.sqrt(3) 
else: 
	OutflowCut = v_rad > 0
#Ghalo *= OutflowCut
Ghalo[np.invert(OutflowCut)] = False

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

GasMassHot = []
GasMassO6 = []
GasMassLowIon = []
GasMassCold = []

GasFluxHot = []
GasFluxO6 = []
GasFluxLowIon = []
GasFluxCold = []

shellprofilear = [Nsnap, z]




#divide things up in bins - this step is done even when divide between ISM and CGM is physical... it's useful to have it this way so that we can do full profile of where hot gas is. 
while (bincount < R_bins):
	BinCuts = ShellCuts[bincount]
	aMetalMassHot = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutHot] * Metallicity[Gclose][Ghalo][BinCuts * TempCutHot]) * 1e10 / h
	aMetalMassO6 = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutO6] * Metallicity[Gclose][Ghalo][BinCuts * TempCutO6]) * 1e10 / h
	aMetalMassLowIon = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutLowIon] * Metallicity[Gclose][Ghalo][BinCuts * TempCutLowIon]) * 1e10 / h
	aMetalMassCold = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutCold] * Metallicity[Gclose][Ghalo][BinCuts * TempCutCold]) * 1e10 / h

	aGasMassHot = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutHot]) * 1e10 / h
	aGasMassO6 = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutO6]) * 1e10 / h
	aGasMassLowIon = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutLowIon]) * 1e10 / h
	aGasMassCold = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutCold]) * 1e10 / h
	
	aGasFluxHot = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutHot] * v_rad[BinCuts * TempCutHot]*1e5*year/kpc ) * 1e10 / h
	aGasFluxO6 = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutO6] * v_rad[BinCuts * TempCutO6]*1e5*year/kpc ) * 1e10 / h
	aGasFluxLowIon = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutLowIon] * v_rad[BinCuts * TempCutLowIon]*1e5*year/kpc ) * 1e10 / h
	aGasFluxCold = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutCold] * v_rad[BinCuts * TempCutCold]*1e5*year/kpc ) * 1e10 / h

	aOxyMassHot = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutHot] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutHot]) * 1e10 / h
	aOxyMassO6 = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutO6] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutO6]) * 1e10 / h
	aOxyMassLowIon = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutLowIon] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutLowIon]) * 1e10 / h
	aOxyMassCold = np.sum(G['m'][Gclose][Ghalo][BinCuts * TempCutCold] * Gas_oxygen[Gclose][Ghalo][BinCuts * TempCutCold]) * 1e10 / h



	MetalMassHot.append(aMetalMassHot)
	MetalMassO6.append(aMetalMassO6)
	MetalMassLowIon.append(aMetalMassLowIon)
	MetalMassCold.append(aMetalMassCold)

	GasMassHot.append(aGasMassHot)
	GasMassO6.append(aGasMassO6)
	GasMassLowIon.append(aGasMassLowIon)
	GasMassCold.append(aGasMassCold)

	GasFluxHot.append(aGasFluxHot)
	GasFluxO6.append(aGasFluxO6)
	GasFluxLowIon.append(aGasFluxLowIon)
	GasFluxCold.append(aGasFluxCold)


	OxyMassHot.append(aOxyMassHot)
	OxyMassO6.append(aOxyMassO6)
	OxyMassLowIon.append(aOxyMassLowIon)
	OxyMassCold.append(aOxyMassCold)


	shellprofilear.append(bincount)
	shellprofilear.append(aGasMassHot)
	shellprofilear.append(aGasMassO6)
	shellprofilear.append(aGasMassLowIon)
	shellprofilear.append(aGasMassCold)
	shellprofilear.append(aGasFluxHot)
	shellprofilear.append(aGasFluxO6)
	shellprofilear.append(aGasFluxLowIon)
	shellprofilear.append(aGasFluxCold)

	bincount+=1

foutname = today+'PhaseShells.txt'
outline3 = SF.line_to_string(shellprofilear, newline='y')
i = open(foutname, 'a')
i.write(outline3)
i.close()	



MetalMassHot = np.array(MetalMassHot)
MetalMassO6 = np.array(MetalMassO6)
MetalMassLowIon = np.array(MetalMassLowIon)
MetalMassCold = np.array(MetalMassCold)

GasMassHot = np.array(GasMassHot)
GasMassO6 = np.array(GasMassO6)
GasMassLowIon = np.array(GasMassLowIon)
GasMassCold = np.array(GasMassCold)

GasFluxHot = np.array(GasFluxHot)
GasFluxO6 = np.array(GasFluxO6)
GasFluxLowIon = np.array(GasFluxLowIon)
GasFluxCold = np.array(GasFluxCold)


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

TotalStellarHaloMass = np.sum(S['m'][Sclose][Shalo])
TotalStellarHaloMassClose = np.sum(S['m'][Sclose][Shalo][S_ShellClose])

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

#baryon fraction	


if (UsePhysicalDivide):
#same as when we were binning but now only separating the particles out
	MetalPartsHot = (G['m'][Gclose][Ghalo][TempCutHot] * Metallicity[Gclose][Ghalo][TempCutHot]) * 1e10 / h
	MetalPartsO6 = (G['m'][Gclose][Ghalo][TempCutO6] * Metallicity[Gclose][Ghalo][TempCutO6]) * 1e10 / h
	MetalPartsLowIon = (G['m'][Gclose][Ghalo][TempCutLowIon] * Metallicity[Gclose][Ghalo][TempCutLowIon]) * 1e10 / h
	MetalPartsCold = (G['m'][Gclose][Ghalo][TempCutCold] * Metallicity[Gclose][Ghalo][TempCutCold]) * 1e10 / h

	GasPartsHot = (G['m'][Gclose][Ghalo][TempCutHot] ) * 1e10 / h
	GasPartsO6 = (G['m'][Gclose][Ghalo][TempCutO6] ) * 1e10 / h
	GasPartsLowIon = (G['m'][Gclose][Ghalo][TempCutLowIon] ) * 1e10 / h
	GasPartsCold = (G['m'][Gclose][Ghalo][TempCutCold] ) * 1e10 / h


	OxyPartsHot = (G['m'][Gclose][Ghalo][TempCutHot] * Gas_oxygen[Gclose][Ghalo][ TempCutHot]) * 1e10 / h
	OxyPartsO6 = (G['m'][Gclose][Ghalo][ TempCutO6] * Gas_oxygen[Gclose][Ghalo][ TempCutO6]) * 1e10 / h
	OxyPartsLowIon = (G['m'][Gclose][Ghalo][ TempCutLowIon] * Gas_oxygen[Gclose][Ghalo][ TempCutLowIon]) * 1e10 / h
	OxyPartsCold = (G['m'][Gclose][Ghalo][ TempCutCold] * Gas_oxygen[Gclose][Ghalo][ TempCutCold]) * 1e10 / h


	MetalMassISM = np.sum(MetalPartsHot[Gvclose[TempCutHot]]) + np.sum(MetalPartsO6[Gvclose[TempCutO6]]) + np.sum(MetalPartsLowIon[Gvclose[TempCutLowIon]]) + np.sum(MetalPartsCold[Gvclose[TempCutCold]])
	OxyMetalMassISM = np.sum(OxyPartsHot[Gvclose[TempCutHot]]) + np.sum(OxyPartsO6[Gvclose[TempCutO6]]) + np.sum(OxyPartsLowIon[Gvclose[TempCutLowIon]]) + np.sum(OxyPartsCold[Gvclose[TempCutCold]])
	GasMassISM = np.sum(GasPartsHot[Gvclose[TempCutHot]]) + np.sum(GasPartsO6[Gvclose[TempCutO6]]) + np.sum(GasPartsLowIon[Gvclose[TempCutLowIon]]) + np.sum(GasPartsCold[Gvclose[TempCutCold]])
	
	G_not_vclose = np.invert(Gvclose)
	
	
	#for just metals
	HaloMassHot = np.sum(MetalPartsHot[G_not_vclose[TempCutHot]])
	HaloMassO6 = np.sum(MetalPartsO6[G_not_vclose[TempCutO6]])
	HaloMassLowIon = np.sum(MetalPartsLowIon[G_not_vclose[TempCutLowIon]])
	HaloMassCold = np.sum(MetalPartsCold[G_not_vclose[TempCutCold]])
	
	#for gas and metals
	HaloMassHotGas = np.sum(GasPartsHot[G_not_vclose[TempCutHot]])
	HaloMassO6Gas = np.sum(GasPartsO6[G_not_vclose[TempCutO6]])
	HaloMassLowIonGas = np.sum(GasPartsLowIon[G_not_vclose[TempCutLowIon]])
	HaloMassColdGas = np.sum(GasPartsCold[G_not_vclose[TempCutCold]])	
	
	
	OxyHaloMassHot = np.sum(OxyPartsHot[G_not_vclose[TempCutHot]])
	OxyHaloMassO6 = np.sum(OxyPartsO6[G_not_vclose[TempCutO6]])
	OxyHaloMassLowIon = np.sum(OxyPartsLowIon[G_not_vclose[TempCutLowIon]])
	OxyHaloMassCold = np.sum(OxyPartsCold[G_not_vclose[TempCutCold]])
	
	CGMMass = np.sum(G['m'][Gclose][Ghalo][G_not_vclose])*1e10/h # dont need to divide by little h. Done already.
	ISMMass = np.sum(G['m'][Gclose][Ghalo][Gvclose])*1e10/h
	TotalGMass = (ISMMass + CGMMass)
	
	ISMMetalMassCold = np.sum(MetalPartsCold)
	
	GasMassCold = np.sum(GasPartsCold)
	
else:
#when we can just use shells
	MetalMassISM = MetalMassHot[0] + MetalMassO6[0] + MetalMassLowIon[0] + MetalMassCold[0]
	OxyMetalMassISM = OxyMassHot[0] + OxyMassO6[0] + OxyMassLowIon[0] + OxyMassCold[0]
	
	ISMMetalMassCold = np.sum(MetalMassCold)
	
	GasMassCold = np.sum(GasMassCold)

	

	HaloMassHot = np.sum(MetalMassHot[1:])
	HaloMassO6 = np.sum(MetalMassO6[1:])
	HaloMassLowIon = np.sum(MetalMassLowIon[1:])
	HaloMassCold = np.sum(MetalMassCold[1:])

	HaloMassHotGas = np.sum(GasMassHot[1:])
	HaloMassO6Gas = np.sum(GasMassO6[1:])
	HaloMassLowIonGas = np.sum(MetalMassLowIon[1:])
	HaloMassColdGas = np.sum(MetalMassCold[1:])


	OxyHaloMassHot = np.sum(OxyMassHot[1:])
	OxyHaloMassO6 = np.sum(OxyMassO6[1:])
	OxyHaloMassLowIon = np.sum(OxyMassLowIon[1:])
	OxyHaloMassCold = np.sum(OxyMassCold[1:])
	G_not_vclose = np.invert(Gvclose)

	CGMMass = np.sum(G['m'][Gclose][Ghalo][G_not_vclose])*1e10/h  
	ISMMass = np.sum(G['m'][Gclose][Ghalo][Gvclose])*1e10/h
	GasMassISM = ISMMass
	TotalGMass = (ISMMass + CGMMass)


coldISMMetallicity = ISMMetalMassCold/GasMassCold


if (use_OuterBoundary and  UsePhysicalDivide):
	#extra radial analysis
	shellprofilear = [Nsnap, z]
	newshellbins_min = np.arange(OuterBoundary_physical/10.0) * 10.0
	newshellbins_max = newshellbins_min + 10.0
	newshellmidpoints = (newshellbins_min + newshellbins_max)/2.0
	bincount = 0
	GasMasser = []
	MetalMasser = []
	OxyMasser = []
	
	while (bincount < len(newshellbins_min)):
		Bincuts = (Gdists[Ghalo]*a/h > newshellbins_min[bincount]) *  (Gdists[Ghalo]*a/h < newshellbins_max[bincount])
		TotalMassInPhysShell = np.sum(G['m'][Gclose][Ghalo][Bincuts])*1e10/h
		TotalMetalInPhysShell = np.sum((G['m'][Gclose][Ghalo][Bincuts] * Metallicity[Gclose][Ghalo][Bincuts])) * 1e10 / h
		TotalOxyInPhysShell = np.sum((G['m'][Gclose][Ghalo][Bincuts] * Gas_oxygen[Gclose][Ghalo][Bincuts])) * 1e10 / h
		TotalMetal6InPhysShell = np.sum((G['m'][Gclose][Ghalo][TempCutO6 * Bincuts] * Metallicity[Gclose][Ghalo][TempCutO6 * Bincuts])) * 1e10 / h
		TotalOxy6InPhysShell =  np.sum((G['m'][Gclose][Ghalo][TempCutO6 * Bincuts] * Gas_oxygen[Gclose][Ghalo][TempCutO6 * Bincuts])) * 1e10 / h

		TotalMetalHotInPhysShell = np.sum((G['m'][Gclose][Ghalo][TempCutHot * Bincuts] * Metallicity[Gclose][Ghalo][TempCutHot * Bincuts])) * 1e10 / h
		TotalOxyHotInPhysShell =  np.sum((G['m'][Gclose][Ghalo][TempCutHot * Bincuts] * Gas_oxygen[Gclose][Ghalo][TempCutHot * Bincuts])) * 1e10 / h


		GasMasser.append(TotalMassInPhysShell) 
		MetalMasser.append(TotalMetalInPhysShell)
		OxyMasser.append(TotalOxyInPhysShell)
		shellprofilear.append(bincount)
		shellprofilear.append(newshellbins_min[bincount])
		shellprofilear.append(TotalMassInPhysShell)
		shellprofilear.append(TotalMetalInPhysShell)
		shellprofilear.append(TotalOxyInPhysShell)
		shellprofilear.append(TotalMetal6InPhysShell)
		shellprofilear.append(TotalOxy6InPhysShell)
		shellprofilear.append(TotalMetalHotInPhysShell)
		shellprofilear.append(TotalOxyHotInPhysShell)
		#can add phase stuff later
		bincount+=1
	foutname = today+'MetalShells.txt'
	outline3 = SF.line_to_string(shellprofilear, newline='y')
	i = open(foutname, 'a')
	i.write(outline3)
	i.close()	
		
#if (UsePhysicalDivide):

TotalMetalMass =  MetalMassStars + MetalMassISM + HaloMassHot + HaloMassO6 + HaloMassLowIon + HaloMassCold

MetalOutLine = [Mstar, MetalMassExpected, MetalMassStars, MetalMassISM, HaloMassHot, HaloMassO6, HaloMassLowIon, HaloMassCold, OxyYield, OxyMassStars, OxyMetalMassISM, OxyHaloMassHot, OxyHaloMassO6, OxyHaloMassLowIon, OxyHaloMassCold, theyield,a, Nsnap]

GasOutLine = [Mstar, Mgas, Mstar, GasMassISM, HaloMassHotGas, HaloMassO6Gas, HaloMassLowIonGas, HaloMassColdGas, theyield,a, Nsnap]


foutname = today+'MetalProfile.txt'
print 'here is your metals'
print SF.line_to_string(MetalOutLine, newline='n')
outline2 = SF.line_to_string(MetalOutLine, newline='y')
g = open(foutname, 'a')
g.write(outline2)
g.close()

foutname = today+'GasProfile.txt'
print 'here is your gasses'
print SF.line_to_string(GasOutLine, newline='n')
outline2 = SF.line_to_string(GasOutLine, newline='y')
g = open(foutname, 'a')
g.write(outline2)
g.close()



baryonfrac = (Mgas/h + Mstar/h)/ (M*0.17/h)
baryonfrac_t2 = (TotalGMass + (Mstar/h))/(M*0.17/h)

metalfrac = TotalMetalMass / MetalMassExpected 
StellarMetallicity = MetalMassStars / (Mstar/h)
ISMMetallicity = MetalMassISM / (ISMMass)
CGMMetallicity = ( HaloMassHot + HaloMassO6 + HaloMassLowIon + HaloMassCold)/ (CGMMass)

Line_for_Baryon_frac = [M, Mstar, Mgas, ISMMass, CGMMass, TotalGMass,TotalMetalMass,MetalMassExpected,  baryonfrac, baryonfrac_t2,metalfrac,StellarMetallicity, ISMMetallicity, CGMMetallicity , a, Nsnap]
foutname = today+'BaryonProfile.txt'
print 'here is line 2'
#print Line_for_Baryon_frac
print SF.line_to_string(Line_for_Baryon_frac, newline='n')
outline2 = SF.line_to_string(Line_for_Baryon_frac, newline='y')
g = open(foutname, 'a')
g.write(outline2)
g.close()


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
foutname = today+'stellar_cum_tableN'+str(halo_to_do)+'.txt'
#indMstar = np.cumsum(S['m'][Sclose][Shalo])
#indMstarClose = np.cumsum(S['m'][Sclose][Shalo][SVclose])
indMstar = np.sum(S['m'][Sclose][Shalo])
indMstarClose = np.sum(S['m'][Sclose][Shalo][SVclose])
theheader = [Mstar, indMstar, indMstarClose]
print 'header ',theheader
header_string = SF.line_to_string(theheader, newline='n')
np.savetxt(foutname, thetable, fmt='%1.6e', header=header_string)

print 'Stellar Metallicity, ISM Metallicity, CGMMetallicity, MetalMassISM',StellarMetallicity, ISMMetallicity, CGMMetallicity, "{0:.3e}".format(MetalMassISM)
print 'Cold ISM Metallicity, cold metal mass',coldISMMetallicity, "{0:.3e}".format(ISMMetalMassCold)

if (do_graphics):
	import graphics_library as GL
	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutCold
	foutname1 = str(today)+'coldgas'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts,  dolabel=False)

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutHot
	foutname1 = str(today)+'hotgas'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts,  dolabel=False)

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutO6
	foutname1 = str(today)+'O6gas'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts, dolabel=False)

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutLowIon
	foutname1 = str(today)+'LowIongas'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, Rvir, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts,  dolabel=False)

	R30 = 30.0*h / a

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutCold
	foutname1 = str(today)+'coldgasR30'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, R30, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts,  dolabel=False)

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutHot
	foutname1 = str(today)+'hotgasR30'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, R30, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts, dolabel=False)

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutO6
	foutname1 = str(today)+'O6gasR30'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, R30, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts, dolabel=False)

	sendCuts = np.copy(Gclose) 
	sendCuts[sendCuts] *= Ghalo
	sendCuts[sendCuts] *= TempCutLowIon
	foutname1 = str(today)+'LowIongasR30'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(G, R30, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts,  dolabel=False)

print 'stellar masses', Mstar, TotalStellarHaloMass, TotalStellarHaloMassClose

if (readDM and plotcontam):
	#sendCuts = np.copy(Gclose) 
	#sendCuts[sendCuts] *= Ghalo
	#sendCuts[sendCuts] *= TempCutLowIon
	if (P2['k'] > 0):
		NP0, NP1, NP2, NP3, NP4, NP5 = P2['header'][6]
		M0, M1, M2, M3, M4, M5 = P2['header'][1]
		print 'P2 ',P2
		if (M2 > 0):
			print 'adjusting P2 mass'
			P2['m'] = np.zeros(NP2)
			P2['m'][0:NP2] = M2

		P2close = OR.prepare_X_close(P2, Rvir, haloX, haloY, haloZ)	
		P2dists = SF.calcdist2(haloX, haloY, haloZ, P2['p'][:,0][P2close], P2['p'][:,1][P2close], P2['p'][:,2][P2close], boxsize)
		P2halo = P2dists < Rvir
		sendCuts = np.copy(P2close) 
		sendCuts[sendCuts] *= P2halo
		print 'total P2 mass ', np.sum(P2['m'][P2close][P2halo])*1e10/h, len(P2['m'][P2close][P2halo])

		[sat1x, sat1y, sat1z] = [51910.3198828,   50545.6714349,   50456.1511434]
		P2sat1 = OR.prepare_X_close(P2, 17.0, sat1x, sat1y, sat1z)	
		P2sat1dists = SF.calcdist2(sat1x, sat1y, sat1z, P2['p'][:,0][P2sat1], P2['p'][:,1][P2sat1], P2['p'][:,2][P2sat1], boxsize)
		P2sat1halo = P2sat1dists < Rvir
		print 'total P2 sat1 mass ', np.sum(P2['m'][P2sat1][P2sat1halo])*1e10/h, len(P2['m'][P2sat1][P2sat1halo])

		foutname1 = str(today)+'P2contam'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
		GL.pluscuts_2d_desnity_hist(P2, Rvir, haloX, haloY, haloZ,  foutname1,  dolabel=False, Cuts=sendCuts, thevmin = 0.1, thevmax = 1e2)	
		GL.pluscuts_2d_desnity_hist(P2, Rvir, haloX, haloY, haloZ,  foutname1,  dolabel=False, Cuts=sendCuts, thevmin = 0.1, thevmax = 1e2, extrapoints = [[51910.3198828, 50545.6714349], [51793.9792632, 50475.6552445], [51896.4042171, 50505.0213578], [ 51928.9067802, 50499.5665531],  [51880.645855, 50558.8252286] ])	

	if (P3['k'] > 0):
		NP0, NP1, NP2, NP3, NP4, NP5 = P3['header'][6]
		M0, M1, M2, M3, M4, M5 = P3['header'][1]
		if (M2 > 0):
			print 'adjusting P3 mass'
			P3['m'] = np.zeros(NP3)
			P3['m'][0:NP3] = M3
		print 'P3 ',P3
		P3close = OR.prepare_X_close(P3, Rvir, haloX, haloY, haloZ)	
		P3dists = SF.calcdist2(haloX, haloY, haloZ, P3['p'][:,0][P3close], P3['p'][:,1][P3close], P3['p'][:,2][P3close], boxsize)
		P3halo = P3dists < Rvir
		sendCuts = np.copy(P3close) 
		sendCuts[sendCuts] *= P3halo
		print 'total P3 mass ', np.sum(P3['m'][P3close][P3halo])*1e10/h, len(P3['m'][P3close][P3halo])

		foutname1 = str(today)+'P3contam'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
		GL.pluscuts_2d_desnity_hist(P3, Rvir, haloX, haloY, haloZ,  foutname1,  dolabel=False, Cuts=sendCuts, thevmin = 0.1, thevmax = 1e2)	
	if (P5['k'] > 0):
		NP0, NP1, NP2, NP3, NP4, NP5 = P5['header'][6]
		M0, M1, M2, M3, M4, M5 = P5['header'][1]
		print 'P5 ',P5
		if (M5 > 0):
			print 'adjusting P5 mass'
			P5['m'] = np.zeros(NP5)
			P5['m'][0:NP5] = M5

		print 'P5 ',P5
		P5close = OR.prepare_X_close(P5, Rvir, haloX, haloY, haloZ)	
		P5dists = SF.calcdist2(haloX, haloY, haloZ, P5['p'][:,0][P5close], P5['p'][:,1][P5close], P5['p'][:,2][P5close], boxsize)
		P5halo = P5dists < Rvir
		sendCuts = np.copy(P5close) 
		sendCuts[sendCuts] *= P5halo
		print 'total P5 mass ', np.sum(P5['m'][P5close][P5halo])*1e10/h, len(P5['m'][P5close][P5halo])

		foutname1 = str(today)+'P5contam'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
		GL.pluscuts_2d_desnity_hist(P5, Rvir, haloX, haloY, haloZ,  foutname1,  dolabel=False, Cuts=sendCuts, thevmin = 0.1, thevmax = 1e2)	
		
	NP0, NP1, NP2, NP3, NP4, NP5 = DM['header'][6]
	M0, M1, M2, M3, M4, M5 = DM['header'][1]
	if (M1 > 0):
		print 'adjusting DM mass'
		DM['m'] = np.zeros(NP1)
		DM['m'][0:NP1] = M1

	print 'DM ',DM
	DMclose = OR.prepare_X_close(DM, Rvir, haloX, haloY, haloZ)	
	DMdists = SF.calcdist2(haloX, haloY, haloZ, DM['p'][:,0][DMclose], DM['p'][:,1][DMclose], DM['p'][:,2][DMclose], boxsize)
	DMhalo = DMdists < Rvir
	sendCuts = np.copy(DMclose) 
	sendCuts[sendCuts] *= DMhalo
	print 'total DM mass ', np.sum(DM['m'][DMclose][DMhalo])*1e10/h
	foutname1 = str(today)+'DM'+Nsnapstring+'N'+str(halo_to_do)+'.pdf'
	GL.pluscuts_2d_desnity_hist(DM, Rvir, haloX, haloY, haloZ,  foutname1, Cuts=sendCuts,  dolabel=False, thevmin = 0.1, thevmax = 1e2)	
