import sys
import numpy as np
import Sasha_functions as SF
import outflow_rates as OR
from tasz import hubble_param
from readsnap import readsnap
from datetime import date
today = str(date.today())
#today = '2016-03-10'

do_graphics = False
year = 3.15569e7
pc = 3.08567758e18
kpc = pc * 1e3

halo_to_do = 0
multifile = True
use_fixed_halos = 5


TempBound1 = 10**4.0
TempBound2 = 10**5.0
TempBound3 = 10**6.0

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

Gclose = OR.prepare_X_close(G, Rvir, haloX, haloY, haloZ)	
Gdists = SF.calcdist2(haloX, haloY, haloZ, G['p'][:,0][Gclose], G['p'][:,1][Gclose], G['p'][:,2][Gclose], boxsize)
Ghalo = Gdists < Rvir

TempCutISM = PhysTemp[Gclose][Ghalo] < TempBound1
TempCutCGMCool = PhysTemp[Gclose][Ghalo] < TempBound2
TempCutCGMWarm = (PhysTemp[Gclose][Ghalo] >= TempBound2) * (PhysTemp[Gclose][Ghalo] < TempBound3)
TempCutCGMHot =  PhysTemp[Gclose][Ghalo] >= TempBound3

ISMDensCut = PhysRho[Gclose][Ghalo] > 0.1 

MISM1 = np.sum(G['m'][Gclose][Ghalo][ISMDensCut]) * 1e10 / h
MISM2 = np.sum(G['m'][Gclose][Ghalo][ISMDensCut * TempCutISM]) * 1e10 / h 

MCGM1 = np.sum(G['m'][Gclose][Ghalo][TempCutCGMCool * np.invert(ISMDensCut)])  * 1e10 / h
MCGM2 = np.sum(G['m'][Gclose][Ghalo][TempCutCGMCool * np.invert(TempCutISM * ISMDensCut)])  * 1e10 / h 
MCGM3 = np.sum(G['m'][Gclose][Ghalo][TempCutCGMWarm]) * 1e10 / h
MCGM4 = np.sum(G['m'][Gclose][Ghalo][TempCutCGMHot]) * 1e10 / h


Line = [M/h, Mstar/h, MISM1, MISM2, MCGM1, MCGM2, MCGM3, MCGM4, a, Nsnap]

foutname = today+'TheramlProfile.txt'
print SF.line_to_string(Line, newline='n')
outline2 = SF.line_to_string(Line, newline='y')
g = open(foutname, 'a')
g.write(outline2)
g.close()
