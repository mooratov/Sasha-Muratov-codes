#take a batch of particles detected as outflows at 0.25 Rvir at snapshot N1 
#follow them to snapshot N2, where there was a strong peak of outflow rate at 1.0 Rvir

import sys
import numpy as np
import outflow_rates as OR
import particle_tracking_functions as PTF
from readsnap import readsnap
from tasz import hubble_param
import Sasha_functions as SF

multifile = True
use_v0 = False 
halo_to_do  = 0
delay_snaps = 4
index1 = 2 
index2 = 9 
use_fixed_halos = 4


if (len(sys.argv) < 3):
	print 'syntax: blah.py Nsnapshot delay_snaps'
	sys.exit()

if (int(sys.argv[1]) < 10):
	Nsnapstring = '00'+str(sys.argv[1])
elif (int(sys.argv[1]) < 100):
	Nsnapstring = '0'+str(sys.argv[1])
else:
	Nsnapstring = str(sys.argv[1])
Nsnap = int(sys.argv[1])
delay_snaps = int(sys.argv[2])

if (len(sys.argv) > 3):
	index1 = int(sys.argv[3])
	index2 = int(sys.argv[4])

the_snapdir = './hdf5/'
if (multifile):
	the_snapdir = './snapdir_'+Nsnapstring+'/'
the_prefix ='snapshot'
the_suffix = '.hdf5'

print 'reading particle type 0 (gas?)' 
G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)



a = G['header'][2]
redshift = G['header'][3]
redshift1 = redshift
boxsize = G['header'][9]
omega_matter = G['header'][10]
omega_L = G['header'][11]
h = G['header'][12]
redshiftstring = "{0:.3f}".format(redshift)
Hubble = hubble_param(a, omega_matter, h)

halostats = SF.find_halo_now(halo_to_do, a, therod=use_fixed_halos)
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


Rvir, Vsig, M = SF.check_Rvir_growth(haloN, a, Rvir, Vsig, M, therod=use_fixed_halos)


#setting the locations of each spherical shell, as a fraction of Rvir
R_bin_array_min = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
R_bin_array_max = R_bin_array_min+0.1

Gclose = OR.prepare_X_close(G, Rvir, haloX, haloY, haloZ)	
dists = SF.calcdist2(haloX, haloY, haloZ, G['p'][:,0][Gclose], G['p'][:,1][Gclose], G['p'][:,2][Gclose], boxsize)
Ghalo = dists < Rvir

#bin particles according to spherical shells
ShellCuts, CrossCuts, ShellThickness_ar = SF.shell_check_custom(dists[Ghalo], Rvir, R_bin_array_min , R_bin_array_max, a, h)

#determine radial velocity
GprelX = G['p'][:,0][Gclose][Ghalo] - haloX
GprelY = G['p'][:,1][Gclose][Ghalo] - haloY
GprelZ = G['p'][:,2][Gclose][Ghalo] - haloZ	
GvrelX = G['v'][:,0][Gclose][Ghalo] * np.sqrt(a) - haloVX + GprelX*a*Hubble/(h*1000)
GvrelY = G['v'][:,1][Gclose][Ghalo] * np.sqrt(a) - haloVY + GprelY*a*Hubble/(h*1000)
GvrelZ = G['v'][:,2][Gclose][Ghalo] * np.sqrt(a) - haloVZ +  GprelZ*a*Hubble/(h*1000)
#rescale positions so that magnitude = 1 (unit vectors)
GprelX_RS = GprelX / dists[Ghalo]
GprelY_RS = GprelY / dists[Ghalo]
GprelZ_RS = GprelZ / dists[Ghalo]
v_rad = GprelX_RS*GvrelX  + GprelY_RS*GvrelY + GprelZ_RS*GvrelZ
print 'vrad ',v_rad
#cut particles that are "outflowing"

OutflowCut = v_rad > Vsig/np.sqrt(3)  
if (use_v0):
	OutflowCut = v_rad > 0

OutflowingParticlesID = np.copy(G['id'][Gclose][Ghalo][OutflowCut])
OutflowingParticleMass = np.copy(G['m'][Gclose][Ghalo][OutflowCut])
OutflowingParticle_vrad = np.copy(v_rad[OutflowCut])

OutflowingParticlesID_025Rvir = np.copy(G['id'][Gclose][Ghalo][OutflowCut*ShellCuts[index1]])
OutflowingParticleMass_025Rvir = np.copy(G['m'][Gclose][Ghalo][OutflowCut*ShellCuts[index1]])
OutflowingParticle_vrad_025Rvir = np.copy(v_rad[OutflowCut*ShellCuts[index1]])
Particles_G_025 = np.copy(G['id'][Gclose][Ghalo][ShellCuts[index1]])


print 'length of 0.25 outflowing particles',len(OutflowingParticlesID_025Rvir)
len25= len(OutflowingParticlesID_025Rvir)

original_outflow_rate = OR.outflow_flux_atshell(G['m'][Gclose][Ghalo], v_rad, OutflowCut, ShellCuts[index1], ShellThickness_ar[index1], h)

del G, GprelX, GprelY, GprelZ, GvrelX, GvrelY, GvrelZ, GprelX_RS, GprelY_RS, GprelZ_RS, v_rad 
#free up memory 





######################################################################on to next snapshot
######################################################################

Nsnap2 = Nsnap + delay_snaps
if (Nsnap2 < 10):
	Nsnapstring = '00'+str(Nsnap2)
elif (Nsnap2 < 100):
	Nsnapstring = '0'+str(Nsnap2)
else:
	Nsnapstring = str(Nsnap2)

the_snapdir = './hdf5/'
if (multifile):
	the_snapdir = './snapdir_'+Nsnapstring+'/'
the_prefix ='snapshot'
the_suffix = '.hdf5'

print 'reading particle type 0 (gas?)' 
G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)

print Nsnap, Nsnap2, 


a = G['header'][2]
redshift = G['header'][3]
redshift2 = redshift
boxsize = G['header'][9]
omega_matter = G['header'][10]
omega_L = G['header'][11]
h = G['header'][12]
redshiftstring = "{0:.3f}".format(redshift)
Hubble = hubble_param(a, omega_matter, h)

halostats = SF.find_halo_now(halo_to_do, a, therod=use_fixed_halos)

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


Rvir, Vsig, M = SF.check_Rvir_growth(haloN, a, Rvir, Vsig, M, therod=use_fixed_halos)


#compute outflow rates
R_bin_array_min = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
R_bin_array_max = R_bin_array_min+0.1

Gclose = OR.prepare_X_close(G, Rvir, haloX, haloY, haloZ)	
dists = SF.calcdist2(haloX, haloY, haloZ, G['p'][:,0][Gclose], G['p'][:,1][Gclose], G['p'][:,2][Gclose], boxsize)
Ghalo = dists < Rvir

#bin particles according to spherical shells
ShellCuts, CrossCuts, ShellThickness_ar = SF.shell_check_custom(dists[Ghalo], Rvir, R_bin_array_min , R_bin_array_max, a, h)

#determine radial velocity
GprelX = G['p'][:,0][Gclose][Ghalo] - haloX
GprelY = G['p'][:,1][Gclose][Ghalo] - haloY
GprelZ = G['p'][:,2][Gclose][Ghalo] - haloZ	
GvrelX = G['v'][:,0][Gclose][Ghalo] * np.sqrt(a) - haloVX + GprelX*a*Hubble/(h*1000)
GvrelY = G['v'][:,1][Gclose][Ghalo] * np.sqrt(a) - haloVY + GprelY*a*Hubble/(h*1000)
GvrelZ = G['v'][:,2][Gclose][Ghalo] * np.sqrt(a) - haloVZ +  GprelZ*a*Hubble/(h*1000)
#rescale positions so that magnitude = 1 (unit vectors)
GprelX_RS = GprelX / dists[Ghalo]
GprelY_RS = GprelY / dists[Ghalo]
GprelZ_RS = GprelZ / dists[Ghalo]
v_rad = GprelX_RS*GvrelX  + GprelY_RS*GvrelY + GprelZ_RS*GvrelZ
print 'vrad ',v_rad
#cut particles that are "outflowing"
OutflowCut = v_rad > Vsig/np.sqrt(3)  
if (use_v0):
	OutflowCut = v_rad > 0

OutflowingParticlesID_G2 = np.copy(G['id'][Gclose][Ghalo][OutflowCut])
OutflowingParticleMass_G2 = np.copy(G['m'][Gclose][Ghalo][OutflowCut])
OutflowingParticle_vrad_G2 = np.copy(v_rad[OutflowCut])

OutflowingParticlesID_G2_095Rvir = np.copy(G['id'][Gclose][Ghalo][OutflowCut*ShellCuts[index2]])
OutflowingParticleMass_G2_095Rvir = np.copy(G['m'][Gclose][Ghalo][OutflowCut*ShellCuts[index2]])
OutflowingParticle_vrad_G2_095Rvir = np.copy(v_rad[OutflowCut*ShellCuts[index2]])
Particles_G2_095_len = len(G['m'][Gclose][Ghalo][ShellCuts[index2]])

count = 0
ShellCutsInner = np.copy(ShellCuts[index2])
while (count < index2):
	ShellCutsInner += ShellCuts[count]
	count+=1 


#Particles_Within_095 = np.copy(G['id'][Gclose][Ghalo][ShellCuts[index2]])

 
Particles_G2_095 = np.copy(G['id'][Gclose][Ghalo][ShellCuts[index2]])

print 'length of 0.95 outflowing particles',len(OutflowingParticlesID_G2_095Rvir)
len95  = len(OutflowingParticlesID_G2_095Rvir)

second_outflow_rate = OR.outflow_flux_atshell(G['m'][Gclose][Ghalo], v_rad, OutflowCut, ShellCuts[index2], ShellThickness_ar[index2], h)


#fraction of outflowing particles at second epoch at shell #2 that are from first burst at epoch #1
here =  np.in1d(OutflowingParticlesID_G2_095Rvir, OutflowingParticlesID_025Rvir)

print 'total number of outflowing particles found at 0.95 Rvir at new snapshot that were outflowing at 0.25 Rvir in old snapshot ',len(OutflowingParticlesID_G2_095Rvir[here])

len95in25 =  len(OutflowingParticlesID_G2_095Rvir[here])

Particles_G_095_len = len(Particles_G2_095)

here_still  =  np.in1d(Particles_G2_095, OutflowingParticlesID_025Rvir)

#fraction of particles at shell #2 at later epoch still from first burst of shell#1
len95stillin25 =  len(Particles_G2_095[here_still])


#1. halo, 2. Nsnap1, 3. Nsnap2, 4. redshift1, 5. redshift2, 6. len25-number of outflowing particles from epoch1, shell1. 7. len95 - number of outflowing particles at epoch2, shell2. 8. len95in25 number of outflowing particles at epoch2, shell2, that were also outflowing at shell1 at epoch1. 9. len95in25/len95.  Fraction of outflowing particles at shell2 epoch2 that were also outflowing in shell1 at epoch1 over the number of outflowing particles in shell2 epoch2. 10. len95in25/len25. Fraction of outflowing particles at shell2 epoch2 that were also outflowing in shell1 at epoch1, over the number of outflowing particles from shell1 epoch1. 11. len95stillin25)/float(len25)- total number of particles in epoch2 shell2 that originated from outflowing particles from epoch1 shell 1. 12. index1. 13. index2. 14. use_v0 (True if vcut =0, False if vcut=vsig) #15 outflowrate1 #16 outflowrate2
print halo_to_do, Nsnap, Nsnap2, redshift1, redshift2, len25, len95, len95in25, len95stillin25, (float(len95in25)/float(len95)), (float(len95in25)/float(len25)), (float(len95stillin25)/float(len25)), index1, index2, use_v0
BigString = [halo_to_do, Nsnap, Nsnap2, redshift1, redshift2, len25, len95, len95in25, len95stillin25, (float(len95in25)/float(len95)), (float(len95in25)/float(len25)), (float(len95stillin25)/float(len25)), index1, index2, int(use_v0), original_outflow_rate, second_outflow_rate ]
thestring = SF.line_to_string(BigString)

foutname = 'entrain_outflow_'+str(Nsnap)+'halo'+str(halo_to_do)+'.txt'
g = open(foutname, 'a')
g.write(thestring)
g.close()

