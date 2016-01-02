#calculates explicitly escape velocity from different positions in simulated halo. 

import matplotlib
matplotlib.use('Agg')
import numpy as np
import Sasha_functions as SF
from readsnap import readsnap
import neo_potential_and_sigmaSFR as potter
#import potential_and_sigmaSFR as potter2
import sys
from tasz import hubble_param

from datetime import date
today = date.today()

do_graphs = True
if (do_graphs):
	import matplotlib.pyplot as plt
	import graphics_library as GL


multifile = False
use_no_Pep = True
use_KDE_cent = False
use_darkKDE_cent = False
adaptive_Rstar = False
use_fixed_halos = 4
the_surface_frac = 0.8
halo_to_do = [0]
pc = 3.08567758e18
Inner_SF_thresh = 0.2
where_do_we_go = 25.0 #distance that stuff needs to propagate to 

if (len(sys.argv) < 2):
	print 'syntax: blah.py Nsnapshot'
	sys.exit()
	
###########
# I/O 
###########

if (int(sys.argv[1]) < 10):
	Nsnapstring = '00'+str(sys.argv[1])
elif (int(sys.argv[1]) < 100):
	Nsnapstring = '0'+str(sys.argv[1])
else:
	Nsnapstring = str(sys.argv[1])
Nsnap = int(sys.argv[1])

#for multiple snapshots

#for HDF5
the_snapdir = './hdf5/'
if (multifile):
	the_snapdir = './snapdir_'+Nsnapstring+'/'
the_prefix ='snapshot'
the_suffix = '.hdf5'

G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)
S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
D = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix)

thetime = S['header'][2]
redshift = S['header'][3]
boxsize = S['header'][9]
omega_matter = S['header'][10]
omega_L = S['header'][11]
h = S['header'][12]

a = float(thetime)
print 'using a = ',a
theredshift = (1.0/a - 1)
redshiftstring = "{0:.3f}".format(theredshift)
print 'z = ',"{0:.3f}".format(theredshift)

SFR_time_range, age_in_Gyr = SF.cut_stellar_age(a, Nsnap, S['age'], omega_matter, h)
stellar_age_cut =  age_in_Gyr < SFR_time_range/1e9
Syoung = stellar_age_cut
Rstars = where_do_we_go * h / a
RstarsPhys = where_do_we_go



workcount = 0 
Hubble = hubble_param(a, omega_matter, h)


###########
#Work function
###########
C = SF.read_halo_catalog(Nsnapstring, redshiftstring)
while (workcount < len(halo_to_do)):
	halostats = SF.find_halo_now(halo_to_do[workcount], a, therod=use_fixed_halos)
	haloN = halostats[1]
	Rvir = halostats[11]
	haloX = halostats[2]
	haloY = halostats[3]
	haloZ = halostats[4]
	Mstar = halostats[13]
	Vsig = halostats[10]
	Mgas = halostats[12]
	M = halostats[4]

	Rvir, Vsig, M = SF.check_Rvir_growth(halo_to_do[workcount], a, Rvir, Vsig, M, therod=use_fixed_halos)

	if (adaptive_Rstar):
		Rstars = Inner_SF_thresh*Rvir
		RstarsPhys = Rstars * a / h

	
	
	
	Sdists = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize)
	Sinner_young =  Sdists[Syoung] < Rstars

	if (use_KDE_cent):
		Sinner = Sdists<Rstars
		S_COM = SF.gaussian_KDE_center(S['p'][:,0][Sinner], S['p'][:,1][Sinner], S['p'][:,2][Sinner], downsample=True)
		[haloX, haloY, haloZ] = [S_COM[0], S_COM[1], S_COM[2]]
		Sdists = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize)
		Sinner_young =  Sdists[Syoung] < Rstars		
	elif (use_darkKDE_cent):
		#dark matter fun
		Ddists = SF.calcdist2(haloX, haloY, haloZ, D['p'][:,0], D['p'][:,1], D['p'][:,2], boxsize)	
		Dinner = Ddists<Rstars
		D_COM = SF.gaussian_KDE_center(D['p'][:,0][Dinner], D['p'][:,1][Dinner], D['p'][:,2][Dinner], downsample=True)
		[haloX, haloY, haloZ] = [D_COM[0], D_COM[1], D_COM[2]]
		Sdists = SF.calcdist2(haloX, haloY, haloZ, S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize)
		Sinner_young =  Sdists[Syoung] < Rstars

	
	[x_cells, y_cells, z_cells, cellmass ]= potter.surface_density_bycells(S, Syoung, Sinner_young, a, h, surface_frac=the_surface_frac)
	Nnewstars = len(Sdists[Syoung][Sinner_young])
	print 'total N newly formed stellar particles in Inner region', Nnewstars
	Mnewstars = np.sum(S['m'][Syoung][Sinner_young])*1e10/h

	print 'total Mass of newly formed stars ',Mnewstars
	
	print 'time since last snap ',SFR_time_range
	

	
	count = 0 
	vesc_cge_ar = []
	count_ar = []
	frac_ar = []
	cellmass*=1e10/h
	totalfrac = 0 

	vesc_cge_ar = potter.calculate_escape_velocity([haloX, haloY, haloZ], [x_cells, y_cells, z_cells], G, S, D, where_do_we_go*h/a)  * np.sqrt(2.0)
	print 'halo is here ',haloX, haloY, haloZ
	while (count < len(x_cells)):
		totalfrac+=cellmass[count]/Mnewstars
		#vesc_cge_beta = potter2.calculate_escape_velocity([haloX, haloY, haloZ], [x_cells[count], y_cells[count], z_cells[count]], G, S, D,   where_do_we_go*h/a)*np.sqrt(2.0)
		print count, x_cells[count], y_cells[count], z_cells[count], vesc_cge_ar[count], cellmass[count],  cellmass[count]/Mnewstars, totalfrac
		frac_ar.append(cellmass[count]/Mnewstars)
		count_ar.append(count)
		count +=1

	vesc_cge_ar = np.array(vesc_cge_ar)
	count_ar = np.array(count_ar)
	frac_ar = np.array(frac_ar)
	vhist_header = [Nsnap, redshift, where_do_we_go ,  Rstars*a/h, Mnewstars, SFR_time_range, Nnewstars, the_surface_frac, haloX, haloY, haloZ]
	vhist_header_string = SF.line_to_string(vhist_header, newline='n')

	outtable = np.column_stack([count_ar, cellmass, vesc_cge_ar, x_cells, y_cells, z_cells, frac_ar ])
	foutname = str(today)+'Xsurface_list'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.txt'
	np.savetxt(foutname, outtable,  fmt='%1.6e', header=vhist_header_string)
	if (do_graphs):
		foutname1a = str(today)+'image_stars'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		foutnameDM = str(today)+'image_DM'+Nsnapstring+'N'+str(halo_to_do[workcount])+'.pdf'
		GL.pluscuts_2d_desnity_hist(S, Rstars, haloX, haloY, haloZ,  foutname1a, dolabel=True, extrapoints = [[x_cells[0], y_cells[0], z_cells[0]]], extrars=[-1])
		GL.pluscuts_2d_desnity_hist(D, Rstars, haloX, haloY, haloZ, foutnameDM, dolabel=True)
	workcount+=1
