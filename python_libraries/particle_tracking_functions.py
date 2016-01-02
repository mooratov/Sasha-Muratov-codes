#functions designed to be used in particle tracking 

import os.path
import numpy as np
from readsnap import readsnap
from Sasha_functions import *
import sys
#def write_outflows(Nsnap, Gid, 

def find_new_outflows(Gid, halocount, Nsnap_prev, Inflows = 'n', particular_shell = -999):
	if (Nsnap_prev < 10):
		Nsnap_prev_string = '00'+str(Nsnap_prev)
	elif (Nsnap_prev < 100):
		Nsnap_prev_string = '0'+str(Nsnap_prev)
	else:
		Nsnap_prev_string = str(Nsnap_prev)
	if (Inflows == 'y'):
		finname = 'Inflowing_Particles'+Nsnap_prev_string+'N'+str(halocount)+'.txt'
	else:
		finname = 'Outflowing_Particles'+Nsnap_prev_string+'N'+str(halocount)+'.txt'
	if not(os.path.isfile(finname)):
		NewOutflows = np.in1d(Gid, Gid)
		print 'no previous outflow/inflow file found, assuming that this is the first one?'
		return NewOutflows
	T = np.loadtxt(finname, ndmin=2)
	if (len(T) < 1):
		if (len(Gid) > 0):
			print 'all outflows are new outflows!'
			NewOutflows = np.in1d(Gid, Gid)
			return NewOutflows
		else:
			print 'no outflows before, and none now'
			NewOutflows = np.array([])
			return NewOutflows
			
	Gid_prev = T[:,0]
		
	if(particular_shell>-1):
		print 'using shell cut ',particular_shell
		Shell_prev = T[:,9]
		ShellCut = (Shell_prev == particular_shell)	
		print len( Gid_prev[ShellCut]),' out of ',len(Gid_prev),' qualify'
		OldOutflows = np.in1d(Gid, Gid_prev[ShellCut])
	else:
		OldOutflows = np.in1d(Gid, Gid_prev)
		
	NewOutflows = np.invert(OldOutflows)
	print 'found some new outflows'
	return NewOutflows
	
def match_particles(PID_list, Nsnap, use_multi='n'):

	if (int(Nsnap) < 10):
		Nsnapstring = '00'+str(Nsnap)
	elif (int(Nsnap) < 100):
		Nsnapstring = '0'+str(Nsnap)
	else:
		Nsnapstring = str(Nsnap)

	the_snapdir = './hdf5/'
	the_prefix ='snapshot'
	if (use_multi == 'y'):
		the_snapdir = './snapdir_'+Nsnapstring+'/'
	the_suffix = '.hdf5'
	
	print 'reading particle type 0 (gas?)' 
	G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)
	here =  np.in1d(G['id'], PID_list)
	
	PID = G['id'][here]
	
	x = G['p'][:,0][here]
	y = G['p'][:,1][here]
	z = G['p'][:,2][here]

	vx = G['v'][:,0][here]
	vy = G['v'][:,1][here]
	vz = G['v'][:,2][here]
	h = G['header'][12]
	PhysTemp, PhysRho = convertTemp(G['u'][here], G['ne'][here], G['rho'][here], h)
	Metallicity = G['z'][:,0][here]
	a =  G['header'][2]
	boxsize = G['header'][9]
	tbl = []
	base = [PID, x,y,z, vx,vy,vz, PhysTemp, PhysRho, Metallicity]
	for part in base:
		tbl.append(part)
	tbl = np.array(tbl)
	Ttbl = tbl.T
	return (Ttbl, a, boxsize)
	
def trace_particles_forward(PID_list, Nsnap, haloN, Nforward, use_multi='n'):
	part_table = []
	#part_table.append (np.sort(PID_list))
	count = 0
	headerstring = []
	while (count <= Nforward):
		print 'going forth ',count,' steps'
		N_to_do = Nsnap + count
		TTbl, a, boxsize = match_particles(PID_list, N_to_do, use_multi)
		redshift = 1.0/a - 1.0
		redshiftstring = "{0:.3f}".format(redshift)
		halostats = find_halo_now(haloN, a, redshiftstring)
		print 'halo stats ',halostats

		PID_current = TTbl[:,0]
		if (len(PID_current) != len(PID_list)):
			print 'oh boy some turned into stars! ',len(PID_current),' ',len(PID_list)
			found = np.in1d(PID_list, PID_current)
			lost = np.invert(found)
			for lost_boy in PID_list[lost]:
				print 'adding in lost boy ',lost_boy
				PID_current = np.append(PID_current, lost_boy)
				base = [lost_boy, 999999.0,999999.0,999999.0, 0.0, 0.0, 0.0, -1, -1, -1]
				base = np.array(base)
				print base.shape
				print TTbl.shape
				#TTbl = np.column_stack([TTbl, base])
				TTbl = np.vstack([TTbl, base])
		haloX = halostats[2]
		haloY = halostats[3]
		haloZ = halostats[4]
		Rvir = halostats[11]
		headerstring.append(Rvir)
		partX = TTbl[:,1]
		partY = TTbl[:,2]
		partZ = TTbl[:,3]
		dists = calcdist2(haloX, haloY, haloZ, partX, partY, partZ, boxsize)
		Sortorder = np.argsort(PID_current)
		print 'for count',count,' appending ', dists[Sortorder]
		part_table.append(PID_current[Sortorder])
		part_table.append(dists[Sortorder])
		count+=1
	part_table = np.array(part_table)
	headerstring = np.array(headerstring)
	return (part_table.T, headerstring)
	
	#for starters lets make it simple - print a data table with particle distances
def trace_particles_backwards(PID_list, Nsnap, haloN, Nback, use_multi='n'):
	part_table = []
	#part_table.append (np.sort(PID_list))
	count = 0
	headerstring = []
	while (count <= Nback):
		print 'going back ',count,' steps'
		N_to_do = Nsnap - count
		TTbl, a, boxsize = match_particles(PID_list, N_to_do, use_multi)
		redshift = 1.0/a - 1.0
		redshiftstring = "{0:.3f}".format(redshift)
		halostats = find_halo_now(haloN, a, redshiftstring)
		print 'halo stats ',halostats
		haloX = halostats[2]
		haloY = halostats[3]
		haloZ = halostats[4]
		Rvir = halostats[11]
		headerstring.append(Rvir)
		partX = TTbl[:,1]
		partY = TTbl[:,2]
		partZ = TTbl[:,3]
		PID_current = TTbl[:,0]
		dists = calcdist2(haloX, haloY, haloZ, partX, partY, partZ, boxsize)
		Sortorder = np.argsort(PID_current)
		print 'for count',count,'  ',N_to_do,' appending ', dists[Sortorder]
		part_table.append(PID_current[Sortorder])
		part_table.append(dists[Sortorder])
		count+=1
	part_table = np.array(part_table)
	headerstring = np.array(headerstring)
	return (part_table.T, headerstring)

def simple_previous_snap_read(Nsnap, use_multi='n'):
	if (int(Nsnap) < 10):
		Nsnapstring = '00'+str(Nsnap)
	elif (int(Nsnap) < 100):
		Nsnapstring = '0'+str(Nsnap)
	else:
		Nsnapstring = str(Nsnap)

	the_snapdir = './hdf5/'
	the_prefix ='snapshot'
	if (use_multi == 'y'):
		the_snapdir = './snapdir_'+Nsnapstring+'/'
	the_suffix = '.hdf5'
	
	print 'reading particle type 0 (gas?)' 
	G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)
	if (len(G['p'][:,0]) < 1 ):
		print 'no previous snapshot detected! ',Nsnap
		sys.exit()
	x = G['p'][:,0]
	y = G['p'][:,1]
	z = G['p'][:,2]
	ID = G['id']
	a = G['header'][2]
	Ttbl = (np.array([ID, x, y, z])).T
	return Ttbl, a
