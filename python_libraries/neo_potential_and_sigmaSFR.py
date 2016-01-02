#functions used to calculate potential differences and star formation surface density 
import numpy as np
import Sasha_functions as SF
import scipy
import matplotlib.pyplot as plt



def calculate_escape_velocity(halo_pos, particle_pos, G, S, D, rmax):
	newtonG = 6.67384e-8
	pc = 3.08567758e18
	kpc = pc * 1e3
	a = G['header'][2]
	h = G['header'][12]
	boxsize = G['header'][9]
	UnitMass_in_g=1.989e43 / h
	do_plots = False

	Gdists = SF.calcdist2(halo_pos[0], halo_pos[1], halo_pos[2], G['p'][:,0], G['p'][:,1], G['p'][:,2], boxsize)
	Sdists =  SF.calcdist2(halo_pos[0], halo_pos[1], halo_pos[2], S['p'][:,0], S['p'][:,1], S['p'][:,2], boxsize)
	Ddists =  SF.calcdist2(halo_pos[0], halo_pos[1], halo_pos[2], D['p'][:,0], D['p'][:,1], D['p'][:,2], boxsize)
	
	
	particle_dists = SF.calcdist2(halo_pos[0], halo_pos[1], halo_pos[2], particle_pos[0], particle_pos[1], particle_pos[2], boxsize)
	
	#print 'particle dist ',particle_dist

	Gcut = Gdists < rmax
	Scut = Sdists < rmax
	Dcut = Ddists < rmax
	
	Gdist_ord = np.argsort(Gdists[Gcut])
	Ddist_ord = np.argsort(Ddists[Dcut])
	Sdist_ord = np.argsort(Sdists[Scut])
	
	
	#print Gdists[Gcut][Gdist_ord]
	
	GMenc_vsr =  np.cumsum(G['m'][Gcut][Gdist_ord]) 
	SMenc_vsr =  np.cumsum(S['m'][Scut][Sdist_ord]) 
	DMenc_vsr =  np.cumsum(D['m'][Dcut][Ddist_ord]) 




#	GM_for_r = interp1d(Gdists[Gcut][Gdist_ord], GMenc_vsr, kind='linear', bounds_error=False)
#	SM_for_r = interp1d(Sdists[Scut][Sdist_ord], SMenc_vsr, kind='linear', bounds_error=False)
#	DM_for_r = interp1d(Ddists[Dcut][Ddist_ord], DMenc_vsr, kind='linear', bounds_error=False)

	rspace = np.linspace(min(particle_dists)/200.0, rmax, num=1000000)
	GM_for_r = np.interp(rspace, Gdists[Gcut][Gdist_ord], GMenc_vsr)
	SM_for_r = np.interp(rspace, Sdists[Scut][Sdist_ord], SMenc_vsr)
	DM_for_r = np.interp(rspace, Ddists[Dcut][Ddist_ord], DMenc_vsr)

	#print 'rspace ',rspace
	#print 'Ddists[Dcut][Ddist_ord] ',Ddists[Dcut][Ddist_ord]

	rspace_physical = rspace*a / h
	theM = GM_for_r + SM_for_r + DM_for_r
	
	
	initpot = (newtonG * theM[0])/rspace_physical[0]
	
	dpot_array = newtonG * theM / rspace_physical**2.0

	if (do_plots):	
		fig1 =  plt.figure(figsize=(10,9))
		plt.plot(rspace, theM, ':r')
		plt.savefig('alpha_test_mass.pdf')
	
		fig2 =  plt.figure(figsize=(10,9))
		plt.plot(rspace, dpot_array, ':k')
		plt.savefig('alpha_test_pot.pdf')
	
	
	dr = (rspace_physical[-1] - rspace_physical[0])/1000000.0
	
	#thepot_ar = scipy.integrate.
	
	thepot_ar = scipy.integrate.cumtrapz(dpot_array, x=rspace_physical, dx=dr, initial=initpot)  #this is 0 to r 
	
	thepot_cge_ar = thepot_ar * UnitMass_in_g / kpc
	
	mypot = thepot_cge_ar[-1] - np.interp(particle_dists, rspace, thepot_cge_ar)
	
	testpot = np.interp(1.0, rspace, thepot_cge_ar)
	print 'test pot alpha ',testpot
	
	vesc_cge = np.sqrt(mypot) / 1e5 #convert to km/s
	return  vesc_cge

	
def surface_density_bycells(S, agecut, poscut, a, h, cell_length=0.3, surface_frac=0.8, density_max_cell_center=True):
	#r = np.random.randn(100,3)
	
	
	r = [S['p'][:,0][agecut][poscut], S['p'][:,1][agecut][poscut], S['p'][:,2][agecut][poscut]]
	
	minx = np.min(S['p'][:,0][agecut][poscut])*a/h
	maxx = np.max(S['p'][:,0][agecut][poscut])*a/h
	
	miny = np.min(S['p'][:,1][agecut][poscut])*a/h
	maxy = np.max(S['p'][:,1][agecut][poscut])*a/h
	
	minz = np.min(S['p'][:,2][agecut][poscut])*a/h
	maxz = np.max(S['p'][:,2][agecut][poscut])*a/h
	
	numxcells = int((maxx-minx)/cell_length)+1
	numycells = int((maxy-miny)/cell_length)+1
	numzcells = int((maxz-minz)/cell_length)+1
	
	xrange = [minx*h/a, (minx*h/a + float(numxcells)*cell_length*h/a)]
	yrange = [miny*h/a, (miny*h/a + float(numycells)*cell_length*h/a)]
	zrange = [minz*h/a, (minz*h/a + float(numzcells)*cell_length*h/a)]
	
#	print 'diagnostic x',(xrange[1] - xrange[0])/(cell_length*h/a), numxcells
#	print 'diagnostic y',(yrange[1] - yrange[0])/(cell_length*h/a), numycells
#	print 'diagnostic x',(zrange[1] - zrange[0])/(cell_length*h/a), numzcells

	#rad should be in physical kpc
	
	
	#num_cells=int(rad*2.0/cell_length)
	H, edges = np.histogramdd(r, bins = (numxcells, numycells,  numzcells), range=(xrange,yrange,zrange), weights=S['m'][agecut][poscut])

	cut = H>0
	ix, iy, iz = np.where(cut)

	count = 0
	Nums = []
	cenx = []
	ceny = []
	cenz = [] 

	#print 'iy ', iy
	#print 'len iy ',len(iy)
	#print 'edges ',edges
	#print 'len edges ',len(edges)
    
	while (count < len(ix)):
		theNum = H[ix[count],iy[count],iz[count]]
		thecenx = (edges[0][ix[count]] + edges[0][ix[count]+1])/2.0
		theceny = (edges[1][iy[count]] + edges[1][iy[count]+1])/2.0
		thecenz = (edges[2][iz[count]] + edges[2][iz[count]+1])/2.0		
		if (density_max_cell_center):
			cellxmin_cut = S['p'][:,0][agecut][poscut] > edges[0][ix[count]]
			cellxmax_cut = S['p'][:,0][agecut][poscut] <  edges[0][ix[count]+1]
			cellymin_cut = S['p'][:,1][agecut][poscut] > edges[1][iy[count]]
			cellymax_cut = S['p'][:,1][agecut][poscut] <  edges[1][iy[count]+1]
			cellzmin_cut = S['p'][:,2][agecut][poscut] > edges[2][iz[count]]
			cellzmax_cut = S['p'][:,2][agecut][poscut] < edges[2][iz[count]+1]
			
			allcut = cellxmin_cut* cellxmax_cut* cellymin_cut* cellymax_cut * cellzmin_cut * cellzmax_cut
			#print 'test test ',len(S['p'][:,0][agecut][poscut][allcut]), S['p'][:,0][agecut][poscut][allcut]

			if (len(S['p'][:,0][agecut][poscut][allcut]) > 5):
				#print 'first center ',[thecenx, theceny, thecenz], theNum
				S_COM = SF.gaussian_KDE_center(S['p'][:,0][agecut][poscut][allcut], S['p'][:,1][agecut][poscut][allcut] , S['p'][:,2][agecut][poscut][allcut] , downsample=False)
				[thecenx, theceny, thecenz] = S_COM
				#print 'adjusting center ',[thecenx, theceny, thecenz]
			
		Nums.append(theNum)
		cenx.append(thecenx)
		ceny.append(theceny)
		cenz.append(thecenz)
		count+=1
	cenx = np.array(cenx)
	ceny = np.array(ceny)
	cenz = np.array(cenz)
	Nums = np.array(Nums)

	order_of_densities = np.argsort(Nums)[::-1]
	#cumsumdens = np.cumsum(Nums[order_of_densities])

	sum_of_densities = np.sum(Nums)

	current_sum = 0
	count = 0 

	final_x = []
	final_y = []
	final_z  = []
	final_Nums = []
	
	#print 'densities in order ',Nums[order_of_densities]

	while (current_sum <= surface_frac * sum_of_densities):
		current_sum += Nums[order_of_densities][count]
		final_x.append(cenx[order_of_densities][count])
		final_y.append(ceny[order_of_densities][count])
		final_z.append(cenz[order_of_densities][count])
		final_Nums.append( Nums[order_of_densities][count])
		count+=1


	return [final_x, final_y, final_z, final_Nums]
	#star_forming_volume = len(H_cut) * volume_of_cell


	#can also include a thing to find maximum density
def surface_density_byclump(S):
	return 0
	#do this histogram 
	
	#find the center , then find the maximum bin near the center
	#traverse to rhalf to see that it stays unclummpy
	#if unclumpy, fit sersic index thing
	#if clumpy, fit out to where it stopsclumpy
	#work iteratively to find second clump
	
     

		
