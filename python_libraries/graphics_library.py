import matplotlib
matplotlib.use('Agg')
from gas_temperature import *
from matplotlib.colors import LogNorm
#from pylab import *
from Sasha_functions import sasha_min
from Sasha_functions import sasha_max
matplotlib.use('Agg')
import matplotlib.pyplot as theplt
matplotlib.use('Agg')
import matplotlib.cm as cm

#simple 2d histogram with optional annulus but with the annulus thing. How about i just write a separate function to draw the annulus - done!
#now I've added cuts in as well

def pluscuts_2d_desnity_hist(P, Rbox, x, y, z, foutname, Cuts=[], focusshell=-1, thecmap = 'CMRmap', numbins=300, thevmin=1, thevmax=1e3, PhysWeight = [], dolabel=False, doylabelL=False, xtoz=False, ytoz=False, extrapoints = [], extrars=[], adaptivemax=False, Ncol_label=False):
	print 'lets graph! ',foutname

	little_h = 0.7

	xx = P['p'][:,0]/little_h
	yy = P['p'][:,1]/little_h
	zz = P['p'][:,2]/little_h

	xx -= x/little_h
	yy -= y/little_h
	zz -= z/little_h


	xmin = 0 - Rbox/little_h
	ymin = 0 - Rbox/little_h
	zmin = 0 - Rbox/little_h
                                                                    
	xmax = 0 + Rbox/little_h
	ymax = 0 + Rbox/little_h
	zmax = 0 + Rbox/little_h
	
	firstmin = xmin
	secondmin = ymin
	firstmax = xmax
	secondmax = ymax
	
	firstone = xx
	secondone = yy
	
	if (xtoz):
		firstone = zz
		firstmin = zmin
		firstmax = zmax
	elif (ytoz):
		secondone = zz
		secondmin = zmin
		secondmax = zmax
	
	theplt.figure(figsize=(11,8))
	if (adaptivemax):
		if ( np.sqrt(float(len(firstone))) < thevmax):
			thevmax = np.sqrt(float(len(firstone)))
	if ((len(Cuts)>1) and (len(PhysWeight) < 1)):
		#just cuts, no weights 
		plt = theplt.hist2d(firstone[Cuts], secondone[Cuts], range=[[firstmin,firstmax],[secondmin,secondmax]], bins=numbins, cmap=cm.get_cmap(thecmap), norm=LogNorm(), vmin=thevmin, vmax=thevmax)
	elif ((len(Cuts)>1) and (len(PhysWeight) > 1)):
		#both cuts and weights
		plt = theplt.hist2d(firstone[Cuts], secondone[Cuts], range=[[firstmin,firstmax],[secondmin,secondmax]], bins=numbins, cmap=cm.get_cmap(thecmap), norm=LogNorm(), weights=PhysWeight[Cuts], vmin=thevmin, vmax=thevmax)
	else:
		#nothin but the points.
		plt = theplt.hist2d(firstone, secondone, range=[[firstmin,firstmax],[secondmin,secondmax]], bins=numbins, cmap=cm.get_cmap(thecmap), norm=LogNorm(), vmin=thevmin, vmax=thevmax)


	for f in plt[0]:
		f += 1e-1 #set the background value for smooth image.
	firstaxis = 'x'
	secondaxis = 'y'
	if (xtoz):
		firstaxis = 'z'
	if (ytoz):
		secondaxis = 'z'
	theplt.xlabel(firstaxis+' (comoving kpc)', fontsize=25) 
	if (doylabelL):
		theplt.ylabel(secondaxis+' (comoving kpc)', fontsize=25) 
	cbar = theplt.colorbar()
	ax = theplt.gca()
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()
	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(25)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(25)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)	
	
	cbartick = cbar.ax.get_yticklabels()
	for label in cbartick:
		label.set_fontsize(28)
		label.set_family('serif')	
	cbar.ax.yaxis.set_tick_params(width=2.5, length=8)

	if (len(extrapoints)>1):
		marker_ar = ['*', 'p', 's', '^', 'v']*10
		color_ar = ['chartreuse','chartreuse','chartreuse','chartreuse','chartreuse']*10
		count = 0
		for pt in extrapoints:
			axis1 = 0 
			axis2 = 1
			tempx = x
			tempy = y 
			if (xtoz):
				axis1 = 2 
				tempx = z
			if (ytoz):
				axis2 = 2 
				tempy = z
			theplt.plot([(pt[axis1] - tempx)/little_h], [(pt[axis2]- tempy)/little_h], marker=marker_ar[count], color=color_ar[count], ms=5, alpha=0.1, fillstyle='none')
			print 'plotted ',(pt[axis1] - tempx)/little_h, (pt[axis2]- tempy)/little_h
			print 'original ',pt
			ang = 0
			if (len(extrars) > 0):
				if (extrars[count] > 0):
					ra = extrars[count] / little_h
					rb = ra
					theX,theY=ellipse(ra,rb,ang,(pt[axis1] - tempx)/little_h,(pt[axis2]- tempy)/little_h)	
					theplt.plot(theX,theY,ls=':', color=color_ar[count] ,ms=1,linewidth=1)

			count += 1
	
	if (dolabel):
		if (Ncol_label):
			cbar.ax.set_ylabel(r'$N ({\rm cm}^{-2})$',fontsize=28)
		else:
			cbar.ax.set_ylabel(r'$\Sigma (M_{\odot} {\rm kpc}^{-2})$',fontsize=28)
	if (focusshell > -1):
		triple_annulus(focusshell, (Rbox/little_h), 0, 0)

	theplt.savefig(foutname)
	theplt.clf()


#lol here is the ellipse function. Copied from somewhere. 
def ellipse(ra,rb,ang,x0,y0,Nb=100):
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang
    co,si=np.cos(an),np.sin(an)
    the=np.linspace(0,2*np.pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y

#just draws an annulus between focusshell and focusshell+1
def draw_the_annulus(focusshell, Rbox, x, y, onlyone=True):
		ang = 0
		ra = (focusshell/10.0)*Rbox
		rb = ra 
	
		ra2 = ((focusshell+1.0)/10.0)*Rbox
		rb2 = ra2
	
		theX,theY=ellipse(ra,rb,ang,x,y)
		theplt.plot(theX,theY,"w:",ms=1,linewidth=3)
			
		if (not onlyone):
			theX,theY=ellipse(ra2,rb2,ang,x,y)
			theplt.plot(theX,theY,"w:",ms=1,linewidth=3)	
			
def triple_annulus(focusshell, Rbox, x, y, onlyone=True):
		ang = 0
		ra = (10.0/10.0)*Rbox
		rb = ra 
		
		theX,theY=ellipse(ra,rb,ang,x,y)
		theplt.plot(theX,theY,"w",ms=1,linewidth=5)			

		ang = 0
		ra = (2.5/10.0)*Rbox
		rb = ra 
		
		theX,theY=ellipse(ra,rb,ang,x,y)
		theplt.plot(theX,theY,"w:",ms=1,linewidth=5)		
	
		ang = 0
		ra = (5.0/10.0)*Rbox
		rb = ra 
		
		theX,theY=ellipse(ra,rb,ang,x,y)
		theplt.plot(theX,theY,"w:",ms=1,linewidth=5)						
    
def phase_diagram(Temp, n, M, foutname, Cuts, zstring='-1', numbins=200, thecmap='hot',  dolabels=False, xmin = -6, ymin = 1, xmax = 3, ymax = 8):
	print 'lets graph! ',foutname
	n = np.log10(n)
	Temp = np.log10(Temp)
	

	
	plt = theplt.hist2d(n[Cuts], Temp[Cuts], range=[[xmin,xmax],[ymin,ymax]], bins=numbins, cmap=cm.get_cmap(thecmap), norm=LogNorm(), weights=M[Cuts], vmin=1e-6, vmax=5e-2)
	


	TempCut = Temp > 5 

	theplt.xlabel(r'log n $({\rm cm}^{-3})$') 
	theplt.ylabel(r'log T (K)') 
	
	if (dolabels )	:
		title('total mass '+"{0:.3e}".format(np.sum(M[Cuts])))
	
		text(0, 6, " T > 5: "+"{0:.3e}".format(np.sum(M[Cuts*TempCut]))+'\n z='+zstring, color='w')

	for f in plt[0]:
		f += 1e-10
	cbar = theplt.colorbar()
	cbar.ax.set_ylabel(r'Mass in pixel $(10^{10} M_{\odot}/h)$')

	theplt.savefig(foutname)
	theplt.clf()

def phase_histogram(Temp, Cuts, foutname, zstring='-1', numbins=200):
	print 'lets graph! ',foutname
	Temp = np.log10(Temp)
	Tmin_lim = 0 
	Tmax_lim = 10
	#n, bins, patches = plt.hist(Temp[Cuts], num_bins, range=(Tmin_lim, Tmax_lim), weights=weightarray , histtype='step',color=thecolors[bincount-1], alpha=0.8)
	if (len(Temp[Cuts]) > 0):
		n, bins, patches = theplt.hist(Temp[Cuts], numbins, range=(Tmin_lim, Tmax_lim), histtype='step', color='black', alpha=0.8)
		theplt.title(zstring)
		theplt.xlabel('log T $(K)$') 
		theplt.yscale('log')
		theplt.ylabel('log N particles')
		theplt.savefig(foutname)
		theplt.clf()
	else: 
		print 'error! there are no particles '

def phase_histogram_delux(Temp, Cuts, OutCuts, foutname, zstring='-1', numbins=200):
	print 'lets graph! ',foutname
	Temp = np.log10(Temp)
	Tmin_lim = 0 
	Tmax_lim = 10
	#n, bins, patches = plt.hist(Temp[Cuts], num_bins, range=(Tmin_lim, Tmax_lim), weights=weightarray , histtype='step',color=thecolors[bincount-1], alpha=0.8)
	if (len(Temp[Cuts]) > 0):
		n, bins, patches = theplt.hist(Temp[Cuts], numbins, range=(Tmin_lim, Tmax_lim), histtype='step', color='black', log=True)
		if (len(Temp[OutCuts]) >0 ):
			n2, bins2, patches2 = theplt.hist(Temp[OutCuts], numbins, range=(Tmin_lim, Tmax_lim), histtype='bar', color='red', alpha=0.8, log=True, rwidth=1.0)
			#n2+=1e-10
		else:
			print 'no outflows'
		theplt.title('z ='+zstring)
		theplt.xlabel('log T $(K)$')
		theymax = sasha_max([1e5, sasha_max(n)])
		theplt.ylim(1, theymax) 
		theplt.yscale('log')
		theplt.ylabel('log N particles')
		theplt.savefig(foutname)
		theplt.clf()
	else: 
		print 'error! there are no particles '		

#not funcitonal as far as i know.
def velocity_density_field(P, Rbox, x, y, z, vx, vy, vz, foutname, focusshell=2):
	print 'lets graph! ',foutname
	xmin = x - Rbox
	ymin = y - Rbox
	zmin = z - Rbox
	if (xmin < 0): xmin =0 
	if (ymin < 0): ymin = 0
	if (zmin < 0): zmin = 0
	xmax = x + Rbox
	ymax = y + Rbox
	zmax = z + Rbox
	
	tempx = np.linspace(xmin, xmax, 50)
	tempy = np.linspace(ymin, ymax, 50)
	
	velx = P['v'][:,0] - vx 
	vely = P['v'][:,1] - vy
	
	X,Y = meshgrid(tempx, tempy)
	U =theplt.hist2d(P['p'][:,0], P['p'][:,1], range=[[xmin,xmax],[ymin,ymax]], bins=50, norm=LogNorm(), vmin=1, weights=velx)
	V = theplt.hist2d(P['p'][:,0], P['p'][:,1], range=[[xmin,xmax],[ymin,ymax]], bins=50, norm=LogNorm(), vmin=1, weights=vely)
	M = theplt.hist2d(P['p'][:,0], P['p'][:,1], range=[[xmin,xmax],[ymin,ymax]], bins=50, norm=LogNorm(), vmin=1)
	Q = theplt.quiver( X[::3, ::3], Y[::3, ::3], U[::3, ::3], V[::3, ::3], M[::3, ::3], units='x', pivot='tip', linewidths=(2,), edgecolors=('k'), headaxislength=5, cmap=cm.get_cmap('hot'))
	theplt.colorbar()
	theplt.savefig(foutname)
	theplt.clf()
	
def dumb_stars(StarsP, NewStarsP, Rbox, x, y, z, foutname, a, xtoz = False, ytoz = False, alphalim = 1e4,  extrapoints = [], extrars=[] ):
	little_h = 0.7

	Sxx = StarsP['p'][:,0]/little_h
	Syy = StarsP['p'][:,1]/little_h
	Szz = StarsP['p'][:,2]/little_h

	Sxx -= x/little_h
	Syy -= y/little_h
	Szz -= z/little_h
	
	xmin = 0 - Rbox/little_h
	ymin = 0 - Rbox/little_h
	zmin = 0 - Rbox/little_h
                                                                    
	xmax = 0 + Rbox/little_h
	ymax = 0 + Rbox/little_h
	zmax = 0 + Rbox/little_h

	
	cut = (Sxx<xmax) * (Sxx>xmin) * (Syy<ymax) * (Syy > ymin) * (Szz < zmax) * (Szz > zmin)
	Sxx = Sxx[cut]
	Syy = Syy[cut]
	Szz = Szz[cut]	
	
	NSxx = StarsP['p'][:,0][NewStarsP * cut]/little_h
	NSyy = StarsP['p'][:,1][NewStarsP * cut]/little_h
	NSzz = StarsP['p'][:,2][NewStarsP * cut]/little_h

	NSxx -= x/little_h
	NSyy -= y/little_h
	NSzz -= z/little_h
	
	firstmin = xmin
	secondmin = ymin
	firstmax = xmax
	secondmax = ymax
	
	firstone = Sxx
	secondone = Syy
	
	firstoneN = NSxx
	secondoneN = NSyy

	
	
	if (xtoz):
		firstone = Szz
		firstoneN = NSzz
		firstmin = zmin
		firstmax = zmax
	elif (ytoz):
		secondone = Szz
		secondoneN = NSzz
		secondmin = zmin
		secondmax = zmax
	
	alpha_of_stars = min ([ alphalim/ float(len(Sxx)), 1.0])
	print 'alpha of stars', alpha_of_stars
	alpha_of_new_stars  =  min ([ alphalim / (20*float(len(NSxx))), 1.0])
	alpha_of_new_stars = max([alpha_of_new_stars, min([2.0*alpha_of_stars, 1.0])])
	print 'alpha of new stars ',alpha_of_new_stars
	theplt.figure(figsize=(11,8))
	theplt.plot(firstone*a, secondone*a, '.k', alpha = alpha_of_stars, rasterized=True)
	theplt.plot(firstoneN*a, secondoneN*a, '.r', alpha = alpha_of_new_stars, rasterized=True)
	theplt.ylim(secondmin*a , secondmax*a)
	theplt.xlim(firstmin*a, firstmax*a)

	if (len(extrapoints)>1):
		marker_ar = ['*', 'p', 's', '^', 'v']*5
		color_ar = ['b','y','g','w','r']*5
		count = 0
		for pt in extrapoints:
			theplt.plot([a*(pt[0] - x)/little_h], [a*(pt[1]- y)/little_h], marker=marker_ar[count], color=color_ar[count], ms=10, alpha=1)
			print 'plotted ',(pt[0] - x)/little_h, (pt[1]- y)/little_h
			print 'original ',pt
			ang = 0
			if (len(extrars) > 0):
				if (extrars[count] > 0):
					ra = a*extrars[count] / little_h
					rb = ra
					theX,theY=ellipse(ra,rb,ang,a*(pt[0] - x)/little_h,a*(pt[1]- y)/little_h)	
					theplt.plot(theX,theY,ls=':', color=color_ar[count] ,ms=1,linewidth=1)
			count+=1
	

	theplt.savefig(foutname)
	
def dumbest_stars(StarsP,  foutname, a, alphalim=1e4):
	Sxx = StarsP['p'][:,0]
	Syy = StarsP['p'][:,1]
	Szz = StarsP['p'][:,2]
	alpha_of_stars = min ([ alphalim/ float(len(Sxx)), 1.0])
	print 'alpha of stars', alpha_of_stars
	theplt.figure(figsize=(11,8))
	theplt.plot(Sxx*a, Syy*a, '.k', alpha = alpha_of_stars, rasterized=True)
	theplt.savefig(foutname)