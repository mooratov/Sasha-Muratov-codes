import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys
import scipy.stats as stats
import scipy.optimize as optimize

today = date.today()
print today
newtonG = 6.67384e-8

plt.rcParams['ps.fonttype'] = 42 
plt.rcParams['axes.linewidth'] = 2.0 #set the value globally


Time_average_ISM = False
use_no_loz = True
use_all_but_weak = True
show_no_zach = False


Inner_label = r'(r=$0.25R_{\rm vir}$)'
Outer_label = r'(r=$R_{\rm vir}$)'

#Inner_label = r'(r=25 kpc)'
#Outer_label = r'(r=50 kpc)'


skip_empty_symbols = True
noWV = True
little_h = 0.7
Nptsforanchor = 1
v_anchor_tol = 10
v_anchor = 60
magicmarker = 10
label_fontsize = 24
axis_fontsize = 22
tick_width = 5
tick_height = 16

show_single_fit_vc = False
show_double_fit_vc = True
show_single_fit_Mh = False
usemagic = False
usemed = True
use_Lprogs_infit = True
velocity_mode_95 = True
plot_inflow = False 

if (len(sys.argv) < 3):
	print 'synatx: blah.py method-index z'
	print '1: flux slope, 6: flux integral, 7: cross integral, 8: cross slope'
	sys.exit()

naked = False
fillnaked = False
dotoplabel = False

mi = int(sys.argv[1])
z = float(sys.argv[2])

def sasha_chisquared(expec, observ):
	tempnum = (expec - observ)**2
	thevar = stats.tvar(observ)
	MeanSquaredError = np.sum(tempnum) / (len(expec)-2)
	RootMeanSquaredError = np.sqrt(MeanSquaredError)
	tempnum /= thevar
	tempdenom = (len(expec)-2)
	tempreturn = np.sum(tempnum)/tempdenom
	return [tempreturn,RootMeanSquaredError]

def sasha_slope_error(expec, observ, x_observ):
	tempnum = (expec - observ)**2
	tempnum = np.sum(tempnum)
	tempnum /=  (len(expec) -2 )
	tempnum = np.sqrt(tempnum)
	thevar = stats.tvar(x_observ)
	thestdev = np.sqrt(thevar)
	tempnum /= thestdev
	return tempnum

	
finname_ar = []
listname = 'list_of_files-ML_halo.txt' 
f = open(listname)
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finname_ar.append(str(xsd[0]))

finnames = finname_ar
fintitle = finnames 

finname2_ar = []
listname = 'list_of_files-N_halo.txt' 
f = open(listname)
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finname2_ar.append(str(xsd[0]))
f.close()

if (not noWV):
	finname3_ar = []
	listname = 'list_of_files-WV_halo.txt' 
	if (velocity_mode_95):
		listname = 'list_of_files-WV_halo-95.txt' 
	f = open(listname)
	dars = f.readlines()
	for line in dars:
		xsd = line.split()
		finname3_ar.append(str(xsd[0]))
	f.close()

finname4_ar = []
listname = 'list_of_files-Zin_halo.txt' 
f = open(listname)
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finname4_ar.append(str(xsd[0]))
f.close()


#fintitle = [r'1e12  $M_{\odot}$', r'1e11  $M_{\odot}$', r'1e12 $M_{\odot}$', r'1e12 $M_{\odot}$', r'1e13 $M_{\odot}$', r'1e11 $M_{\odot}$', r'1e10 $M_{\odot}$', r'5e12 $M_{\odot}$']

ML_ar = []
M_ar = []
v_ar = []
Mstar_ar = [] 
v_ext_ar = []
ML_ext_ar = []
Mstar_ext_ar = []
M_ext_ar = []
Ulim_ar = []
veloc_ext_ar = []
veloc_ar = []
outflowmet_ar = []
outflowmet_ext_ar = []
ISMmet_ar = []
ISMmet_ext_ar = []

outflowmet9_ar = []
outflowmet9_ext_ar = []

ML9_ar = []
ML9_ext_ar = []

z_rel_ar = []
z_rel_ext_ar = []

MG_ar = []
MG_ext_ar = []

MG9_ar = []
MG9_ext_ar = []

ZG_ar = []
ZG_ext_ar = []

ZG9_ar = []
ZG9_ext_ar = []


#marker_ar = ['^', 'D', 'v', 's', 'h']
marker_ar = ['^', '^', '^', '^', '^']
marker_ar.extend(marker_ar)
marker_ar.extend(marker_ar)
marker_ar.extend(marker_ar)
marker_ar.extend(marker_ar)

color_ar = ['k', 'k', 'k', 'k', 'k']
color_ar.extend(color_ar)
color_ar.extend(color_ar)
color_ar.extend(color_ar)
color_ar.extend(color_ar)

count = 0
empcuts = []
fillcuts = []

ff = open('ML_from_slope_info.txt', 'w')
if (usemagic):
	ff.write('using magic \n')
elif (usemed):
	ff.write('using med \n')
else:
	ff.write('using endpoint \n')

ff.write('v_anchor, v_anchor_tol '+str(v_anchor)+'  '+str(v_anchor_tol)+'\n')

for name in finnames:
	print 'doing ',name
	f = open(name)
	g = open(finname2_ar[count])
	if (not noWV):
		i = open(finname3_ar[count])
	j = open(finname4_ar[count])

	if (not naked):
		if ( ('m10' in name) or('m11' in name) or ('m09' in name)):
			marker_ar[count] = marker_ar[count]
			Ulim_ar.append(False)
		elif 'dm' in name:
			marker_ar[count] = 's'
			Ulim_ar.append(False)
		elif 'TK' in name:
			marker_ar[count] = '>'
			Ulim_ar.append(False)
		elif 'Lprogs' in name:
			marker_ar[count] = '_'
			Ulim_ar.append(True)
		elif (('m12c' in name) or ('m12d' in name)):
			marker_ar[count] = 'v'
			print name, ' is going to a *'
			if (show_no_zach):
				marker_ar[count] = 'None'


			Ulim_ar.append(False)
		#elif 'm12q' in name:
		#	marker_ar[count] = 'x'
		#	Ulim_ar.append(False)
		else: 
			marker_ar[count] = 'v'
			Ulim_ar.append(False)

	color_ar[count] = 'r' 
	#print name, 'changing to blue'
	color_ar[count] = 'b'

	dars = np.loadtxt(f, ndmin=2)
	MLs = dars[:,mi]
	Ms = dars[:,2] / little_h
	vcs = dars[:,3]
	Mstar = dars[:,4] / little_h
	#z_ends = Mstar*0 + z
	#z_ends = dars[:,15]
	
	zend = dars[:,15]
	minzend = min(zend)
	
	if (('hiz' in name) and minzend >2.01):
		minzend = 2.01
	if (('medz' in name) and minzend >0.51):
		minzend = 0.51
	if (('loz' in name) and minzend >0.01):
		minzend = 0.01
	
	
	main = (zend <= minzend)
	#print zend, (z+1e-8)
	z_ends = zend
	fillcut = main
	empcut = np.invert(fillcut)
	fillcuts.append(fillcut)
	empcuts.append(empcut)
	
	dars2 = np.loadtxt(g,ndmin=2)
	if (not noWV):
		dars3 = np.loadtxt(i, ndmin=2)
	dars4 = np.loadtxt(j, ndmin=2)

	if (not noWV):
		veloc = dars3[:,1]
	outflowmet = dars4[:,2]
	if (Time_average_ISM):
		ISMmet = dars4[:,4]
	else: 
		ISMmet = dars4[:,3]
	ML9 = dars4[:,6]
	outflowmet9 = dars4[:,7]
	MG = dars4[:,8]
	MG9 = dars4[:,9]
	ZG =dars4[:,10]
	ZG9 = dars4[:,11]
	print len(ZG9), len(vcs)
	if (len(ZG9) != len(vcs)):
		print 'something be wrong '
		print name 
		sys.exit()
	#print '2nd dars',dars2
	
	magicz = dars2[:,1]
	magicMstar = dars2[:,2]/little_h
	magicvc = dars2[:,3]
	magicMh = dars2[:,4]/little_h
	medz = dars2[:,5]
	medMstar = dars2[:,6]/little_h
	medvc = dars2[:,7] 
	medMh = dars2[:,8]	/little_h
	#print finname2_ar[count]
	#print 'magic z ',magicz
	#print 'magc Mstar ',magicMstar
	#print 'magic vc ',magicvc
	#print 'magic Mh ',magicMh	
	ML_ar.append(MLs)
	if (not noWV):
		veloc_ar.append(veloc)
		veloc_ext_ar.extend(veloc)
	
	outflowmet_ar.append(outflowmet)
	outflowmet_ext_ar.extend(outflowmet)
	ISMmet_ar.append(ISMmet)
	ISMmet_ext_ar.extend(ISMmet)	
	outflowmet9_ar.append(outflowmet9)
	outflowmet9_ext_ar.extend(outflowmet9)	
	ML9_ar.append(ML9)
	ML9_ext_ar.extend(ML9)
	
	MG_ar.append(MG)
	MG9_ar.append(MG9)
	ZG_ar.append(ZG)
	ZG9_ar.append(ZG9)
	
	MG_ext_ar.extend(MG)
	MG9_ext_ar.extend(MG9)
	ZG_ext_ar.extend(ZG)
	ZG9_ext_ar.extend(ZG9)
	
	
	
	if ((not Ulim_ar[count]) or use_Lprogs_infit):
		ML_ext_ar.extend(MLs)
	
	if (medz[0] > 2):
		color_ar[count] =  'k'
	elif (medz[0] > 1):
		color_ar[count] = 'b'
	else:
		color_ar[count] = 'r'
		print 'lo redshift check ', name, MLs
	
	#countB = 0
	#while (countB < len(medz)):
	#	if (medz < 0.5 and Ms > 5e11):
	#		marker_ar[count] = '-'
	
	if (usemagic):
		M_ar.append(magicMh)
		v_ar.append(magicvc)
		Mstar_ar.append(magicMstar)
		z_rel_ar.append(magicz)

		if ((not Ulim_ar[count])  or use_Lprogs_infit):
			M_ext_ar.extend(magicMh)
			v_ext_ar.extend(magicvc)
			Mstar_ext_ar.extend(magicMstar)
			z_rel_ext_ar.extend(magicz)
		
		#print 'lala ', M_ar, M_ext_ar
	elif (usemed):
		M_ar.append(medMh)
		v_ar.append(medvc)
		Mstar_ar.append(medMstar)
		z_rel_ar.append(medz)
		if ((not Ulim_ar[count])  or use_Lprogs_infit):
			M_ext_ar.extend(medMh)
			v_ext_ar.extend(medvc)
			Mstar_ext_ar.extend(medMstar)
			z_rel_ext_ar.extend(medz)
	else:		
		M_ar.append(Ms)
		v_ar.append(vcs)
		Mstar_ar.append(Mstar)
		z_rel_ar.append(z_ends)
		if ((not Ulim_ar[count]) or use_Lprogs_infit):
			M_ext_ar.extend(Ms)
			v_ext_ar.extend(vcs)
			Mstar_ext_ar.extend(Mstar)
			z_rel_ext_ar.extend(z_ends)
	
	

	
	f.close()
	g.close()
	if (not noWV):
		i.close()
	j.close()
	count +=1

count = 0


ML_ar = np.array(ML_ar)
M_ar = np.array(M_ar)
Mstar_ar = np.array(Mstar_ar)
v_ar = np.array(v_ar)

if (not noWV):
	veloc_ar = np.array(veloc_ar)
	veloc_ext_ar = np.array(veloc_ext_ar)

outflowmet_ar = np.array(outflowmet_ar)
outflowmet_ext_ar = np.array(outflowmet_ext_ar)

ISMmet_ar = np.array(ISMmet_ar)
ISMmet_ext_ar = np.array(ISMmet_ext_ar)

outflowmet9_ar = np.array(outflowmet9_ar)
outflowmet9_ext_ar = np.array(outflowmet9_ext_ar)


ML9_ar = np.array(ML9_ar)
ML9_ext_ar = np.array(ML9_ar)

metrat_ar = outflowmet_ar/ML_ar
metrat9_ar = outflowmet9_ar / ML9_ar

MG_ar = np.array(MG_ar) * -1 
MG_ext_ar = np.array(MG_ext_ar) * -1 

print 'mg ar ',MG_ar

MG9_ar = np.array(MG9_ar) * -1 
MG9_ext_ar = np.array(MG9_ext_ar) * -1 

ZG_ar = np.array(ZG_ar) * -1 
ZG_ext_ar = np.array(ZG_ext_ar) * -1 

ZG9_ar = np.array(ZG9_ar) * -1 
ZG9_ext_ar = np.array(ZG9_ext_ar) * -1 

metratG_ar = ZG_ar/MG_ar
metratG9_ar = ZG9_ar / MG9_ar

print 'metrat G',metratG_ar

v_ext_ar = np.array(v_ext_ar)
ML_ext_ar = np.array(ML_ext_ar)
M_ext_ar = np.array(M_ext_ar)
Mstar_ext_ar = np.array(Mstar_ext_ar)
z_rel_ar = np.array(z_rel_ar)
z_rel_ext_ar = np.array(z_rel_ext_ar)

hiz_cut = z_rel_ext_ar > 2
medz_cut = (z_rel_ext_ar < 2) * (z_rel_ext_ar > 0.5)
loz_cut = (z_rel_ext_ar < 0.5)

R_vir_derived_ar_ext = 1.0 /  (v_ext_ar * 1e5)**2.0
R_vir_derived_ar_ext *= little_h / 3.08e21
R_vir_derived_ar_ext *= (newtonG * M_ext_ar *  2e33 / little_h)

#print 'Rvirs derived ',R_vir_derived_ar
magic_rat = ( M_ext_ar*2e33 / (R_vir_derived_ar_ext*3.08e21)**3.0)
#print 'magic number? ',magic_rat

magic_rat_mean = np.mean(magic_rat)
magic_rat_mean_hiz = np.mean(magic_rat[hiz_cut])
magic_rat_mean_medz = np.mean(magic_rat[medz_cut])
magic_rat_mean_loz = np.mean(magic_rat[loz_cut])



#print 'omega lengths ', len(v_ext_ar), len(Mstar_ext_ar)
ff.write( 'omega lengths ' + str( len(v_ext_ar)) +'  '+str( len(Mstar_ext_ar)) + '\n')

v_ar *= ((1.0+z_rel_ar))**0.5
v_ext_ar *= ((1.0+z_rel_ext_ar))**0.5
#only fixing for screwed up thing vcirc 3-18-14

big_cut = v_ext_ar > v_anchor
little_cut = v_ext_ar <= v_anchor
if (use_no_loz and (not use_all_but_weak)):
	big_cut *= (hiz_cut + medz_cut)
	little_cut  *= (hiz_cut + medz_cut)


print 'little guys ',len(v_ext_ar[little_cut])
print 'big guys ',len(v_ext_ar[big_cut])

ML_ar_log = np.log10(ML_ext_ar)
v_ar_log = np.log10(v_ext_ar)
M_ar_log = np.log10(M_ext_ar)
Mstar_ar_log = np.log10(Mstar_ext_ar)
if (not noWV):
	veloc_ar_log = np.log10(veloc_ext_ar)
outflowmet_ar_log = np.log10(outflowmet_ext_ar)

vorder = np.argsort(v_ext_ar)
big_order_cut = v_ext_ar[vorder] > v_anchor
little_order_cut = v_ext_ar[vorder] <= v_anchor


count = 0
ML_avg = 0
while (count < Nptsforanchor):
	ML_avg += ML_ar_log[vorder][little_order_cut][-1*count - 1]
	ML_avg += ML_ar_log[vorder][big_order_cut][count]

	#print 'using ML ',ML_ext_ar[vorder][big_order_cut][count], ML_ext_ar[vorder][little_order_cut][-1*count - 1]
	#print 'using v',v_ext_ar[vorder][big_order_cut][count], v_ext_ar[vorder][little_order_cut][-1*count - 1]
	count +=1

#dividing by 2*Npoints for anchor. Would be 6 for Nptsperanchor=3 because using 3 values from each side 
ML_avg /= (2.0*Nptsforanchor)

near_vanchor = np.fabs(v_ext_ar - v_anchor) < v_anchor_tol
#print 'gonna use ',len(ML_ext_ar[near_vanchor]),' poins for anchor'
#ff.write('gonna use '+str(len(ML_ext_ar[near_vanchor]))+' poins for anchor \n')
ML_avg_hiz = np.mean(ML_ar_log[near_vanchor*hiz_cut])
ML_avg_medz = np.mean(ML_ar_log[near_vanchor*medz_cut])
ML_avg_loz = np.mean(ML_ar_log[near_vanchor*loz_cut])
ML_avg = np.mean(ML_ar_log[near_vanchor])

ML_avg = 10**ML_avg
ML_avg_hiz = 10** ML_avg_hiz
ML_avg_medz = 10**ML_avg_medz
ML_avg_loz = 10** ML_avg_loz


print 'ML avg at hiz, medz, loz, allz ',ML_avg_hiz,ML_avg_medz,ML_avg_loz, ML_avg


#ML_avg = ML_avg_hiz
fixpoint = ML_avg_hiz
#'ML avg'
#print ML_avg


def residuals(b , y, x):
    err = y - (b[1]*x + b[0])
    return err
    
def residualsZ(b , y, x, zz):
    err = y - (b[1]*x + b[0] + b[2]*(1+zz))
    return err

def residualsZZ(b , y, x, zz):
    err = y - (b[1]*x + b[0] + b[2]*np.log10(1+zz))
    return err

    
p0 = [-5, 1]
p0z = [-5, 1, -1]
p1 = [-5, -1]
p2 = [-5, 1, -1]
p3 = [-1, 1, 1]

stoopcut = ML_ar_log > -999
if (use_no_loz):
	stoopcut  *= (hiz_cut + medz_cut)
if (use_all_but_weak):
	stoopcut = ML_ar_log > 0

#with 1+z dependence
Mstarfit = optimize.leastsq(residualsZ,p3[:], args=(ML_ar_log[stoopcut], Mstar_ar_log[stoopcut], z_rel_ext_ar[stoopcut] ))


Mstar_int = Mstarfit[0][0]
Mstar_slope = Mstarfit[0][1]
Mstar_z_slope = Mstarfit[0][2]


#without 1+z dependence
Mstarfit = optimize.leastsq(residuals,p0[:], args=(ML_ar_log[stoopcut], Mstar_ar_log[stoopcut] ))

Mstar_int = Mstarfit[0][0]
Mstar_slope = Mstarfit[0][1]



#Mfit = optimize.leastsq(residuals,p0[:], args=(ML_ar_log[stoopcut], M_ar_log[stoopcut]))


#with 1+z dependence
Mstarfa = optimize.leastsq(residualsZ,p3[:], args=(ML_ar_log[stoopcut], Mstar_ar_log[stoopcut], z_rel_ext_ar[stoopcut]), full_output=True)
#without 1+z dependence
Mstarfa = optimize.leastsq(residuals,p0[:], args=(ML_ar_log[stoopcut], Mstar_ar_log[stoopcut]), full_output=True)



#print 'Mstarfa ',Mstarfa
Mstarfacov = Mstarfa[1]
print 'Mstarfa cov ',Mstarfacov
perrS = np.sqrt(np.diag(Mstarfacov))
print 'full perrS', perrS
print 'perrS ',perrS[0]


print 'Mstar int, Mstar slope', Mstar_int, Mstar_slope, Mstar_z_slope

ff.write('Mstar int, Mstar slope ' + str(Mstar_int) + ' ' + str(Mstar_slope) + '\n')
#M_int = Mfit[0][0]
#M_slope = Mfit[0][1]

#print 'Mh int, Mh slope ',M_int, M_slope

#ff.write('Mh int, Mh slope '+str( M_int) + ' '+str(M_slope)+'\n')


logMspace = np.linspace(8, 14 ,num=50)
logMstarspace = np.linspace(4, 12, num=50)
Mspace = 10**logMspace
Mstarspace = 10**logMstarspace

#M1fit = pow(Mspace, M_slope) * 10**M_int
Mstar1fit = pow(Mstarspace, Mstar_slope) * 10**Mstar_int


#predicted_Out_M1 = pow(M_ext_ar, M_slope)  * 10**M_int
#predicted_Out_M1  = np.log10(predicted_Out_M1)
#[chisquared_M1, RMSE_aM] = sasha_chisquared(predicted_Out_M1, ML_ar_log)

#print 'chi squared, RMSE Mh ',chisquared_M1, RMSE_aM
#ff.write('chi squared, RMSE Mh  '+str( chisquared_M1) + ' '+ str(RMSE_aM)+'\n')


predicted_Out_Mstar1 = pow(Mstar_ext_ar, Mstar_slope)  * 10**Mstar_int
predicted_Out_Mstar1  = np.log10(predicted_Out_Mstar1)
[chisquared_Mstar1, RMSE_aMstar] = sasha_chisquared(predicted_Out_Mstar1, ML_ar_log)
print 'chi squared, RMSE, error1 Mstar ',chisquared_Mstar1, RMSE_aMstar, perrS[0]
ff.write('chi squared, RMSE, error1 Mstar  '+str( chisquared_Mstar1) + ' '+ str(RMSE_aMstar)+'  '+str(perrS[0])+'\n')



allfit = optimize.leastsq(residuals,p0[:], args=(ML_ar_log[stoopcut], v_ar_log[stoopcut]))
allfit 


print allfit
all_yint = allfit[0][0]
all_slope = allfit[0][1]

print 'vc two-param int, slope ',all_yint, all_slope
ff.write('vc two-param int, slope   '+str( all_yint) + ' '+ str(all_slope)+'\n')

vc2fa = optimize.leastsq(residuals,p0[:], args=(ML_ar_log[stoopcut], v_ar_log[stoopcut]), full_output=True)
#print 'vc2fa ',vc2fa
vc2facov = vc2fa[1]
print 'vc2fa cov ',vc2facov
perrVc2 = np.sqrt(np.diag(vc2facov))
print 'full perrVc2', perrVc2
print 'perrVc2 ',perrVc2[0]

if (not noWV):
	veloc_fit = optimize.leastsq(residuals,p0[:], args=(veloc_ar_log[stoopcut], v_ar_log[stoopcut]), full_output=True)

	print veloc_fit
	veloc_yint = veloc_fit[0][0]
	veloc_slope = veloc_fit[0][1]
	veloc_cov = veloc_fit[1]
	perr_veloc = np.sqrt(np.diag(veloc_cov))
	print 'perr wind velocity ',perr_veloc

	print 'wind velocity int, slope ',veloc_yint, veloc_slope
	ff.write('wind velocity int, slope    '+str( veloc_yint) + ' '+ str(veloc_slope)+'\n')


#fixpoint = v_anchor**all_slope * 10** all_yint
print 'fixpoint ', fixpoint
ff.write('fixpoint   '+str( fixpoint) +'\n')


def residuals2(b, y, x):
	err = y - (b[0]*(x) + np.log10(fixpoint))
	return err

def residuals3(b, y, x):
	err = y - np.log10(x**b[0] + fixpoint)
	return err
	
def residuals4(b, y, x):
	err = y -   (np.log10((fixpoint / v_anchor**b[0])) + b[0]*x) 
	return err

def residuals5(b, y, x):
	#print 'lalala 'fixpoint
	err = y - ((np.log10(fixpoint) - b[0]*np.log10(v_anchor)) + b[0]*x)
	return err
	
def residuals6(b, y, x,zz):
	#print 'lalala 'fixpoint
	err = y - ((np.log10(get_fixpoint(zz)) - b[0]*np.log10(v_anchor) - b[1]*(1.0+3.0)) + b[0]*x + b[1]*(1+zz))
	return err
def get_fixpoint(zz):
	thereturn = zz*0
	countzz = 0
	while (countzz < len(zz)):	
		if (zz[countzz] > 2):
			thereturn[countzz] = ML_avg_hiz
		elif (zz[countzz] < 2 and zz[countzz] > 1):
			thereturn[countzz] = ML_avg_medz
		else:
			thereturn[countzz] =  ML_avg_loz
		countzz +=1
	return thereturn

def residuals7(b, y, x,zz):
	err = y - ((np.log10(fixpoint) - b[0]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[0]*x + b[1]*np.log10(1+zz))
	return err
	
def residuals8(b, y, x,zz):
	err = y 
	little_fit_cut =  (x <= v_anchor)
	big_fit_cut = (x > v_anchor)
	err[little_fit_cut] = y[little_fit_cut] - ((np.log10(fixpoint) - b[2]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[0]*x[little_fit_cut] + b[1]*np.log10(1+zz[little_fit_cut]))
	err[big_fit_cut] = y[big_fit_cut] - ((np.log10(fixpoint) - b[0]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[2]*x[big_fit_cut] + b[1]*np.log10(1+zz[big_fit_cut]))
	
def residuals10(b, y, x,zz):
	err = y 
	little_fit_cut =  (x <= v_anchor)
	big_fit_cut = (x > v_anchor)
	err = y - ((np.log10(fixpoint) - b[2]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[0]*x + b[1]*np.log10(1+zz))
	err[big_fit_cut] = y[big_fit_cut] - ((np.log10(fixpoint) - b[0]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[2]*x[big_fit_cut] + b[1]*np.log10(1+zz[big_fit_cut]))


	#err[big_fit_cut] = y[big_fit_cut] - ((np.log10(fixpoint) - b[2]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[2]*x[big_fit_cut] + b[1]*np.log10(1+zz[big_fit_cut]))
	#print len(err), b[0], b[1], b[2]
	return err
	
def residuals9(b, y, x, yy, xx, zz, zzA):
	#err = np.arange(len(y)+len(yy))
	little_fit_cut =  (x <= v_anchor)
	big_fit_cut = (x > v_anchor)
	err  = y - ((np.log10(fixpoint) - b[0]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[0]*x + b[1]*np.log10(1.0+zz)) 
	err2 = yy - ((np.log10(fixpoint) - b[2]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[2]*xx + b[1]*np.log10(1.0+zzA))
	err = np.vstack((err, err2))
	#err[little_fit_cut] = y - ((np.log10(fixpoint) - b[0]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[0]*x + b[1]*np.log10(1.0+zz))
	#err[big_fit_cut] = yy - ((np.log10(fixpoint) - b[2]*np.log10(v_anchor) - b[1]*np.log10(1.0+3.0)) + b[2]*xx + b[1]*np.log10(1.0+zzA))
	return err
#def residuals8(b, y, x,zz):

#bigfit = optimize.leastsq(residuals3,p0[0], args=(ML_ar_log[big_cut], v_ext_ar[big_cut])) - archaic


bigfit = optimize.leastsq(residuals7,p1[:], args=(ML_ar_log[big_cut], v_ar_log[big_cut], z_rel_ext_ar[big_cut]))
bigfa =  optimize.leastsq(residuals7,p1[:], args=(ML_ar_log[big_cut], v_ar_log[big_cut], z_rel_ext_ar[big_cut]), full_output=True)

#bigfit = optimize.leastsq(residuals5,p0[0], args=(ML_ar_log[big_cut], v_ar_log[big_cut]))
#bigfa = optimize.leastsq(residuals5,p0[0], args=(ML_ar_log[big_cut], v_ar_log[big_cut]), full_output=True)

#print 'bigfa ',bigfa
bigfacov = bigfa[1]
if ( len(bigfacov)>0):
 perr = np.sqrt(np.diag(bigfacov))
else: perr = [0]
print 'full perr ',perr
print 'perr ',perr[0]

print 'vc big fit slope, yint(fixed) ',bigfit, fixpoint
ff.write('vc big fit slope, yint(fixed) '+str(bigfit[0][0])+'  '+ str(bigfit[0][1])+'  '+ str(fixpoint)+'\n')
#big_yint = bigfit[0][0]
#big_slope = bigfit[0][1]
big_slope = bigfit[0][0]
big_yint = fixpoint
big_yint = np.log10(fixpoint)
big_yint = fixpoint / v_anchor**big_slope
#big_yint = 10**(np.log10(fixpoint) - big_slope*np.log10(75))

#littlefit = optimize.leastsq(residuals3,p0[0], args=(ML_ar_log[little_cut], v_ext_ar[little_cut])) - archaic

#littlefit = optimize.leastsq(residuals5,p0[0], args=(ML_ar_log[little_cut], v_ar_log[little_cut]))
#littlefa = optimize.leastsq(residuals5,p0[0], args=(ML_ar_log[little_cut], v_ar_log[little_cut]), full_output=True)

littlefit = optimize.leastsq(residuals7,p1[:], args=(ML_ar_log[little_cut], v_ar_log[little_cut],  z_rel_ext_ar[little_cut]))
littlefa = optimize.leastsq(residuals7,p1[:], args=(ML_ar_log[little_cut], v_ar_log[little_cut],  z_rel_ext_ar[little_cut]), full_output=True)

print littlefa
#bothfits =  optimize.leastsq(residuals9,p2[:], args=(ML_ar_log[little_cut], v_ar_log[little_cut], ML_ar_log[big_cut], v_ar_log[big_cut], z_rel_ext_ar[little_cut],  z_rel_ext_ar[big_cut]), full_output=True)
bothfits = optimize.leastsq(residuals10,p2[:], args=(ML_ar_log, v_ar_log ,z_rel_ext_ar), full_output=True)

#print 'littlefa ',littlefa
lilfacov = littlefa[1]
if ( len(lilfacov) > 0):
 perrl = np.sqrt(np.diag(lilfacov))
else: perrl = [0]
print 'full perrl ',perrl
print 'perr ',perrl[0]


print 'vc little fit, yint(fixed) ',littlefit, fixpoint
ff.write('vc little fit slope, yint(fixed) '+str(littlefit[0][0])+'  '+str(littlefit[0][1])+'  '+ str(fixpoint)+'\n')

print bothfits
bothfacov = bothfits[1]
if ( len(bothfacov) > 0):
 perrboth = np.sqrt(np.diag(bothfacov))
else: perrboth = [0]
print 'full perr both ',perrboth
print 'perr both',perrboth[0]


print 'vc both fit ',perrboth
ff.write('vc little fit slope, yint(fixed) '+str(littlefit[0][0])+'  '+str(littlefit[0][1])+'  '+ str(fixpoint)+'\n')

#little_yint = littlefit[0][0]
#little_slope = littlefit[0][1]

little_slope = littlefit[0][0]
little_yint = fixpoint
little_yint = np.log10(fixpoint)
little_yint = fixpoint / v_anchor**little_slope

yint_diffrat_medz = (1.0 + 1.25)**1.25 / (1.0 + 3.0) ** 1.25
yint_diffrat_loz = ( 1.0 + 0.25)**1.25 / (1.0 + 3.0) ** 1.25

Manchor_hiz = 1.77e10
Manchor_medz = 4.03e10
Manchor_loz = 8.27e10


Mlittle_yint_hiz = fixpoint / Manchor_hiz**(little_slope/3.0)
Mbig_yint_hiz = fixpoint / Manchor_hiz**(big_slope/3.0)

Mlittle_yint_medz = fixpoint / Manchor_medz**(little_slope/3.0)  * yint_diffrat_medz
Mbig_yint_medz = fixpoint / Manchor_medz**(big_slope/3.0) * yint_diffrat_medz

Mlittle_yint_loz = fixpoint / Manchor_loz**(little_slope/3.0)  * yint_diffrat_loz
Mbig_yint_loz = fixpoint / Manchor_loz**(big_slope/3.0) * yint_diffrat_loz


#little_yint = 10**(np.log10(fixpoint) - little_slope*np.log10(75))
print 'little y_int ', little_yint

#fit by eye for -1 and -2 power law
xspace = np.linspace(10, 1000 ,num=50)
v1fit = pow(xspace, -1) * 10**2.75
v2fit = pow(xspace, -2) * 10**4.75 

#fit for low-mass end
#l1fit = pow(xspace, little_slope) * 10**little_yint
l1space = np.linspace(10, v_anchor, num=50)


amean = 0

zmean = np.mean(z_rel_ext_ar)
if (usemed):
	zmean = stats.mstats.mode(z_rel_ext_ar[hiz_cut])[0][0]
print 'using zmean ',zmean
amean = 1.0 / (1.0 + zmean)
Mv_l1space = (l1space *1e5 )**3.0 * amean**(3.0/2.0) / (newtonG**(3.0/2.0) * magic_rat_mean_hiz**0.5)
Mv_l1space /= 2e33


Mv_xspace = (xspace *1e5 )**3.0 * amean**(3.0/2.0) / (newtonG**(3.0/2.0) * magic_rat_mean_hiz**0.5)
Mv_xspace /= 2e33

truemagic = Mv_xspace / xspace**3.0
#print 'true magic ',truemagic
ff.write('truemagic   '+str( truemagic[0]) +'\n')


Manchor = (v_anchor *1e5 )**3.0 * amean**(3.0/2.0) / (newtonG**(3.0/2.0) * magic_rat_mean_hiz**0.5) 
Manchor /= 2e33


Ml1space_hiz_log = np.linspace(np.log10(4e8), np.log10(Manchor_hiz), num=50)
Mb1space_hiz_log = np.linspace(np.log10(Manchor_hiz), np.log10(3e12), num=50)

Ml1space_medz_log = np.linspace(np.log10(4e8), np.log10(Manchor_medz), num=50)
Mb1space_medz_log = np.linspace(np.log10(Manchor_medz), np.log10(3e12), num=50)

Ml1space_loz_log = np.linspace(np.log10(4e8), np.log10(Manchor_loz), num=50)
Mb1space_loz_log = np.linspace(np.log10(Manchor_loz), np.log10(3e12), num=50)

Ml1space_hiz = 10**Ml1space_hiz_log
Mb1space_hiz = 10**Mb1space_hiz_log 

Ml1space_medz = 10**Ml1space_medz_log
Mb1space_medz = 10**Mb1space_medz_log

Ml1space_loz = 10**Ml1space_loz_log
Mb1space_loz = 10**Mb1space_loz_log 

print 'Manchor ', Manchor
ff.write('Manchor   '+str( Manchor) +'\n')

print 'v_c fit: little y slope, int no need for 10** for int', little_slope, little_yint

#fit for high-mass end
#b1fit = pow(xspace, big_slope) * 10**big_yint
b1space = np.linspace(v_anchor, 1000, num=50)
Mv_b1space = (b1space *1e5 )**3.0 * amean**(3.0/2.0) / (newtonG**(3.0/2.0) * magic_rat_mean_hiz**0.5)
Mv_b1space /= 2e33

b1fit = pow(b1space, big_slope) *  big_yint
b1fit_medz =  pow(b1space, big_slope) *  big_yint * (1.0 + 1.25)**1.25 / (1.0+3.0)**1.25
b1fit_loz =  pow(b1space, big_slope) *  big_yint * (1.0 + 0.25)**1.25 / (1.0+3.0)**1.25

l1fit = pow(l1space, little_slope) *  little_yint

l1fit_medz =  pow(l1space, little_slope) *  little_yint * (1.0 + 1.25)**1.25 / (1.0+3.0)**1.25
l1fit_loz =  pow(l1space, little_slope) *  little_yint * (1.0 + 0.25)**1.25 / (1.0+3.0)**1.25

Ml1fit  = Ml1space_hiz**(little_slope/3.0) * Mlittle_yint_hiz  
Mb1fit = Mb1space_hiz**(big_slope/3.0) * Mbig_yint_hiz

Ml1fit_medz = Ml1space_medz**(little_slope/3.0) *  Mlittle_yint_medz  
Mb1fit_medz = Mb1space_medz**(big_slope/3.0) * Mbig_yint_medz

Ml1fit_loz = Ml1space_loz**(little_slope/3.0) *  Mlittle_yint_loz  
Mb1fit_loz = Mb1space_loz**(big_slope/3.0) * Mbig_yint_loz



print 'v_c fit: big y slope, int no need for 10** for int', big_slope, big_yint

#l1fit = pow(xspace, little_slope )+ little_yint
#b1fit = xspace ** big_slope + big_yint

a1fit = pow(xspace, all_slope) * 10**all_yint


predicted_Out_a1 = pow(v_ext_ar, all_slope)  * 10**all_yint
predicted_Out_a1  = np.log10(predicted_Out_a1)
[chisquared_a1, RMSE_a1] = sasha_chisquared(predicted_Out_a1, ML_ar_log)


predicted_Out_lb1 = pow(v_ext_ar, little_slope)  * little_yint
predicted_Out_lb1[big_cut] =  pow(v_ext_ar[big_cut], big_slope)  * big_yint

predicted_Out_lb1  = np.log10(predicted_Out_lb1)


[chisquared_lb1, RMSE_lb1] = sasha_chisquared(predicted_Out_lb1, ML_ar_log)



predicted_Out_l1 = pow(v_ext_ar[little_cut], little_slope)  * little_yint
predicted_Out_l1  = np.log10(predicted_Out_l1)

[chisquared_l1, RMSE_l1] = sasha_chisquared(predicted_Out_l1, ML_ar_log[little_cut])
ErrOnSlope_l1 = sasha_slope_error(predicted_Out_l1,ML_ar_log[little_cut],v_ar_log[little_cut])

predicted_Out_b1 = pow(v_ext_ar[big_cut], big_slope)  * big_yint
predicted_Out_b1  = np.log10(predicted_Out_b1)

[chisquared_b1, RMSE_b1] = sasha_chisquared(predicted_Out_b1, ML_ar_log[big_cut])
ErrOnSlope_b1 = sasha_slope_error(predicted_Out_b1,ML_ar_log[big_cut],v_ar_log[big_cut])

print 'vc single two-parmater fit chi squared, RMSE ',chisquared_a1, RMSE_a1
ff.write('vc single two-parmater fit chi squared, RMSE '+str( chisquared_a1) + ' '+str(RMSE_a1)+'\n')

print 'vc combined low and high one-parmater anchored fit chi squared, RMSE ',chisquared_lb1, RMSE_lb1
ff.write('vc combined low and high one-parmater anchored fit chi squared, RMSE  '+str( chisquared_lb1) + ' '+ str(RMSE_lb1)+'\n')


print 'vc low velocity anchored fit chi squared, RMSE, error on slope , error from residuals',chisquared_l1, RMSE_l1, ErrOnSlope_l1, perrl[0]
ff.write('vc low velocity anchored fit chi squared, RMSE, error on slope, error from residuals  '+str( chisquared_l1) + ' '+str(RMSE_l1)+' '+str(ErrOnSlope_l1)+' '+str(perrl[0])+'\n')

print 'vc high velocity anchored fit chi squared, RMSE, error on slope, error from residuals ',chisquared_b1, RMSE_b1,ErrOnSlope_b1, perr[0]
ff.write('vc high velocity anchored fit chi squared, RMSE, error on slope, error from residuals   '+str( chisquared_b1) + ' '+ str(RMSE_b1)+' '+str(ErrOnSlope_b1)+' '+str(perr[0])+'\n')


thar = fintitle 
fig =  plt.figure(figsize=(10,9))


maxplot = len(M_ar)
count = 0
while (count < maxplot):

	if (fillnaked):
		plt.plot(M_ar[count], ML_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(M_ar[count][fillcuts[count]], ML_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(M_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(M_ar[count], ML_ar[count], yerr=(ML_ar[count]*0.4, ML_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

		
	#ERRORBAR EXPERIMENTS
	
	#plt.errorbar(M_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], xerr=[0, 2], yerr=[2,0])
	#plt.errorbar(M_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], xerr=2.0, yerr=2.0)
	#plt.errorbar(M_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], xerr=[(M_ar[count][empcuts[count]] -1), ( M_ar[count][empcuts[count]] +2)], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none')
	thexerr = 0.5*M_ar[count][empcuts[count]]*1
	thexerr2 = 0 * M_ar[count][empcuts[count]]*1
	thexerr3 = -1 *  M_ar[count][empcuts[count]]*1
	#plt.errorbar(M_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], xerr=[thexerr, thexerr], yerr =[thexerr, thexerr], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none')
	#plt.errorbar(M_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], xerr=[thexerr,thexerr2], yerr =[thexerr, thexerr], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none')
    


	count+=1
#print Mspace, len(Mspace)
#print M1fit, len(M1fit)
if (show_single_fit_Mh):
	#plt.plot(Mspace, M1fit, color='k', ls= '--', lw =2)
	plt.plot(Mv_xspace, a1fit, color = 'k', ls ='--', lw=2)
if (show_double_fit_vc):
	#plt.plot(Mv_l1space, l1fit, color='b', ls='--', lw=2)
	#plt.plot(Mv_b1space, b1fit, color='r', ls ='--', lw=2)
	plt.plot(Ml1space_hiz, Ml1fit, color='k', ls=':', lw=2)
	plt.plot(Mb1space_hiz, Mb1fit, color='k', ls=':', lw=2)

	plt.plot(Ml1space_medz, Ml1fit_medz, color='b', ls=':', lw=2)
	plt.plot(Mb1space_medz, Mb1fit_medz, color='b', ls=':', lw=2)
	
	plt.plot(Ml1space_loz, Ml1fit_loz, color='r', ls=':', lw=2)
	plt.plot(Mb1space_loz, Mb1fit_loz, color='r', ls=':', lw=2)

plt.ylim(1e-1 , 1e3)
plt.xlim(5e8, 2e12)
plt.xscale('log')
plt.yscale('log')




xlabel_str = ''
if (usemagic):
	xlabel_str = r'Halo Mass $(M_{\odot})$ at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'Halo Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'Halo Mass $(M_{\odot})$ at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)
	
plt.ylabel(r'$\eta$ = $\dot{M}_{out}$/SFR',fontsize=label_fontsize)

#plt.text(10**11.0, 10**2.2, r'$\eta$ = $10^{%s}$'%"{0:.3g}".format(M_int)+r'$\times  M_{h}^{%s}$'%"{0:.3g}".format(M_slope)+"\n"+ r'$\chi^2$ = '+"{0:.3e}".format(chisquared_M1) )

#plt.legend(thar, loc = 'best', prop={'size':10})


ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/12.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/12.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/12.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


plt.savefig('AllML'+str(mi)+'_log_int_'+ str(today) +'z'+str(int(z))+'.pdf')
#plt.show()
plt.clf()
fig =  plt.figure(figsize=(10,9))

count = 0
while (count < maxplot):

	if (fillnaked):
		plt.plot(v_ar[count], ML_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(v_ar[count], MG_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)

	else:
		plt.plot(v_ar[count][fillcuts[count]], ML_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(v_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(v_ar[count][fillcuts[count]], MG_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(v_ar[count][empcuts[count]], MG_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(v_ar[count], ML_ar[count], yerr=(ML_ar[count]*0.4, ML_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1

#plt.plot(v_ar[0], ML_ar[0], 'Dk')
#plt.plot(v_ar[1], ML_ar[1], 'Db')
#plt.plot(v_ar[2], ML_ar[2], 'sk')
#plt.plot(v_ar[3], ML_ar[3], 'sb')
#plt.plot(v_ar[4], ML_ar[4], '^k')
#plt.plot(xspace, v1fit, color='k', ls='--', markevery=1, lw=2)
#plt.plot(xspace, v2fit, color='k', ls='--', markevery=1, lw=2)
if (show_single_fit_vc):
	plt.plot(xspace, a1fit, color='k', ls= '--', lw =2)
if (show_double_fit_vc):
	plt.plot(l1space, l1fit, color='k', ls=':', lw=2)
	plt.plot(b1space, b1fit, color='k', ls =':', lw=2)
	plt.plot (l1space, l1fit_medz, color ='b', ls=':', lw=2)
	plt.plot (b1space, b1fit_medz, color ='b', ls=':', lw=2)
	plt.plot (l1space, l1fit_loz, color ='r', ls=':', lw=2)
	plt.plot (b1space, b1fit_loz, color ='r', ls=':', lw=2)	

#plt.plot(v_ar[7], ML_ar[7], 'hb')

plt.xlim(1e1, 3e2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.xscale('log')

ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')
x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/3.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/3.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/3.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')
xlabel_str = ''
if (usemagic):
	xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
else:
	xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.xlabel(r'$V_c$ (km/s) at z='+str(int(z)), fontsize=20)
plt.ylabel(r'$\eta$ = $\dot{M}_{out}$/SFR',fontsize=label_fontsize)
if (dotoplabel):
	plt.title('z='+str(int(z)), fontsize=label_fontsize)
#plt.legend(thar, loc = 'best', prop={'size':10})

#text fun
#plt.text(10**1.5, 10**2.2, "v < "+str(v_anchor)+" km/s slope = "+"{0:.3e}".format(little_slope)+r'$\pm$ '+"{0:.3e}".format(RMSE_l1)+"\n v > "+str(v_anchor)+" km/s slope = "+"{0:.3e}".format(big_slope)+r'$\pm$ '+"{0:.3e}".format(RMSE_b1)+"\n"+r'$\eta$('+str(v_anchor)+' km/s) = '+"{0:.3e}".format(fixpoint)+"\n"+r'$\chi^2$ one slope = '+"{0:.3e}".format(chisquared_a1)+"\n"+r'$\chi^2$ two slopes = '+"{0:.3e}".format(chisquared_lb1))



plt.savefig('AllMLvc'+str(mi)+'_log_int_'+ str(today) +'z'+str(int(z))+'.pdf')
plt.clf()


fig2 =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):

	if (fillnaked):
		plt.plot(Mstar_ar[count], ML_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(Mstar_ar[count][fillcuts[count]], ML_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(Mstar_ar[count][empcuts[count]], ML_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(Mstar_ar[count], ML_ar[count], yerr=(ML_ar[count]*0.4, ML_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)


	count+=1

plt.plot(Mstarspace, Mstar1fit, color='k', ls= ':', lw =2)

plt.ylim(1e-1, 1e3)
plt.xlim(1e4, 1e11)
plt.xscale('log')
plt.yscale('log')
#plt.ylim(0, 10)

xlabel_str = ''
if (usemagic):
	xlabel_str = r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.xlabel(r'Stellar Mass $(M_{\odot})$ at z='+str(int(z)),fontsize=20)
plt.ylabel(r'$\eta$ = $\dot{M}_{out}$/SFR',fontsize=label_fontsize)
#plt.legend(thar, loc = 'best', prop={'size':10})

#text fun
#plt.text(10**9.0, 10**2.2, r'$\eta$ = $10^{%s}$'%"{0:.3g}".format(Mstar_int)+r'$\times  M_{star}^{%s}$'%"{0:.3g}".format(Mstar_slope)+"\n"+ r'$\chi^2$ = '+"{0:.3e}".format(chisquared_Mstar1) )


ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/192.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/192.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/192.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')

plt.savefig('AllMstar'+str(mi)+'_log_int_'+ str(today) +'z'+str(int(z))+'.pdf')
ff.close()

if (not noWV):
	fig4 =  plt.figure(figsize=(10,9))
	count = 0
	while (count < maxplot):
		if (fillnaked):
			plt.plot(v_ar[count], veloc_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		else:
			if (len(v_ar[count][empcuts[count]]) > 0):
				plt.plot(v_ar[count][fillcuts[count]], veloc_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
			if (len(v_ar[count][empcuts[count]]) > 0):
				if (not skip_empty_symbols):
					plt.plot(v_ar[count][empcuts[count]], veloc_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (Ulim_ar[count]):
			plt.errorbar(v_ar[count], veloc_ar[count], yerr=(veloc_ar[count]*0.4, veloc_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
		count+=1
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1e1, 3e2)
	plt.ylim(1e1, 1e3)
	veloc_fit = pow(xspace, veloc_slope) * 10**veloc_yint
	plt.plot(xspace, veloc_fit, ls=':', color ='k', lw=2)
	xlabel_str = ''
	if (usemagic):
		xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
	elif (usemed):
		xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
	else:
		xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))
	plt.xlabel(xlabel_str, fontsize=label_fontsize)

	if (velocity_mode_95):
		plt.ylabel(r'$95^{\rm th}$ percentile wind velocity (km/s)',fontsize=label_fontsize)
	else:
		plt.ylabel(r'median wind velocity (km/s)',fontsize=label_fontsize)
	ax = plt.gca()
	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(axis_fontsize)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(axis_fontsize)
		label.set_family('serif')
		#label.set_weight('bold')                                                                         
	ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
	ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
	ax.xaxis.set_tick_params(width=2, length=12, which='minor')
	ax.yaxis.set_tick_params(width=2, length=12, which='minor')
	
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	titlefont = 20
	plt.text(x_range[1]/3.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
	plt.text(x_range[1]/3.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
	plt.text(x_range[1]/3.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')
	if (velocity_mode_95):
		plt.savefig('95_veloc_vc'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')
	else:
		plt.savefig('med_veloc_vc'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')




ff.close()

plt.clf()
fig5b =  plt.figure(figsize=(10.5,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(vc_ar[count], outflowmet9_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(vc_ar[count], ZG9_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
	else:
		plt.plot(v_ar[count][fillcuts[count]], outflowmet9_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(v_ar[count][empcuts[count]], outflowmet9_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(v_ar[count][fillcuts[count]], ZG9_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(v_ar[count][empcuts[count]], ZG9_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(v_ar[count], outflowmet9_ar[count], yerr=(outflowmet9_ar[count]*0.4, outflowmet9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
		#plt.errorbar(v_ar[count], outflowmet9_ar[count], yerr=(outflowmet9_ar[count]*0.2, outflowmet9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
	count+=1
plt.xscale('log')
plt.yscale('log')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
plt.xlim(1e1, 3e2)
plt.ylim(5e-5, 1e-1)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
else:
	xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))

plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.ylabel(r'$\eta_z$ = $\dot{M}_{\rm out}Z_{\rm out}$/SFR (r=25 kpc)',fontsize=label_fontsize)
plt.ylabel(r'$\eta_z$ = $\dot{M}_{\rm out}Z_{\rm out}$/SFR '+Outer_label,fontsize=label_fontsize)

#plt.ylabel(r'$\eta_z$ = $\dot{M}_{\rm out}Z_{\rm out}$/SFR (r=$R_{\rm vir}$)',fontsize=label_fontsize)
#plt.ylabel(r'$\eta (at 1.0 $R_{\rm vir}$)',fontsize=label_fontsize)

ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[0]*1.2, y_range[0]*6.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[0]*1.2, y_range[0]*3.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[0]*1.2, y_range[0]*1.5,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')
xfarag = np.linspace(5e0, 3.5e2, num=100)
yfarag = xfarag*0 + 0.02
plt.plot(xfarag, yfarag, ls=':', color = 'k', lw=2)

plt.savefig('metal9_vc'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig5 =  plt.figure(figsize=(10.5,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(vc_ar[count], outflowmet_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(vc_ar[count], ZG_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)

	else:
		plt.plot(v_ar[count][fillcuts[count]], outflowmet_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(v_ar[count][empcuts[count]], outflowmet_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(v_ar[count][fillcuts[count]], ZG_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(v_ar[count][empcuts[count]], ZG_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(v_ar[count], outflowmet_ar[count], yerr=(outflowmet_ar[count]*0.4, outflowmet_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e1, 3e2)
plt.ylim(5e-5, 1e-1)

#xfarag = plt.linspace(5e0, 3.5e2, 100)
#yfarag = xfarag + 0.02
plt.plot(xfarag, yfarag, ls=':', color = 'k', lw=2)


#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
else:
	xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.ylabel(r'<$\eta \times Z_{\rm out}$> (at 0.25 $R_{\rm vir}$)',fontsize=label_fontsize)

plt.ylabel(r'$\eta_z$ = $\dot{M}_{\rm out}Z_{\rm out}$/SFR '+Inner_label,fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[0]*1.2, y_range[0]*6.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[0]*1.2, y_range[0]*3.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[0]*1.2, y_range[0]*1.5,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('metal_vc'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')
plt.clf()
fig6 =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(ISMmet_ar[count], metrat_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(ISMmet_ar[count], metratG_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
		
	else:
		plt.plot(ISMmet_ar[count][fillcuts[count]], metrat_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(ISMmet_ar[count][empcuts[count]], metrat_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(ISMmet_ar[count][fillcuts[count]], metratG_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(ISMmet_ar[count][empcuts[count]], metratG_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
		if (Ulim_ar[count]):
			plt.errorbar(ISMmet_ar[count], metrat_ar[count], yerr=(metrat_ar[count]*0.4, metrat_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
			if (plot_inflow):
				plt.errorbar(ISMmet_ar[count], metratG_ar[count], yerr=(metratG_ar[count]*0.4, metratG_ar[count]*0), marker='_', ms=magicmarker, color = 'g', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('log')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'$Z_{\rm ISM}$'
elif (usemed):
	xlabel_str =r'ISM metallicity'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
xlabel_str = r'$Z_{\rm ISM}$'
plt.xlabel(xlabel_str, fontsize=label_fontsize)
plt.ylabel(r'$Z_{\rm out}$ '+Inner_label,fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/24.0, y_range[0]*6.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/24.0, y_range[0]*3.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/24.0, y_range[0]*1.5,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


plt.savefig('metalrat_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig6b =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(ISMmet_ar[count], metrat9_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(ISMmet_ar[count], metratG9_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
	else:
		plt.plot(ISMmet_ar[count][fillcuts[count]], metrat9_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(ISMmet_ar[count][empcuts[count]], metrat9_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(ISMmet_ar[count][fillcuts[count]], metratG9_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(ISMmet_ar[count][empcuts[count]], metratG9_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(ISMmet_ar[count], metrat9_ar[count], yerr=(metrat9_ar[count]*0.4, metrat9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
		if (plot_inflow):
			plt.errorbar(ISMmet_ar[count], metratG9_ar[count], yerr=(metratG9_ar[count]*0.4, metratG9_ar[count]*0), marker='_', ms=magicmarker, color = 'g', ls='none', lolims=True)
	count+=1
plt.xscale('log')
plt.yscale('log')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))

xlabel_str = r'$Z_{\rm ISM}$'
plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm out}$ '+Outer_label,fontsize=label_fontsize)

#plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.ylabel(r'$Z_{\rm out} (at Rvir)$',fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/24.0, y_range[0]*6.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/24.0, y_range[0]*3.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/24.0, y_range[0]*1.5,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('metalrat9_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig7 =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(vc_ar[count], metrat_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(vc_ar[count], metratG_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)

	else:
		plt.plot(v_ar[count][fillcuts[count]], metrat_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(v_ar[count][empcuts[count]], metrat_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(v_ar[count][fillcuts[count]], metratG_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(v_ar[count][empcuts[count]], metratG_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
		if (Ulim_ar[count]):
			plt.errorbar(v_ar[count], metrat_ar[count], yerr=(metrat_ar[count]*0.4, metrat_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
			if (plot_inflow):
				plt.errorbar(v_ar[count], metratG_ar[count], yerr=(metratG_ar[count]*0.2, metratG_ar[count]*0), marker='_', ms=magicmarker, color = 'g', ls='none', lolims=True)
	count+=1
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
else:
	xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm out}$ '+Inner_label,fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/3.0, y_range[0]*6.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/3.0, y_range[0]*3.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/3.0, y_range[0]*1.5,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


plt.savefig('metalrat_vc'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig7b =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(vc_ar[count], metrat9_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(v_ar[count][fillcuts[count]], metrat9_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(v_ar[count][empcuts[count]], metrat9_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(v_ar[count][fillcuts[count]], metratG9_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(v_ar[count][empcuts[count]], metratG9_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(v_ar[count], metrat9_ar[count], yerr=(metrat9_ar[count]*0.4, metrat9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
		if (plot_inflow):
			plt.errorbar(v_ar[count], metratG9_ar[count], yerr=(metratG9_ar[count]*0.2, metratG9_ar[count]*0), marker='_', ms=magicmarker, color = 'g', ls='none', lolims=True)
	count+=1
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
else:
	xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.ylabel(r'$Z_{\rm out}$ (r=$R_{\rm vir}$)',fontsize=label_fontsize)
plt.ylabel(r'$Z_{\rm out}$ '+Outer_label,fontsize=label_fontsize)

ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/3.0, y_range[0]*6.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/3.0, y_range[0]*3.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/3.0, y_range[0]*1.5,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('metalrat9_vc'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')


plt.clf()
fig8 =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(ISMmet_ar[count], (metrat9_ar[count]/metratG9_ar[count]), marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(ISMmet_ar[count][fillcuts[count]], metrat9_ar[count][fillcuts[count]]/metratG9_ar[count][fillcuts[count]] , marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(ISMmet_ar[count][empcuts[count]], metrat9_ar[count][empcuts[count]]/metratG9_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker+5)
	if (Ulim_ar[count]):
		plt.errorbar(ISMmet_ar[count], metrat9_ar[count]/metratG9_ar[count], yerr=((metrat9_ar[count]/metratG9_ar[count])*0.4, (metrat9_ar[count]/metratG9_ar[count])*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)
#		print 'sasha duck duck duck ',ISMmet_ar[count],metrat9_ar[count]
	count+=1
plt.xscale('log')
plt.yscale('log')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
#plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
plt.ylim(1e-1, 1e2)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
xlabel_str = r'$Z_{\rm ISM}$'

plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm out} / Z_{\rm in}$ '+Outer_label,fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/24.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/24.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/24.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


plt.savefig('outflowZ_ISM9'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig8b =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(ISMmet_ar[count], (metrat_ar[count]/metratG_ar[count]), marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(ISMmet_ar[count][fillcuts[count]], metrat_ar[count][fillcuts[count]]/metratG_ar[count][fillcuts[count]] , marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(ISMmet_ar[count][empcuts[count]], metrat_ar[count][empcuts[count]]/metratG_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(ISMmet_ar[count], metrat_ar[count]/metratG_ar[count], yerr=((metrat_ar[count]/metratG_ar[count])*0.4, (metrat_ar[count]/metratG_ar[count])*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('log')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
#plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
plt.ylim(1e-1, 1e2)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
xlabel_str = r'$Z_{\rm ISM}$'

plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm out} / Z_{\rm in}$ '+Inner_label,fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/24.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/24.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/24.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('outflowZ_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig9 =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(Mstar_ar[count], (metrat_ar[count]/ISMmet_ar[count]), marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(Mstar_ar[count][fillcuts[count]], metrat_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(Mstar_ar[count][empcuts[count]], metrat_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		if (plot_inflow):
			plt.plot(Mstar_ar[count][fillcuts[count]], metratG_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
			if (not skip_empty_symbols):
				plt.plot(Mstar_ar[count][empcuts[count]], metratG_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)

	if (Ulim_ar[count]):
		plt.errorbar(Mstar_ar[count], metrat_ar[count]/ISMmet_ar[count], yerr=((metrat_ar[count]/ISMmet_ar[count])*0.4, ML9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('linear')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
#plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
plt.ylim(0, 2)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
#xlabel_str = r'$Z_{\rm ISM}$'

plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm out}$ ' +Inner_label +r' / $Z_{\rm ISM}$',fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/100.0, (y_range[1]-y_range[0])*0.90+y_range[0],   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.83+y_range[0],   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.76+y_range[0],   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


plt.savefig('MstarZ_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()
fig9b =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(Mstar_ar[count], (metrat_ar[count]/ISMmet_ar[count]), marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(Mstar_ar[count][fillcuts[count]], metrat9_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(Mstar_ar[count][empcuts[count]], metrat9_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		#plt.plot(Mstar_ar[count][fillcuts[count]], metratG9_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
		#plt.plot(Mstar_ar[count][empcuts[count]], metratG9_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(Mstar_ar[count], metrat9_ar[count]/ISMmet_ar[count], yerr=((metrat9_ar[count]/ISMmet_ar[count])*0.4, ML9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('linear')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
#plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
plt.ylim(0,2)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
#xlabel_str = r'$Z_{\rm ISM}$'

plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm out}$ ' +Outer_label+r' / $Z_{\rm ISM}$',fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/100.0, (y_range[1]-y_range[0])*0.90+y_range[0],   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.83+y_range[0],   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.76+y_range[0],   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('MstarZ9_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()

fig9c =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(Mstar_ar[count], (metrat_ar[count]/ISMmet_ar[count]), marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		#plt.plot(Mstar_ar[count][fillcuts[count]], metrat9_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		#plt.plot(Mstar_ar[count][empcuts[count]], metrat9_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		plt.plot(Mstar_ar[count][fillcuts[count]], metratG9_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count],  color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(Mstar_ar[count][empcuts[count]], metratG9_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count],  color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(Mstar_ar[count], metratG9_ar[count]/ISMmet_ar[count], yerr=((metratG9_ar[count]/ISMmet_ar[count])*0.4, ML9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('linear')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace
#plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
plt.ylim(0,1.0)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
#xlabel_str = r'$Z_{\rm ISM}$'

plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm in}$ ' +Outer_label+r' / $Z_{\rm ISM}$',fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/100.0, (y_range[1]-y_range[0])*0.90+y_range[0],   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.83+y_range[0],   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.76+y_range[0],   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('inflowMstarZ9_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()


fig9d =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):
	if (fillnaked):
		plt.plot(Mstar_ar[count], (metrat_ar[count]/ISMmet_ar[count]), marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		#plt.plot(Mstar_ar[count][fillcuts[count]], metrat9_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		#plt.plot(Mstar_ar[count][empcuts[count]], metrat9_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		plt.plot(Mstar_ar[count][fillcuts[count]], metratG_ar[count][fillcuts[count]]/ISMmet_ar[count][fillcuts[count]] , marker=marker_ar[count],  color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(Mstar_ar[count][empcuts[count]], metratG_ar[count][empcuts[count]]/ISMmet_ar[count][empcuts[count]], marker=marker_ar[count],  color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(Mstar_ar[count], metratG_ar[count]/ISMmet_ar[count], yerr=((metratG_ar[count]/ISMmet_ar[count])*0.4, ML9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1
plt.xscale('log')
plt.yscale('linear')
xspace = np.arange(-5, -2, 0.1)
xspace = 10**xspace
yspace = xspace

#plt.plot(xspace,yspace, ls=':')
#plt.xlim(1e1, 3e2)
#plt.ylim(1e1, 1e3)
plt.ylim(0,1.0)
xlabel_str = ''
if (usemagic):
	xlabel_str = r'ISM metallicity'
elif (usemed):
	xlabel_str =r'ISM metallicity'
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'ISM metallicity'+str(int(z))
#xlabel_str = r'$Z_{\rm ISM}$'

plt.xlabel(xlabel_str, fontsize=label_fontsize)

plt.ylabel(r'$Z_{\rm in}$ '+Inner_label+r' / $Z_{\rm ISM}$',fontsize=label_fontsize)
ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')

x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/100.0, (y_range[1]-y_range[0])*0.90+y_range[0],   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.83+y_range[0],   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/100.0,  (y_range[1]-y_range[0])*0.76+y_range[0],   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')



plt.savefig('inflowMstarZ_ISM'+str(mi)+'int'+ str(today) +'z'+str(int(z))+'.pdf')

plt.clf()





fig =  plt.figure(figsize=(10,9))

count = 0
while (count < maxplot):

	if (fillnaked):
		plt.plot(v_ar[count], ML9_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		#plt.plot(v_ar[count], MG_ar[count], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)

	else:
		plt.plot(v_ar[count][fillcuts[count]], ML9_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(v_ar[count][empcuts[count]], ML9_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
		#plt.plot(v_ar[count][fillcuts[count]], MG_ar[count][fillcuts[count]], marker=marker_ar[count], color='g', ls='none', ms=magicmarker)
		#plt.plot(v_ar[count][empcuts[count]], MG_ar[count][empcuts[count]], marker=marker_ar[count], color='g', ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(v_ar[count], ML9_ar[count], yerr=(ML9_ar[count]*0.4, ML9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)

	count+=1

#plt.plot(v_ar[0], ML_ar[0], 'Dk')
#plt.plot(v_ar[1], ML_ar[1], 'Db')
#plt.plot(v_ar[2], ML_ar[2], 'sk')
#plt.plot(v_ar[3], ML_ar[3], 'sb')
#plt.plot(v_ar[4], ML_ar[4], '^k')
#plt.plot(xspace, v1fit, color='k', ls='--', markevery=1, lw=2)
#plt.plot(xspace, v2fit, color='k', ls='--', markevery=1, lw=2)
#if (show_single_fit_vc):
#	plt.plot(xspace, a1fit, color='k', ls= '--', lw =2)
#if (show_double_fit_vc):
#	plt.plot(l1space, l1fit, color='k', ls=':', lw=2)
#	plt.plot(b1space, b1fit, color='k', ls =':', lw=2)
#	plt.plot (l1space, l1fit_medz, color ='b', ls=':', lw=2)
#	plt.plot (b1space, b1fit_medz, color ='b', ls=':', lw=2)
#	plt.plot (l1space, l1fit_loz, color ='r', ls=':', lw=2)
#	plt.plot (b1space, b1fit_loz, color ='r', ls=':', lw=2)	

#plt.plot(v_ar[7], ML_ar[7], 'hb')

plt.xlim(1e1, 3e2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.xscale('log')

ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')
x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/3.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/3.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/3.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


xlabel_str = ''
if (usemagic):
	xlabel_str = r'$V_c$ (km/s) at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'$V_c$ (km/s) at ${\rm z}_{med}$'
else:
	xlabel_str =r'$V_c$ (km/s) at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.xlabel(r'$V_c$ (km/s) at z='+str(int(z)), fontsize=20)
plt.ylabel(r'$\eta$ = $\dot{M}_{out}$'+Outer_label+'/SFR',fontsize=label_fontsize)
if (dotoplabel):
	plt.title('z='+str(int(z)), fontsize=label_fontsize)
#plt.legend(thar, loc = 'best', prop={'size':10})

#text fun
#plt.text(10**1.5, 10**2.2, "v < "+str(v_anchor)+" km/s slope = "+"{0:.3e}".format(little_slope)+r'$\pm$ '+"{0:.3e}".format(RMSE_l1)+"\n v > "+str(v_anchor)+" km/s slope = "+"{0:.3e}".format(big_slope)+r'$\pm$ '+"{0:.3e}".format(RMSE_b1)+"\n"+r'$\eta$('+str(v_anchor)+' km/s) = '+"{0:.3e}".format(fixpoint)+"\n"+r'$\chi^2$ one slope = '+"{0:.3e}".format(chisquared_a1)+"\n"+r'$\chi^2$ two slopes = '+"{0:.3e}".format(chisquared_lb1))



plt.savefig('AllMLvc9'+str(mi)+'_log_int_'+ str(today) +'z'+str(int(z))+'.pdf')
plt.clf()


plt.clf()


fig2 =  plt.figure(figsize=(10,9))
count = 0
while (count < maxplot):

	if (fillnaked):
		plt.plot(Mstar_ar[count], ML9_ar[count], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
	else:
		plt.plot(Mstar_ar[count][fillcuts[count]], ML9_ar[count][fillcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', ms=magicmarker)
		if (not skip_empty_symbols):
			plt.plot(Mstar_ar[count][empcuts[count]], ML9_ar[count][empcuts[count]], marker=marker_ar[count], color=color_ar[count], ls='none', fillstyle='none', ms=magicmarker)
	if (Ulim_ar[count]):
		plt.errorbar(Mstar_ar[count], ML9_ar[count], yerr=(ML9_ar[count]*0.4, ML9_ar[count]*0), marker='_', ms=magicmarker, color = 'r', ls='none', lolims=True)


	count+=1


plt.ylim(1e-1, 1e3)
plt.xlim(1e4, 1e11)
plt.xscale('log')
plt.yscale('log')
#plt.ylim(0, 10)

xlabel_str = ''
if (usemagic):
	xlabel_str = r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{out,half}$'
elif (usemed):
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at ${\rm z}_{med}$'
else:
	xlabel_str =r'Stellar Mass $(M_{\odot})$ at z='+str(int(z))
plt.xlabel(xlabel_str, fontsize=label_fontsize)

#plt.xlabel(r'Stellar Mass $(M_{\odot})$ at z='+str(int(z)),fontsize=20)
plt.ylabel(r'$\eta$ = $\dot{M}_{out}$' +Outer_label+'/SFR',fontsize=label_fontsize)
#plt.legend(thar, loc = 'best', prop={'size':10})

#text fun
#plt.text(10**9.0, 10**2.2, r'$\eta$ = $10^{%s}$'%"{0:.3g}".format(Mstar_int)+r'$\times  M_{star}^{%s}$'%"{0:.3g}".format(Mstar_slope)+"\n"+ r'$\chi^2$ = '+"{0:.3e}".format(chisquared_Mstar1) )


ax = plt.gca()
ticklabels = ax.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
ticklabels = ax.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(axis_fontsize)
    label.set_family('serif')
    #label.set_weight('bold')                                                                         
ax.xaxis.set_tick_params(width=tick_width, length=tick_height)
ax.yaxis.set_tick_params(width=tick_width, length=tick_height)
ax.xaxis.set_tick_params(width=2, length=12, which='minor')
ax.yaxis.set_tick_params(width=2, length=12, which='minor')
x_range = ax.get_xlim()
y_range = ax.get_ylim()

titlefont = 20
plt.text(x_range[1]/192.0, y_range[1]/3.0,   "4.0 > z > 2.0", fontsize=titlefont, fontweight='bold', color ='k')
plt.text(x_range[1]/192.0, y_range[1]/6.0,   "2.0 > z > 0.5", fontsize=titlefont, fontweight='bold', color ='b')
plt.text(x_range[1]/192.0, y_range[1]/12.0,   "0.5 > z", fontsize=titlefont, fontweight='bold', color ='r')


plt.savefig('AllMstar9'+str(mi)+'_log_int_'+ str(today) +'z'+str(int(z))+'.pdf')
ff.close()


metrat_ar_ext = np.log10(np.hstack(metrat_ar))
fillcut_ar_ext = np.hstack(fillcuts)
ISMmet_ar_ext = np.log10(np.hstack(ISMmet_ar))
metrat9_ar_ext = np.log10(np.hstack(metrat9_ar))

print 'asdf asdf ',metrat_ar
print 'dfsa dfsa', metrat_ar_ext
print 'dofa dofa' ,ISMmet_ar_ext

garthcut = np.isfinite(ISMmet_ar_ext)
print 'dofa dofa' ,ISMmet_ar_ext[garthcut]

p3 = [0, 1, -1]


ISM_Outfit = optimize.leastsq(residualsZ,p3[:], args=(metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut], ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut] ))


ISMmet_int = ISM_Outfit[0][0]
ISMmet_slope = ISM_Outfit[0][1]
ISMmet_z_slope = ISM_Outfit[0][2]

print 'heres ISM metrat fit ',ISMmet_int,ISMmet_slope,ISMmet_z_slope

print 'noz'
ISM_Outfit = optimize.leastsq(residuals,p2[:], args=(metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut], ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut] ))


ISMmet_int = ISM_Outfit[0][0]
ISMmet_slope = ISM_Outfit[0][1]

print 'heres ISM metrat fit ',ISMmet_int,ISMmet_slope

print 'hiz only'
ISM_Outfit = optimize.leastsq(residuals,p2[:], args=(metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut * hiz_cut], ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut * hiz_cut] ), full_output=True)

ISMmet_int = ISM_Outfit[0][0]
ISMmet_slope = ISM_Outfit[0][1]
print 'heres ISM metrat fit ',ISMmet_int,ISMmet_slope


print 'proper z fit at 0.25 Rvir' 
ISM_Outfit = optimize.leastsq(residualsZZ,p3[:], args=(metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut], ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut] ), full_output=True)


ISMmet_int = ISM_Outfit[0][0]
ISMmet_slope = ISM_Outfit[0][1]
ISMmet_z_slope = ISM_Outfit[0][2]

ISM_Outfitcov = ISM_Outfit[1]
#print 'ISM_Outfit cov ',ISM_Outfitcov
perrS = np.sqrt(np.diag(ISM_Outfitcov))
#print 'error is here ', perrS #TSIS IS WRONG NEED TO MULTIPLY BY VARIANCE!


print 'heres ISM metrat fit ',ISMmet_int,ISMmet_slope,ISMmet_z_slope

print 'proper z fit at 1.0 Rvir' 
ISM_Outfit = optimize.leastsq(residualsZZ,p3[:], args=(metrat9_ar_ext[stoopcut * fillcut_ar_ext * garthcut], ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut] ), full_output=True)


ISMmet_int = ISM_Outfit[0][0]
ISMmet_slope = ISM_Outfit[0][1]
ISMmet_z_slope = ISM_Outfit[0][2]

ISM_Outfitcov = ISM_Outfit[1]
#print 'ISM_Outfit cov ',ISM_Outfitcov
perrS = np.sqrt(np.diag(ISM_Outfitcov))
#print 'error is here ', perrS TSIS IS WRONG NEED TO MULTIPLY BY VARIANCE!


print 'heres ISM metrat fit at outer radius',ISMmet_int,ISMmet_slope,ISMmet_z_slope


#def residualsZZ(b , y, x, zz):
#    err = y - (b[1]*x + b[0] + b[2]*np.log10(1+zz))
#    return err


def fitfuncStar((x, myz),a1,a2, a3):
	#fixpoint already assumes z=3
	y = a1*x + a2 + a3*np.log10(1.0+myz) 
	
	return y

popt,pcov = optimize.curve_fit(fitfuncStar, ( ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut]),  metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut], p0=(-1.0,1.0, 0.0)) 

print 'here we go again popt for 0.25 Rvir', popt
[a1, a2, a3]  = popt
print  'parameters ', a1, a2, a3
perrS = np.sqrt(np.diag(pcov))
print 'errors ',perrS




popt2,pcov2 = optimize.curve_fit(fitfuncStar, ( ISMmet_ar_ext[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut]),  metrat9_ar_ext[stoopcut * fillcut_ar_ext * garthcut], p0=(-1.0,1.0, 0.0)) 

print 'here we go again popt for 1.0 Rvir', popt2
[a1, a2, a3]  = popt2
print 'parameters ',a1, a2, a3
perrS2 = np.sqrt(np.diag(pcov2))
print 'errors ',perrS2


popt2,pcov2 = optimize.curve_fit(fitfuncStar, ( v_ext_ar[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut]),  metrat9_ar_ext[stoopcut * fillcut_ar_ext * garthcut], p0=(-1.0,1.0, 0.0)) 

print 'here we go again v_c popt for 1.0 Rvir', popt2
[a1, a2, a3]  = popt2
print 'parameters ',a1, a2, a3
perrS2 = np.sqrt(np.diag(pcov2))
print 'errors ',perrS2


popt2,pcov2 = optimize.curve_fit(fitfuncStar, ( v_ext_ar[stoopcut * fillcut_ar_ext * garthcut], z_rel_ext_ar[stoopcut*fillcut_ar_ext * garthcut]),  metrat_ar_ext[stoopcut * fillcut_ar_ext * garthcut], p0=(-1.0,1.0, 0.0)) 

print 'here we go again v_c popt for 0.25 Rvir', popt2
[a1, a2, a3]  = popt2
print 'parameters ',a1, a2, a3
perrS2 = np.sqrt(np.diag(pcov2))
print 'errors ',perrS2

