import numpy as np
import matplotlib.pyplot as plt
import sasha_simple_plots as ssp 
from datetime import date
import Sasha_functions as SF
today = date.today()
print today

unnormed = False
area_of_cell = 0.3*0.3

finname_ar = []
f = open('list_of_files.txt', 'r')
dars = f.readlines()
for line in dars:
	xsd = line.split()
	finname_ar.append(str(xsd[0]))
f.close()	

mass_ar = []
vesc_ar = []
sig_ar = []
mass_ext_ar = []
sig_ext_ar = []
vesc_ext_ar = []
temp_mass_ar = [] 
temp_vesc_ar = []
temp_sig_ar = []
ep_name_ar = []
duration_ar = []
Neps_ar = []
SFinterval_ar = []
maxSF_ar = []
totMass_ar = []
maxsig_ar = []
maxweightsig_ar = []

for name in finname_ar:
	print name
	if ((name!='break') and (name!='')):
		f = open(name)
		header = f.readline()
		xsd = header.split()
		SFR_time_range = float(xsd[6])
		print SFR_time_range
		SFinterval_ar.append(SFR_time_range)
		f.close()

validcount = 0
lastname = ''
Neps = 0
duration = 0
maxSF = 0
totMass = 0
maxsig = 0
maxweightedsig = 0

for name in finname_ar:
	print name
	if ((name!='break') and (name!='')):
		f = open(name)
		dars = np.loadtxt(f, ndmin=2)
		mass = dars[:,1]
		vesc = dars[:,2] 
		duration +=  SFinterval_ar[validcount]
		totMass += np.sum(mass)
		totSFR = np.sum(mass)/ SFinterval_ar[validcount]
		if (totSFR > maxSF):
			maxSF = totSFR
		temp_sig_ar.extend(mass/ SFinterval_ar[validcount])
		sig = mass/ (SFinterval_ar[validcount])
		weighted_avg_sig = np.average(sig, weights=sig)/area_of_cell
		print 'weighted average sig ',weighted_avg_sig
		unweighted_avg_sig = np.mean(sig_ar)/area_of_cell
		print 'unweighted average sig ',unweighted_avg_sig		
		if (weighted_avg_sig > maxweightsig):
			maxweightedsig = weighted_avg_sig	
		if (unweighted_avg_sig > maxsig):
			maxsig = unweighted_avg_sig	
		temp_mass_ar.extend(mass)
		temp_vesc_ar.extend(vesc)
		validcount+=1
		Neps+=1 
		lastname = name
		f.close()
	else:
		name_of_ep = lastname.split('/')[0]
		ep_name_ar.append(lastname)
		mass_ar.append(temp_mass_ar)
		vesc_ar.append(temp_vesc_ar)
		sig_ar.append(np.array(temp_sig_ar))
		mass_ext_ar.extend(temp_mass_ar)
		vesc_ext_ar.extend(temp_vesc_ar)
		sig_ext_ar.extend(temp_sig_ar)
		Neps_ar.append(Neps)
		duration_ar.append(duration)
		maxSF_ar.append(maxSF)
		totMass_ar.append(totMass)
		maxsig_ar.append(maxsig)
		maxweightedsig_ar.append(maxweightedsig)

		maxSF = 0
		duration = 0		
		totMass = 0
		Neps = 0
		temp_mass_ar = []
		temp_vesc_ar = []
		temp_sig_ar = []

mass_ar = np.array(mass_ar)
vesc_ar = np.array(vesc_ar)
sig_ar = np.array(sig_ar)



figure0 = plt.figure(figsize=(10, 8))
plt.xlabel('escape velocity ')
plt.ylabel('stellar mass formed ')

num_bins = 100

#ep 253 - B1, z=1.36, Mstar =8.5e+09 , eta=7, sigsfr=5.7
#ep 300 - B1, z= 0.9, Mstar = 9.9e9, eta = 3, sigsfr = 32
#ep 318 - B1, z = 0.72, Mstar = 1.1e10, eta=1.5, sigsfr=18
#ep 340 - B1, z=0.5, Mstar = 1.3e10 , eta= 0.1, sigsfr=19

my_special_legend = [' B1, z=1.36, Mstar =8.5e+09 , eta=7, sigsfr=4', 'B1, z= 0.9, Mstar = 9.9e9, eta = 3, sigsfr = 3', 
' z = 0.72, Mstar = 1.1e10, eta=1.5, sigsfr=3', ' z=0.5, Mstar = 1.3e10 , eta= 0.3, sigsfr=3', 'dm506, z=2.16, Mstar=1.1e+10, eta=5, sigsfr=7']

foutname = str(today)+'selected_eps.txt'
f = open(foutname, 'a')
color_ar = ['k', 'r', 'b', 'g', 'm', 'y', 'DarkSlateBlue', 'LightGreen']
color_ar.extend(color_ar)
color_ar.extend(color_ar)
ls_ar = [':', ':', ':', '-', '-', '-', '-', '-']
ls_ar.extend(ls_ar)
ls_ar.extend(ls_ar)

count = 0
for mass in mass_ar:
	totalmass = np.sum(mass)
	if (unnormed):
		print totalmass
		totalmass =1.0
	sig_squared = sig_ar[count] * sig_ar[count]
	sig_square_sum = np.sum(sig_squared)
	sig_sum = np.sum(sig_ar[count])
	average_weighted_sig = sig_square_sum/sig_sum
	
	vesc_ar[count] = np.array(vesc_ar[count])
	sig_ar[count] = np.array(sig_ar[count])
	vesc_order = np.argsort(vesc_ar[count])
	
	mass_sorted_cumsum = np.cumsum(sig_ar[count][vesc_order])
	mass_sorted_cumsum /= mass_sorted_cumsum[-1]
	vesc_med = np.interp(0.5, mass_sorted_cumsum,vesc_ar[count][vesc_order] )
	vesc_25 = np.interp(0.25, mass_sorted_cumsum,vesc_ar[count][vesc_order] )
	vesc_75 = np.interp(0.75, mass_sorted_cumsum,vesc_ar[count][vesc_order] )
	print 'name of ep ',ep_name_ar[count]
	print 'num of snaps in ep', Neps_ar[count]
	print 'vesc ',len(vesc_ar[count])
	print 'sig ',len(sig_ar[count])
	weight_avg_vesc = np.average(vesc_ar[count], weights=sig_ar[count])
	print 'weighted average vesc ',weight_avg_vesc
	print '25th, 50th, 75th percentile vesc ',vesc_25, vesc_med, vesc_75
	weighted_avg_sig = np.average(sig_ar[count], weights=sig_ar[count])/area_of_cell
	print 'weighted average sig ',weighted_avg_sig
	unweighted_avg_sig = np.mean(sig_ar[count])/area_of_cell
	print 'unweighted average sig ',unweighted_avg_sig
	
	myline = [count, ep_name_ar[count], Neps_ar[count], weight_avg_vesc, vesc_25, vesc_med, vesc_75, weighted_avg_sig, unweighted_avg_sig, maxSF_ar[count], totMass_ar[count], duration_ar[count] ]
	mystring = SF.line_to_string(myline)
	f.write(mystring)
	print count, average_weighted_sig
	shell_histresults = np.histogram(vesc_ar[count], num_bins, range=(0, 1000), weights=np.array(mass_ar[count])/totalmass)
	[histx, histy] = ssp.create_plottable_hist(shell_histresults)
	plt.plot(histx,histy, color=color_ar[count], ls = ls_ar[count], lw=2)
	count+=1
f.close()
plt.legend(my_special_legend, loc='best')
plt.savefig(str(today) +'normedeps'+'.pdf')
