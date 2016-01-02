import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import date
import sys

adaptive_filter_rat = 1e-10
delay = 0
SFR_filter = -1
use_flow = True
kpc25mode = True
doplots = False
#percentile = 90

if (len(sys.argv) < 6):
	print 'syntax: blah.py Nhalo outflowsfile vstats-base percentile zstart zend [extra percentile]'
	sys.exit()
percentile = str(sys.argv[4])
name = str(sys.argv[7])
if (len(sys.argv) > 8):
	print 'modifying percentile '
	percentileA = percentile+str(sys.argv[8])
	print percentileA
	#name = str(sys.argv[7])
else:
	percentileA = percentile

finname = name+'/tasz.txt'
f = open(finname)
dars = np.loadtxt(f)
theas = dars[:,0]
thets = dars[:,1]
f.close()
tfora = interp1d(theas, thets)

haloN = str(sys.argv[1])

#outflows_finname = name+'/'+str(sys.argv[2])
outflows_finname = str(sys.argv[2])


f = open(outflows_finname)
dars = np.loadtxt(f)
f.close()
Zs = dars[:,1]
As = 1.0 / (1.0 + Zs)
Ts = tfora(As)
Ts = Ts * 1000.0

Ns = dars[:,0]
TotalMass = dars[:,6]
NewSFR_timescale = dars[:,11]
SFRInner = dars[:,3]

zstart = 4.0
zend = 2.0
zstart = float(sys.argv[5])
zend = float(sys.argv[6])
little_h = 0.7

cut1 = Zs >= zend
cut2 = Zs <= zstart
cuts = cut1*cut2

Zs = Zs[cuts]
As = As[cuts]
Ts = Ts[cuts]
Ns = Ns[cuts]
SFRInner = SFRInner[cuts]
TotalMass = TotalMass[cuts]
NewSFR_timescale = NewSFR_timescale[cuts]
getSFR =  interp1d(Ts, SFRInner, kind='linear')


f = open(outflows_finname)
dars_string =np.loadtxt(f, dtype='str')
Nlist = dars_string[:,0][cuts]

count = 0
for theN in Nlist:
	if (int(theN) < 100):
		Nlist[count] = '0'+Nlist[count]
	count += 1

Nsnaps = len(Nlist)

f.close()



NumShells = 20
Flow_Shells = []
Cross_Shells = []




starter = -5 
count = 0
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow = dars[:,riteindex][cuts]
	#print 'The Flow of Shell ', count,' is index ',riteindex
	Flow_Shells.append(TheFlow)
	count += 1


starter = -1 
count = 0
while (count < NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter
	TheCross = (dars[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	#print 'The Cross of Shell ', count,' is index ',riteindex
	Cross_Shells.append(TheCross)
	count+=1

Flow_in = Cross_Shells[2]
if (use_flow):
	Flow_in = Flow_Shells[2]
Outflow_filter = np.max(Flow_in) * adaptive_filter_rat 


vbase = str(sys.argv[3])
#vsup = str(sys.argv[4])

#vsup = name+'/results-p'+percentileA+'/vtable_N'+haloN+'p'+percentile+'.txt'
if (kpc25mode):
	vsup = name+'/results25kpc-p'+percentileA+'/vtable_N'+haloN+'p'+percentile+'.txt'
else:
	vsup = name+'/results-p'+percentileA+'/vtable_N'+haloN+'p'+percentile+'.txt'


ff = open(vsup)
supdars = np.loadtxt(ff)
supvels = supdars[:,1]
supzs = supdars[:,4]
supcuts = (supzs <= zstart) * (supzs >= zend) 
supvels = supvels[supcuts]
ff.close()




if (len(supvels) != len(Cross_Shells[9])):
	print 'your files arent the same lengths ',len(supvels),len(Cross_Shells[9])
	sys.exit()

count = 0

med_ar = []
max_ar = []
bulk_ar = []
time_ar = [] 
ar90 = []


Stable = []
SFRtable = []
Tbegin = Ts[0]

while(count < len(Nlist)):
	epoch = As[count]
	thetime = tfora(epoch)
	theoutflow = Flow_in[count]
	time_ar.append(thetime*1000)
	print 'doing ',Nlist[count]
	print 'epoch ',epoch
	Nsnapstring = Nlist[count]
	finname1=vbase+Nsnapstring+'N'+haloN+'.txt'	
	
	ShellUseInner = '2'
	ShellUseOuter = '9'
	#finnames = [finname1, finname2, finname3]

	f = open(finname1)
	dars = np.loadtxt(f)
	Rs = dars[:,0]
	Med_out = dars[:,3]
	Max_out = dars[:,5]
	v90 = dars[:,7]
	Bulk = dars[:,2]
	f.close()
	
	f = open(finname1)
	header = f.readline()
	xsd = header.split()
	f.close()

	Mhalo = float(xsd[7])

	print 'median ',Med_out[3], supvels[count], theoutflow
	
	med_ar.append(Med_out[3])
	max_ar.append(np.max(Max_out))
	bulk_ar.append(Bulk[3])
	#ar90.append(v90[3])
	ar90.append(supvels[count])
	
	if ((Ts[count] - delay) > Tbegin):
		theSFR = getSFR(Ts[count])
	else:
		theSFR = 0

	if (theoutflow > Outflow_filter):
		Stable.append([Med_out[3], np.max(Max_out), Bulk[3], supvels[count], thetime*1000, Mhalo, Zs[count]])
	if ( (theoutflow > Outflow_filter) and (theSFR > SFR_filter)):
		SFRtable.append([Med_out[3], np.max(Max_out), Bulk[3], supvels[count], thetime*1000, Mhalo, theSFR, theoutflow, Zs[count]])
	count = count +1
	

	
Stable = np.array( Stable)
SFRtable = np.array(SFRtable)

thevs = SFRtable[:,3]
flow_weight = SFRtable[:,7]
weighted_v = np.sum(thevs * flow_weight) / np.sum(flow_weight)
print 'weighted halo median velocity ',weighted_v

Vline = haloN+'  '+str(weighted_v) + '  '+percentile+'\n'
g = open(name+'/WV_integrated.txt', 'a')
g.write(Vline)
g.close()


if (doplots):
	#print SFRtable
	np.savetxt('SFRtable_adj'+haloN+'p'+percentile+'.txt', SFRtable)


	thar = ['median 0.25 Rvir', percentile+'th percentile', 'bulk velocity (0.25 Rvir)']
	#plt.plot(time_ar, InnerMed, 'Dk')
	plt.plot(time_ar, med_ar, 'k-')
	plt.plot(time_ar, max_ar, 'b--')
	plt.plot(time_ar, bulk_ar, 'r--')
	plt.plot(time_ar, ar90, 'g--')
	plt.legend(thar, loc = 'best', prop={'size':10})
	plt.xlabel ('cosmic time Myr')
	plt.ylabel ('flow velocity (km/s)')
	#plt.show()
	plt.savefig('velocity_time_Nadj'+haloN+'p'+percentile+'.pdf')
	plt.clf() 

	#plt.plot(time_ar, InnerMed, 'Dk')
	fig2 =plt.plot(time_ar, med_ar, 'k-')
	plt.plot(time_ar, max_ar, 'b--')
	plt.plot(time_ar, ar90, 'g--')
	plt.plot(time_ar, bulk_ar, '.r')
	plt.yscale('log')
	plt.ylim(10, 1e4)
	plt.legend(thar, loc = 'best', prop={'size':10})
	plt.xlabel ('cosmic time Myr ')
	plt.ylabel ('flow velocity (km/s)')
	plt.savefig('velocity_time_log_Nadj'+haloN+'p'+percentile+'.pdf')

	plt.clf()

	plt.yscale('linear')
	figA = plt.plot(SFRtable[:,5], SFRtable[:,3], '.k')
	plt.xscale('log')
	plt.xlabel(r'halo mass (Msun/h)')
	plt.ylabel('90th percentile velocity (km/s)')
	plt.savefig('M_vs_v_Nadj'+haloN+'p'+percentile+'.pdf')

	plt.clf()

	figB = plt.plot(SFRtable[:,6], SFRtable[:,3], '.k')
	plt.xlabel(r'SFR (Msun/yr)')
	plt.ylabel('90th percentile velocity (km/s)')
	plt.xscale('log')
	plt.savefig('SFR_vs_v_Nadj'+haloN+'p'+percentile+'.pdf')
	#plt.show()
	plt.clf()

