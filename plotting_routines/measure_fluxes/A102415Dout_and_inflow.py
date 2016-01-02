import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import date
import sys
import scipy.optimize as optimize

today = date.today()
print today

plot_markers = 'n'
naked = True
jacket = False
specialshell = True #use specialval Extra Shell instead of whatever input was
doublespecial = True
specialval = 2 #for 1.0 Rvir, specialval = 5
Dspecialval = 1  #0 for 10 kpc, 1 for 25 kpc, 2 for 50 kpc - confusing but this is inner region
SFR_correct_factor = 1.15  #to correct for stars formed between snapshots 
doplots = False 

#doublespeciallabel = '10 kpc'
#speciallabel = '25 kpc'
speciallabel = '50 kpc'
doublespeciallabel = '25 kpc'


#speciallabel = '1.0'

#these parameters should be set by command prompt
zstart =999
zend = 999
stagger = 999 #the stagger between start of integration for SFR and for outflow rate.

little_h = 0.7 
paramfac = 1e5
adaptive_paramfac = 'y'

#z_wanted = 


if (len(sys.argv) < 2):
	print 'syntax  blah.py Nhalo filename filename2  [zstart zend stagger name]'
	sys.exit()

haloN = str(sys.argv[1])
filename = str(sys.argv[2])
filename2 = str(sys.argv[3])


#sigh i suck at coding 
uselims = 'n'  #this is only for the y limits on the plots. Defaults set below 
subtract_delay = 'y'  #this is to "go back in time" after the SFR is computed. 
slack_time = 'n' #this one applies to whether the outflow episode can be prolonged by the "Outflow_1st_par" things. 
slack_delay = 'y'   #this is to prolong an SF episode after the fixed "delay" given. THe only meaning of "delay" is how long the SF episode is integrate for
stagger_outflow = 'y' # this is to officially start the integration for outflow total at a time after the SF integration. Given by "stagger"

the_y_lims_zoomed = (-15, 15) 
#default for 0.25 and 0.95 Rvir
ShellUseInner = '2'
ShellUseOuter = '9'

ShellUseInnerorig = ShellUseInner

if (len(sys.argv) >= 6):
	zstart = float(sys.argv[4])
	zend = float(sys.argv[5])
	stagger = float(sys.argv[6])
	name = str(sys.argv[7])
	print 'got zstart, zend, stagger name from command prompt'





if (zstart > 3.5):
	z_wanted = np.array([4.0, 3.5, 3.0, 2.5, 2.0])
	zmedpoint = 3.0
if (zstart > 1 and zstart < 3.5):
	z_wanted = np.array([2.0, 1.5, 1.25, 1.0, 0.75, 0.5])
	zmedpoint = 1.25
if (zstart < 1):
	z_wanted = np.array([0.5, 0.25, 0.1, 0.0])
	zmedpoint = 0.25

zmedpoint = max(zmedpoint, zend)

z_wanted_str = []
for q in z_wanted:
	z_wanted_str.append(str(q))
a_wanted = 1.0 / (1 + z_wanted)

print 'using shells ', ShellUseInner, ShellUseOuter
SFR_smoothing = 20 #Myr
delay = 60 # minimal length of SF interval
interval = 60 #minimal length of outflow interval lets keep interval = delay for now
max_interval = 250

#delay = stagger
#interval = stagger
# i think this thing serves no function


#for initiating


def line_to_string(line, newline = 'y'):
	linestring = ''
	for q in line:
		linestring += "{0:.6g}".format(q) + '  '
	if (newline == 'y'):
		linestring+=' \n'
	return linestring

finname = name+'/tasz.txt'
f = open(finname)
dars = np.loadtxt(f)
theas = dars[:,0]
thets = dars[:,1]
f.close()

finname = filename
f = open(finname)
dars = np.loadtxt(f)
Zs = dars[:,1]
Ns = dars[:,0]
cut1 = Zs >= zend
cut2 = Zs <= zstart
cuts = cut1*cut2


Zs = Zs[cuts] 
Ns = Ns[cuts]
zend = max(zend, np.min(Zs))

Rvirused = dars[:,2][cuts]
SFRInner = dars[:,3][cuts]
SFRouter = dars[:,4][cuts]
#SFR1 = SFRInner + SFRouter
#for now excluding SF outside 0.2 Rvir
SFR1 = SFRInner * SFR_correct_factor
DenseGasMass = dars[:,5][cuts]
TotalMass = dars[:,6][cuts]
Mgas = dars[:,7][cuts]
Mstar = dars[:,8][cuts]

fbar = (Mgas + Mstar)/TotalMass
fgas = Mgas/TotalMass
fstar = Mstar/TotalMass




totalSmass = dars[:,9][cuts]
Flow1 = dars[:,10][cuts]
NewSFR_timescale = dars[:,11][cuts]


newtonG = 6.67384e-8
sSFR = SFR1/ (Mstar/little_h) 

Vcirc = np.sqrt((newtonG * TotalMass *  2e33 / little_h) / (Rvirused/little_h * 3.08e21))
Vcirc = Vcirc / 1e5
#Vcirc *= np.sqrt(1.0 + Zs)


getvcirc = interp1d(Zs[::-1], Vcirc[::-1], kind='linear')
getmass = interp1d(Zs[::-1], TotalMass[::-1], kind='linear')
getRvir  = interp1d(Zs[::-1], Rvirused[::-1], kind='linear')
getMstar = interp1d(Zs[::-1], Mstar[::-1], kind='linear')


print 'vcirc at z=',zend,getvcirc(zend),' km/s'
print 'mass at z=',zend,getmass(zend),' Msun/h'
print 'rvir at z=',zend,getRvir(zend),' kpc/h'
print 'Mstar at z=2',zend,getMstar(zend), ' Msun/h'

mend = float(getmass(zend))
vend = float(getvcirc(zend))
mstarend = float(getMstar(zend))

#mstart = float(getmass(zstart))
#mstarstart = float(getMstar(zstart))
#vstart = float(getvcirc(zstart))

mstart = TotalMass[0]
mstarstart = Mstar[0]
vstart = Vcirc[0]

if (adaptive_paramfac == 'y'):
	#paramfac =  20 *  (1e11/mend)
	#paramfac = 100 * 1e10/mstarend
	paramfac = max(1 * (1e10/mstarend)**2.0, 1.0)

SFR_1stpar = 1.0/paramfac
SFR_2ndpar = 1.0
SFR_3rdpar = 5.0/paramfac

#for slack
SFR_4thpar = 1.0/paramfac
SFR_5thpar = 1.0
SFR_6rhpar = 2.0/paramfac



massstring = "{0:.4g}".format(float(getmass(zend))/little_h)
#massstring = '%.3E' % Decimal(massstring)
#Flow1 = dars[:,7][cuts]
#Flow_Shells = []

#Cross_Shells = [] 


#needs to be one less than the actual number used - 0 counts.
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

count = 0
NumExtraShells = 6
Extra_FlowShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow= dars[:,riteindex][cuts]
	Extra_FlowShells.append(TheFlow)
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

count = 0
Extra_CrossShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheCross= (dars[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	Extra_CrossShells.append(TheCross)
	count += 1



NPartCross_ar = []
NPartFlux_ar = []
NpartShell_ar = []

starter1 = -2
starter2 = -3
starter3 = -4 
count = 0

while (count <= NumShells):
	anticount = NumShells - count - 1
	riteindex = -6*anticount+starter1
	riteindex2 = -6*anticount+starter2
	riteindex3 = -6*anticount+starter3
	NpartInOutflow_Cross = dars[:,riteindex][cuts]
	NpartInOutflow_Flux = dars[:,riteindex2][cuts]
	NpartInShell = dars[:,riteindex3][cuts]
	
	NPartCross_ar.append(NpartInOutflow_Cross)
	NPartFlux_ar.append(NpartInOutflow_Flux)
	NpartShell_ar.append(NpartInShell)
	#print 'The Flow of Shell ', count,' is index ',riteindex
	count += 1

NPartCross_ar=np.array(NPartCross_ar)
NPartFlux_ar=np.array(NPartFlux_ar)
NpartShell_ar=np.array(NpartShell_ar)


partmass = 2.7e4

Outflow_fraction = NPartFlux_ar/NpartShell_ar
Total_Outflow_Flux = (NPartFlux_ar* partmass / little_h) / NewSFR_timescale

TotalOutflow_fraction = sum(NPartFlux_ar[0:10])/sum(NpartShell_ar[0:10])

TheTotalOutflow = sum(Total_Outflow_Flux[1:10])
#TheTotalOutflow = np.copy(Total_Outflow_Flux[1])
#count = 2 
#while (count <= 9):
#	TheTotalOutflow += Total_Outflow_Flux[count]
#	count+=1 


f.close()

	

Flow_in = Flow_Shells[int(ShellUseInner)]
Flow_out = Flow_Shells[int(ShellUseOuter)]

Cross_in = Cross_Shells[int(ShellUseInner)]
Cross_out = Cross_Shells[int(ShellUseOuter)]

if (specialshell):
	Flow_out = Extra_FlowShells[int(specialval)]
	Cross_out = Extra_CrossShells[int(specialval)]
if (doublespecial):
	Flow_in = Extra_FlowShells[int(Dspecialval)]
	Cross_in = Extra_CrossShells[int(Dspecialval)]	

maxy = min((max(np.max(SFR1), np.max(Flow_in))), 100)

finname2 = filename2
f2 = open(finname2)
dars2 = np.loadtxt(f2)

InFlow_Shells = []
InCross_Shells = []


starter = -5 
count = 0
while (count < NumShells):
	anticount = NumShells - count  -1 
	riteindex = -6*anticount+starter
	TheFlow = dars2[:,riteindex][cuts]
	#print 'The InFlow of Shell ', count,' is index ',riteindex
	InFlow_Shells.append(TheFlow)
	count +=1

count = 0
Extra_InFlowShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheFlow= dars2[:,riteindex][cuts]
	Extra_InFlowShells.append(TheFlow)
	count += 1



starter = -1 
count = 0
while (count < NumShells):
	anticount = NumShells - count  -1
	riteindex = -6*anticount+starter
	TheCross = (dars2[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale * -1 
	#print 'The Cross of Shell ', count,' is index ',riteindex
	InCross_Shells.append(TheCross)
	count+=1
	
count = 0
Extra_InCrossShells = []
while (count < (NumShells + NumExtraShells)):
	anticount = NumShells + NumExtraShells - count - 1
	riteindex = -6*anticount+starter
	TheCross= (dars2[:,riteindex][cuts]* 1e10 / little_h) / NewSFR_timescale
	Extra_InCrossShells.append(TheCross)
	count += 1


Flow_in2 = InFlow_Shells[int(ShellUseInner)]
Flow_out2 = InFlow_Shells[int(ShellUseOuter)]

Cross_in2 = InCross_Shells[int(ShellUseInner)] 
Cross_out2 = InCross_Shells[int(ShellUseOuter)]


if (specialshell):
	Flow_out2 = Extra_InFlowShells[int(specialval)]
	Cross_out2 = Extra_InCrossShells[int(specialval)]
	ShellUseOuter = 'S'+str(specialval)
if (doublespecial):
	Flow_in2 = Extra_InFlowShells[int(Dspecialval)]
	Cross_in2 = Extra_InCrossShells[int(Dspecialval)]
	ShellUseInner = 'S'+str(Dspecialval)


count = 0
#print 'showing infall rates'
while (count < len(Flow_in2)):
	#print count, Zs[count], Cross_in2[count], Flow_in2[count] 
	count +=1

#Flow_in2 = dars2[:,-31][cuts]
#Flow_out2 = dars2[:,-3][cuts]
f2.close()


miny = max(-50, np.min(Flow_in2))
the_y_lims = (miny, maxy)


As = 1.0 / (1.0 + Zs)

tfora = interp1d(theas, thets)
afort = interp1d(thets, theas)
Ts = tfora(As)


Ts *= 1000.0

T_wanted = tfora(a_wanted)
T_wanted *= 1000


total_expelled = np.trapz(Flow_in, Ts, 0.1)
total_formed =  np.trapz(SFR1, Ts, 0.1)
total_expelled_cross = np.trapz(Cross_in, Ts, 0.1)
total_gained = np.trapz(Flow_in2, Ts, 0.1)
total_gained_cross = np.trapz(Cross_in2, Ts, 0.1)

totalmasses = []
cumtotalmasses = []

half_total_expelled = total_expelled/2.0

for shell in Flow_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thetotalmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			mass_so_far =np.trapz(shell[0:count+1], Ts[0:count+1], 0.1)
			cumtotalmass.append(mass_so_far)			
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	totalmasses.append(thetotalmass)
	cumtotalmasses.append(cumtotalmass)

count = 0 
magiccount = 0 
magicZ = Zs[count]
for p in cumtotalmasses[int(ShellUseInnerorig)]:
	#print 'doo be doo ',count, Ns[count], p, half_total_expelled
	if (p >= half_total_expelled):
		print 'half outflow point reached!' 
		print count, Ns[count], Zs[count]
		magiccount = count
		magicZ = Zs[count]
		break
	count +=1

incrosstotalmasses = []
incrosscumtotalmasses = []
for shell in InCross_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thetotalmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			cumtotalmass.append(np.trapz(shell[0:count+1], Ts[0:count+1], 0.1))
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	incrosstotalmasses.append(thetotalmass)
	incrosscumtotalmasses.append(cumtotalmass)

intotalmasses = []
incumtotalmasses = []
for shell in InFlow_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thetotalmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			cumtotalmass.append(np.trapz(shell[0:count+1], Ts[0:count+1], 0.1))
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	intotalmasses.append(thetotalmass)
	incumtotalmasses.append(cumtotalmass)



#for shell in Cross_Shells:
#	if (len(shell) > 1): 
#		thecrossmass = np.trapz(shell, Ts, 0.1)
#	else: thecrossmass = 0
#	crossmasses.append(thecrossmass)

crossmasses = [] 
crosscumtotalmasses = []
for shell in Cross_Shells:
	cumtotalmass = []
	if (len(shell) > 1): 
		thecrossmass = np.trapz(shell, Ts, 0.1)
		count = 0 
		while (count < len(shell)):
			cumtotalmass.append(np.trapz(shell[0:count+1], Ts[0:count+1], 0.1))
			count+=1
	else: thetotalmass = 0
	cumtotalmass = np.array(cumtotalmass)
	crossmasses.append(thecrossmass)
	crosscumtotalmasses.append(cumtotalmass)


totalmasses = np.array(totalmasses)
shellarray = np.arange(NumShells)
crossmasses = np.array(crossmasses)
intotalmasses = np.array(intotalmasses)
incrosstotalmasses = np.array(incrosstotalmasses)
masscolumns = np.column_stack([shellarray, totalmasses, crossmasses, intotalmasses, incrosstotalmasses])


foutname = (name+'/' + str(today) + 'total_mass_per_shell'+ShellUseInner+ ShellUseOuter +'N'+haloN+'.txt')
np.savetxt(foutname, masscolumns, fmt='%1.6g')

ML_from_total = total_expelled/total_formed
ML_from_total_cross = total_expelled_cross/total_formed

SFRs = interp1d(Ts, SFR1, kind='linear')
SFRsOuter = interp1d(Ts, SFRouter, kind='linear')
Outflows_Inner = interp1d(Ts, Flow_in,  kind='linear')
Outflows_Outer = interp1d(Ts, Flow_out, kind='linear')
Cross_Inner = interp1d(Ts, Cross_in, kind='linear')
Cross_Outer = interp1d(Ts, Cross_out, kind='linear')


Flow_Shells_Int = [] 

for shell in Flow_Shells:
	fs = interp1d(Ts, shell, kind='linear')
	Flow_Shells_Int.append(fs)

Cross_Shells_Int = []
for shell in Cross_Shells:
	fs = interp1d(Ts, shell, kind='linear')
	Cross_Shells_Int.append(fs)

#New_Outflow = interp1d(Ts, NewOutflowRate, kind='linear')


Inflows_Inner = interp1d(Ts, Flow_in2,  kind='linear')
Inflows_Outer = interp1d(Ts, Flow_out2, kind='linear')

count = 0
Tbegin = Ts[0]
Tend = Ts[len(Ts)-1]

Mformed  = 0
Mexpelled = 0
Mcrossed = 0
Mnew = 0
AvgSFR =0 
AvgSFR_for_delay = 0
AvgOutflowInner = 0
AvgOutflowOuter = 0
AvgInflowInner = 0
AvgInflowOuter = 0
AvgNewOutflows = 0

AvgSFR_ar = []
AvgSFR_d_ar = []
AvgOutflowInner_ar = []
AvgOutflowOuter_ar = []
AvgInflowInner_ar = []
AvgInflowOuter_ar = []
AvgNewOutflows_ar = [] 
Mformed_ar = []
Mexpelled_ar = []
Mcrossed_ar = []
Mnew_ar = []

Mvir_ar = []
vcirc_ar = []
Mstar_ar = []
z_ar = []
t_ar = []

Flow_total_interval_ar = []
interval_start_ar =[]
flow_interval_start_ar = []
delay_ar = []
interval_end_ar =[]
zero_ar = []
tstep = 1

interval_count = 0

if (interval <= SFR_smoothing):
	print 'interval must be larger than SFR_smoothing!'
	sys.exit()

the_time = Tbegin
print Tbegin, Tend
SFRprev = 0



while (the_time < (Tend-max(delay,(interval+stagger)))):
	if (((SFRs(the_time) > SFR_1stpar) and (SFRs(the_time) > SFRprev*SFR_2ndpar)) or SFRs(the_time )> SFR_3rdpar):
		Flow_total_interval = np.zeros(10)
		interval_count += 1
		print '###############################################'
		print 'beginning interval ', interval_count, the_time
		interval_start_ar.append(the_time)
		delay_end_time = the_time + delay
		delay_actual = 0
		if (delay_end_time < Tend):
			while (the_time <= delay_end_time):
				#AvgSFR += SFRs(the_time)*tstep
				AvgSFR_for_delay += SFRs(the_time)*tstep
				Mformed += SFRs(the_time)*tstep
				the_time+=tstep
				delay_actual += tstep
				#Mformed += SFRs(the_time)*tstep 
			the_delay_slack = 0
			if (slack_delay == 'y'):
				while (the_time < Tend and ((SFRs(the_time) > SFR_4thpar and SFRs(the_time) > SFR_5thpar * SFRprev) or SFRs(the_time)>SFRs(the_time-tstep) or SFRs(the_time) > SFR_6rhpar) and (the_delay_slack < max_interval)):
					#AvgSFR += SFRs(the_time)*tstep
					AvgSFR_for_delay += SFRs(the_time)*tstep
					Mformed += SFRs(the_time)*tstep
					the_time+=tstep
					the_delay_slack += tstep
					delay_actual += tstep
					#Mformed += SFRs(the_time)*tstep 
				print 'used delay slack ',the_delay_slack
				print 'delay went until ',the_time
			delay_ar.append(the_time)
		else:
			print 'skipping delay'
		if (subtract_delay == 'y'):
			the_time -= delay_actual 
			the_time += 0.00001
			print 'subtracted interval, time is now ',the_time
		if (stagger_outflow  == 'y'):
			print 'staggering outflow'
			the_time += stagger			
		print 'beginning interval'
		interval_begin_time = the_time
		flow_interval_start_ar.append(the_time)
		actual_interval = max(interval, delay_actual)
		interval_end_time = the_time+actual_interval

		if (interval_end_time > Tend):
			interval_end_time = Tend
			actual_interval = Tend - the_time
		while (the_time < interval_end_time):
			AvgSFR += SFRs(the_time)*tstep
			AvgOutflowInner += Outflows_Inner(the_time)*tstep
			AvgOutflowOuter += Outflows_Outer(the_time)*tstep
			AvgInflowInner += Inflows_Inner(the_time)*tstep
			AvgInflowOuter += Inflows_Outer(the_time)*tstep 
			#AvgNewOutflows += New_Outflow(the_time)*tstep
			#Mformed += SFRs(the_time)*tstep   
			#only do it this way if you aren't using the delay option
			Mexpelled += Outflows_Inner(the_time)*tstep
			Mcrossed += Cross_Inner(the_time)*tstep
			#Mnew +=  New_Outflow(the_time)*tstep
			
			pocount = 0 
			for po in Flow_total_interval:
				Flow_total_interval[pocount] += Flow_Shells_Int[pocount](the_time)*tstep
				pocount+=1 
				#print 'po ', po
			the_time += tstep
		the_slack = 0		
		interval_end_ar.append(interval_end_time)
		zero_ar.append(0.0)
		AvgSFR /= (actual_interval+delay_actual + the_slack)
		AvgOutflowInner /= (actual_interval+the_slack)
		AvgInflowInner /= (actual_interval+the_slack)
		AvgOutflowOuter /= (actual_interval+the_slack)
		AvgInflowOuter /= (actual_interval+the_slack)
		AvgNewOutflows /= (actual_interval+the_slack)
		if (delay > 0):
			AvgSFR_for_delay /= (delay_actual)
			if (AvgSFR > AvgSFR_for_delay):
				print 'the delay must be wrong coz the average SFR is bigger'
				AvgSFR_for_delay = AvgSFR		
		print 'for interval t=',(interval_begin_time),' to ',the_time
		print 'Avg SFR ', AvgSFR
		print 'Avg SFR during delay', AvgSFR_for_delay
		print 'Avg Outflow Inner ',AvgOutflowInner
		print 'Avg Outflow Outer ',AvgOutflowOuter
		print 'Avg Inflow Inner ',AvgInflowInner
		print 'Avg Inflow Outer ',AvgInflowOuter
		print 'M formed ',Mformed
		print 'M expelled ',Mexpelled
		print 'M crossed ',Mcrossed
		#print 'M new ',Mnew

		AvgSFR_ar.append(AvgSFR)
		AvgSFR_d_ar.append(AvgSFR_for_delay)
		AvgOutflowInner_ar.append(AvgOutflowInner)
		AvgOutflowOuter_ar.append(AvgOutflowOuter)
		AvgInflowInner_ar.append(AvgInflowInner)
		AvgInflowOuter_ar.append(AvgInflowOuter)
		
		
		
		Mformed_ar.append(Mformed)
		Mexpelled_ar.append(Mexpelled)
		Mcrossed_ar.append(Mcrossed)
		#Mnew_ar.append(Mnew)
		thea = afort(the_time/1000.0 - (the_time - interval_begin_time)/2000 )
		thez = 1.0 / thea  - 1.0

		Mvir_ar.append(getmass(thez))
		vcirc_ar.append(getvcirc(thez))
		Mstar_ar.append(getMstar(thez))
		z_ar.append(thez)
		#t_ar.append(the_time/1000.0 - (the_time - interval_begin_time)/2000)
		t_ar.append(the_time)
		print 'fti ',Flow_total_interval
		Flow_total_interval_ar.append(Flow_total_interval)

		
		AvgSFR =0 
		AvgSFR_for_delay = 0
		AvgOutflowInner = 0
		AvgOutflowOuter = 0
		AvgInflowInner = 0
		AvgOutflowInner = 0
		Mformed = 0
		Mexpelled = 0
		Mcrossed = 0
		#Mnew = 0
		del(Flow_total_interval)
		if (stagger_outflow  == 'y'):
			print 'subtracting out stagger'
			the_time -= stagger			
	if ((the_time - Tbegin) > SFR_smoothing):
		SFRprev = SFRs(the_time-SFR_smoothing)
	the_time += tstep

print 'done, doing plots'
print 'total expelled via trapz ',total_expelled
tot_exp_ep = np.sum(np.array(Mexpelled_ar))
print 'total expelled found in episodes', tot_exp_ep
print 'total crossed via trapz', total_expelled_cross
tot_cross_ep =  np.sum(np.array(Mcrossed_ar))
print 'total crossed in episodes ', tot_cross_ep
print 'total formed via trapz ',total_formed
tot_formed_ep = np.sum(np.array(Mformed_ar))
print 'total formed in episodes ',tot_formed_ep

ep_rat_flux = tot_exp_ep/total_expelled
ep_rat_cross = tot_cross_ep/total_expelled_cross
ep_rat_form = tot_formed_ep/total_formed

thecolors = ['r', 'g', 'k','#a52a2a', '#ffd700', '#7fff00', '#006400', '#00ced1', '#483d8b', '#8b7765']
thecolors.extend(thecolors)
thecolors.extend(thecolors)
thecolors.extend(thecolors)

Shells = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]

innerstring = '0.'+ShellUseInner+r'5'
outerstring = '0.'+ShellUseOuter+r'5'
if (specialshell):
	outerstring = speciallabel
if (doublespecial):
	innerstring = doublespeciallabel
thar_10 =[r'SFR$\times$10', r'outer SFR$\times$10', innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow', r'M_{gas} with n > 1 ${\rm cm}^{-3}$ ($10^7 M_{\odot}$)']

thar22 = ['SFR',  innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow', innerstring+r' $R_{\rm vir}$ outflow traced',outerstring+r' $R_{\rm vir}$ outflow traced', innerstring+r' $R_{\rm vir}$ inflow traced',  outerstring+r' $R_{\rm vir}$ inflow traced']
thar22_10 = [r'SFR$\times$10',  innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow', innerstring+r' $R_{\rm vir}$ outflow traced',outerstring+r' $R_{\rm vir}$ outflow traced', innerstring+r' $R_{\rm vir}$ inflow traced',  outerstring+r' $R_{\rm vir}$ inflow traced']



#thar22b =['SFR', '0.2 $R_{\rm vir}$ outflow', '0.9 $R_{\rm vir}$ outflow', '0.2 $R_{\rm vir}$ inflow', '0.9 $R_{\rm vir}$ inflow', '0.2 $R_{\rm vir}$ crossed']

thar22b =['SFR', innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow']

thar22b_10 =[r'SFR$\times$10', innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow']



thar22c  = ['SFR',  innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow', innerstring+r' $R_{\rm vir}$ outflow traced',outerstring+r' $R_{\rm vir}$ outflow traced']
thar22c_10 =  [r'SFR$\times$10',  innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow', innerstring+r' $R_{\rm vir}$ outflow traced',outerstring+r' $R_{\rm vir}$ outflow traced']

thar2 =[innerstring+r' $R_{\rm vir}$ outflow', outerstring+r' $R_{\rm vir}$ outflow', innerstring+r' $R_{\rm vir}$ inflow', outerstring+r' $R_{\rm vir}$ inflow']

if (doplots):	###################################################################
	#Figure 1
	#outandinflow.pdf 
	###################################################################
	fig1 = plt.figure(figsize=(10,8))

	plt.plot(Ts, SFR1*10 , color='k', markevery=1, lw=2)
	plt.plot(Ts, SFRouter*10, color='g',markevery=1, lw=1)
	plt.plot(Ts, Flow_in , color='r', ls='--', markevery=1, lw=2)
	plt.plot(Ts, Flow_out , color='r', ls=':',  markevery=1, lw=2)
	plt.plot(Ts, Flow_in2 , color='b', ls='--', markevery=1, lw=2)
	plt.plot(Ts, Flow_out2 , color='b', ls=':',  markevery=1, lw =2)
	plt.plot(Ts, (DenseGasMass/1e7), color='#483d8b', ls='-.', lw=2)
	if (plot_markers == 'y'):
		plt.plot(interval_start_ar, zero_ar, color='b', marker='x', markevery=1)
		plt.plot(flow_interval_start_ar, zero_ar, color ='m', marker='x', markevery=1)
		plt.plot(interval_end_ar, zero_ar, color='r', marker='v', markevery=1)
		plt.plot(delay_ar, zero_ar, color='g', marker='^', markevery=1)


	if (uselims == 'y'):
		plt.ylim(-15, 15)
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	plt.ylabel(r'SFR$\times$10, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)', fontsize=20)
	plt.xlabel('cosmic time (Myr)',fontsize=20)
	ax = plt.gca()

	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	thextitle = x_range[0] + 0.05*(x_range[1]-x_range[0])
	theytitle = y_range[0] + 0.05*(y_range[1]-y_range[0])
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'

	#plt.text(thextitle,theytitle, figure_title)
	#plt.legend(thar_10, loc = 'best', prop={'size':12})
	plt.legend(thar_10, loc = 'best', prop={'size':12}, title=figure_title)

	#plt.title(figure_title, y=1.08, x=0.83, fontsize=15)
	#plt.tight_layout()

	ax2 = ax.twiny()




	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)

	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(T_wanted)
	ax2.set_xticklabels(z_wanted_str)
	ax2.set_xlabel('z', fontsize=20)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax2.xaxis.set_tick_params(width=2.5, length=8)
	ax2.yaxis.set_tick_params(width=2.5, length=8)    



	#xrange = ax2.get_xlim()
	#titlex = xrange[0]
	#print xrange

	#plt.title(figure_title, y=1.08, x=0.83, fontsize=15)
	#plt.tight_layout()


	plt.savefig(name+'/' + str(today) + 'outandinflow' +ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')

	###################################################################
	#Figure 2
	#outandinflow_skinny.pdf
	###################################################################

	#skinny_axes = [min(Ts), max(Ts), the_y_lims[0], the_y_lims[1]]
	#figure 2 skinny version of cosmic time
	fig22 = plt.figure(figsize=(10,8))
	plt.plot(Ts, SFR1*10 , color='k', markevery=1, lw=3)
	plt.plot(Ts, Flow_in , color='r', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out , color='r', ls=':',  markevery=1, lw=3)
	plt.plot(Ts, Flow_in2 , color='b', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out2 , color='b', ls=':',  markevery=1, lw=3)
	if (not naked):
		plt.plot(Ts, Cross_in , color='#483d8b', ls='--', markevery=1, lw=1.5)
		plt.plot(Ts, Cross_out , color='#483d8b', ls=':',markevery=1, lw=1.5)
		if (jacket):
			plt.plot(Ts, Cross_in2, color='#8B0000', ls ='--', lw=1.5)
			plt.plot(Ts, Cross_out2, color='#8B0000', ls =':',lw=1.5)

	#plt.plot(Ts, NewOutflowRate, color = 'g',ls=':', lw=3)
	#plt.axhline(linewidth=4, color="g")        # inc. width of x-axis and color it green
	#plt.axvline(linewidth=4, color="r")        # inc. width of y-axis and color it red

	#plt.plot(Ts, (DenseGasMass/1e7), color='#483d8b', ls='-.')
	if (plot_markers == 'y'):
		plt.plot(interval_start_ar, zero_ar, 'xb')
		plt.plot(flow_interval_start_ar, zero_ar, color ='m', marker='x', markevery=1)
		plt.plot(interval_end_ar, zero_ar, 'vr')
		plt.plot(delay_ar, zero_ar, '^g')
	if (uselims == 'y'): #disabling for now
		plt.ylim(the_y_lims)
		plt.xlim((1200, 3500))
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	#plt.ylabel('SFR, Outflow, Inflow  $M_{\odot}$${\rm yr}^{-1}$',fontsize=20, weight='semibold')
	plt.ylabel(r'SFR$\times$10, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
	plt.xlabel('cosmic time (Myr)',fontsize=20)
	if (naked):
		plt.legend(thar22b_10, loc = 'best',  prop={'size':12}, title=figure_title)
	else:
		if (jacket):
			plt.legend(thar22_10, loc = 'best',  prop={'size':12}, title=figure_title)
		else:
			plt.legend(thar22c_10, loc = 'best',  prop={'size':12}, title=figure_title)
	ax = plt.gca()
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	thextitle = x_range[0] + 0.05*(x_range[1]-x_range[0])
	theytitle = y_range[0] + 0.05*(y_range[1]-y_range[0])
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'

	#plt.text(thextitle,theytitle, figure_title)

	ax2 = ax.twiny()
	#print ax.xrange
	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)
	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(T_wanted)
	ax2.set_xticklabels(z_wanted_str)
	ax2.set_xlabel('z', fontsize=20)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax2.xaxis.set_tick_params(width=2.5, length=8)
	ax2.yaxis.set_tick_params(width=2.5, length=8)    
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'


	#print plt.get_xlim()
	plt.draw()
	plt.savefig(name+'/' + str(today) + 'outandinflow_skinny' +ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')

	#garbage code but keeping here in case i need it sometime
	#xrange = ax.get_xlim()
	#xtitle = (xrange[0] + xrange[1])/2
	#yrange = ax.get_ylim()
	#ytitle = yrange[1]*1.15
	#plt.tight_layout(pad=1.3, h_pad = 2.0)
	#plt.text(xtitle, ytitle, figure_title,       horizontalalignment='center',  fontsize=20)
	#plt.title(figure_title, y=1.08, fontsize = 15, loc='left')

	###################################################################
	#Figure 3
	#outandinflow_skinnyB.pdf
	###################################################################


	fig22b =  plt.figure(figsize=(10,8))
	plt.plot(Ts, SFR1 , color='k', markevery=1, lw=3)
	#plt.plot(Ts, SFRouter, color='g',markevery=1)
	plt.plot(Ts, Flow_in , color='r', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out , color='r', ls=':',  markevery=1, lw=3)
	plt.plot(Ts, Flow_in2 , color='b', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out2 , color='b', ls=':',  markevery=1, lw=3)
	if (not naked):
		plt.plot(Ts, Cross_in , color='#483d8b', ls='--', markevery=1, lw=1.5)
		plt.plot(Ts, Cross_out , color='#483d8b', ls=':',markevery=1, lw=1.5)
		if (jacket):
			plt.plot(Ts, Cross_in2, color='#8B0000', ls ='--', lw=1.5)
			plt.plot(Ts, Cross_out2, color='#8B0000', ls =':',lw=1.5)

	if (plot_markers == 'y'):
		plt.plot(interval_start_ar, zero_ar, 'xb')
		plt.plot(flow_interval_start_ar, zero_ar, color ='m', marker='x', markevery=1)
		plt.plot(interval_end_ar, zero_ar, 'vr')
		plt.plot(delay_ar, zero_ar, '^g')
	plt.ylim(the_y_lims)
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	#plt.ylabel('SFR, Outflow, Inflow  $M_{\odot}$${\rm yr}^{-1}$',fontsize=20, weight='semibold')
	plt.ylabel(r'SFR, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
	plt.xlabel('cosmic time (Myr)',fontsize=20)
	#plt.title(r''+str(zstart)+' > z > '+str(zend)+' $M_h$ ='+massstring+'$M_\odot / h$ (at z='+str(zend)+')', fontsize=15)
	if (naked):
		plt.legend(thar22b, loc = 'best',  prop={'size':12}, title=figure_title)
	else:
		if (jacket):
			plt.legend(thar22, loc = 'best',  prop={'size':12}, title=figure_title)
		else:
			plt.legend(thar22c, loc = 'best',  prop={'size':12}, title=figure_title)
	ax = plt.gca()
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	thextitle = x_range[0] + 0.05*(x_range[1]-x_range[0])
	theytitle = y_range[0] + 0.05*(y_range[1]-y_range[0])
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	ax2 = ax.twiny()

	#plt.text(thextitle,theytitle, figure_title)


	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)

	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(T_wanted)
	ax2.set_xticklabels(z_wanted_str)
	ax2.set_xlabel('z', fontsize=20)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax2.xaxis.set_tick_params(width=2.5, length=8)
	ax2.yaxis.set_tick_params(width=2.5, length=8)    
	#figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	#plt.title(figure_title, y=1.08, x=0.83, fontsize=15)
	#plt.tight_layout()
	plt.draw()
	plt.savefig(name+'/' + str(today) + 'outandinflow_skinnyB_' +ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')
	plt.clf()

	###################################################################
	#Figure 3c
	#outandinflow_skinnyC.pdf
	###################################################################


	fig22c =  plt.figure(figsize=(10,8))
	plt.plot(Ts, SFR1 , color='k', markevery=1, lw=3)
	plt.plot(Ts, Flow_in , color='r', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out , color='r', ls=':',  markevery=1, lw=3)
	plt.plot(Ts, Flow_in2 , color='b', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out2 , color='b', ls=':',  markevery=1, lw=3)
	if (not naked):
		plt.plot(Ts, Cross_in , color='#483d8b', ls='--', markevery=1, lw=1.5)
		plt.plot(Ts, Cross_out , color='#483d8b', ls=':',markevery=1, lw=1.5)
		if (jacket):
			plt.plot(Ts, Cross_in2, color='#8B0000', ls ='--', lw=1.5)
			plt.plot(Ts, Cross_out2, color='#8B0000', ls =':',lw=1.5)


	if (plot_markers == 'y'):
		plt.plot(interval_start_ar, zero_ar, 'xb')
		plt.plot(flow_interval_start_ar, zero_ar, color ='m', marker='x', markevery=1)
		plt.plot(interval_end_ar, zero_ar, 'vr')
		plt.plot(delay_ar, zero_ar, '^g')
	plt.ylim([-60, 120])
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	#plt.ylabel('SFR, Outflow, Inflow  $M_{\odot}$${\rm yr}^{-1}$',fontsize=20, weight='semibold')
	plt.ylabel(r'SFR, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
	plt.xlabel('cosmic time (Myr)',fontsize=20)
	#plt.title(r''+str(zstart)+' > z > '+str(zend)+' $M_h$ ='+massstring+'$M_\odot / h$ (at z='+str(zend)+')', fontsize=15)
	if (naked):
		plt.legend(thar22b, loc = 'best',  prop={'size':12}, title=figure_title)
	else:
		if (jacket):
			plt.legend(thar22, loc = 'best',  prop={'size':12}, title=figure_title)
		else:
			plt.legend(thar22c, loc = 'best',  prop={'size':12}, title=figure_title)
	ax = plt.gca()
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	thextitle = x_range[0] + 0.05*(x_range[1]-x_range[0])
	theytitle = y_range[0] + 0.05*(y_range[1]-y_range[0])
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	ax2 = ax.twiny()

	#plt.text(thextitle,theytitle, figure_title)


	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)

	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(T_wanted)
	ax2.set_xticklabels(z_wanted_str)
	ax2.set_xlabel('z', fontsize=20)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax2.xaxis.set_tick_params(width=2.5, length=8)
	ax2.yaxis.set_tick_params(width=2.5, length=8)    
	#figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	#plt.title(figure_title, y=1.08, x=0.83, fontsize=15)
	#plt.tight_layout()
	plt.draw()
	plt.savefig(name+'/' + str(today) + 'outandinflow_skinnyC_' +ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')
	plt.clf()

	hugecuts = Cross_in > 10 
	count = 0
	print 'the redshifts of the biggest outflows:'
	while (count < len(Zs[hugecuts])):
		print Zs[hugecuts][count], Cross_in[hugecuts][count], Flow_in[hugecuts][count]
		count+=1

	print 'the redshifts of the glitchiest inflows: '
	hugecuts = ( Cross_in2 < 2*  Flow_in2) * (Cross_in2 < -5)
	count = 0
	while (count < len(Zs[hugecuts])):
		print Zs[hugecuts][count], Cross_in2[hugecuts][count], Flow_in2[hugecuts][count]
		count+=1

	###################################################################
	#Figure 3c
	#outandinflow_skinnyD.pdf
	###################################################################


	fig22d =  plt.figure(figsize=(10,8))
	plt.plot(Ts, SFR1*10 , color='k', markevery=1, lw=3)
	plt.plot(Ts, Flow_in , color='r', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out , color='r', ls=':',  markevery=1, lw=3)
	plt.plot(Ts, Flow_in2 , color='b', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out2 , color='b', ls=':',  markevery=1, lw=3)
	if (not naked):
		plt.plot(Ts, Cross_in , color='#483d8b', ls='--', markevery=1, lw=1.5)
		plt.plot(Ts, Cross_out , color='#483d8b', ls=':',markevery=1, lw=1.5)
		if (jacket):
			plt.plot(Ts, Cross_in2, color='#8B0000', ls ='--', lw=1.5)
			plt.plot(Ts, Cross_out2, color='#8B0000', ls =':',lw=1.5)


	if (plot_markers == 'y'):
		plt.plot(interval_start_ar, zero_ar, 'xb')
		plt.plot(flow_interval_start_ar, zero_ar, color ='m', marker='x', markevery=1)
		plt.plot(interval_end_ar, zero_ar, 'vr')
		plt.plot(delay_ar, zero_ar, '^g')
	plt.ylim([-60, 120])
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	#plt.ylabel('SFR, Outflow, Inflow  $M_{\odot}$${\rm yr}^{-1}$',fontsize=20, weight='semibold')
	plt.ylabel(r'SFR$\times$10, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
	plt.xlabel('cosmic time (Myr)',fontsize=20)
	#plt.title(r''+str(zstart)+' > z > '+str(zend)+' $M_h$ ='+massstring+'$M_\odot / h$ (at z='+str(zend)+')', fontsize=15)
	if (naked):
		plt.legend(thar22b_10, loc = 'best',  prop={'size':12}, title=figure_title)
	else:
		if (jacket):
			plt.legend(thar22_10, loc = 'best',  prop={'size':12}, title=figure_title)
		else:
			plt.legend(thar22c_10, loc = 'best',  prop={'size':12}, title=figure_title)
	ax = plt.gca()
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	thextitle = x_range[0] + 0.05*(x_range[1]-x_range[0])
	theytitle = y_range[0] + 0.05*(y_range[1]-y_range[0])
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	ax2 = ax.twiny()

	#plt.text(thextitle,theytitle, figure_title)


	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)

	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(T_wanted)
	ax2.set_xticklabels(z_wanted_str)
	ax2.set_xlabel('z', fontsize=20)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax2.xaxis.set_tick_params(width=2.5, length=8)
	ax2.yaxis.set_tick_params(width=2.5, length=8)    
	#figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	#plt.title(figure_title, y=1.08, x=0.83, fontsize=15)
	#plt.tight_layout()
	plt.draw()
	plt.savefig(name+'/' + str(today) + 'outandinflow_skinnyD_' +ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')
	plt.clf()


	###################################################################
	#Figure 3c
	#outandinflow_skinnyE.pdf
	###################################################################


	fig22e =  plt.figure(figsize=(10,8))
	plt.plot(Ts, SFR1 , color='k', markevery=1, lw=3)
	plt.plot(Ts, Flow_in , color='r', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out , color='r', ls=':',  markevery=1, lw=3)
	plt.plot(Ts, Flow_in2 , color='b', ls='--', markevery=1, lw=3)
	plt.plot(Ts, Flow_out2 , color='b', ls=':',  markevery=1, lw=3)
	if (not naked):
		plt.plot(Ts, Cross_in , color='#483d8b', ls='--', markevery=1, lw=1.5)
		plt.plot(Ts, Cross_out , color='#483d8b', ls=':',markevery=1, lw=1.5)
		if (jacket):
			plt.plot(Ts, Cross_in2, color='#8B0000', ls ='--', lw=1.5)
			plt.plot(Ts, Cross_out2, color='#8B0000', ls =':',lw=1.5)


	if (plot_markers == 'y'):
		plt.plot(interval_start_ar, zero_ar, 'xb')
		plt.plot(flow_interval_start_ar, zero_ar, color ='m', marker='x', markevery=1)
		plt.plot(interval_end_ar, zero_ar, 'vr')
		plt.plot(delay_ar, zero_ar, '^g')
	#plt.ylim([-60, 120])
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	#plt.ylabel('SFR, Outflow, Inflow  $M_{\odot}$${\rm yr}^{-1}$',fontsize=20, weight='semibold')
	plt.ylabel(r'SFR, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
	plt.xlabel('cosmic time (Myr)',fontsize=20)
	#plt.title(r''+str(zstart)+' > z > '+str(zend)+' $M_h$ ='+massstring+'$M_\odot / h$ (at z='+str(zend)+')', fontsize=15)
	if (naked):
		plt.legend(thar22b, loc = 'best',  prop={'size':12}, title=figure_title)
	else:
		if (jacket):
			plt.legend(thar22, loc = 'best',  prop={'size':12}, title=figure_title)
		else:
			plt.legend(thar22c, loc = 'best',  prop={'size':12}, title=figure_title)
	ax = plt.gca()
	x_range = ax.get_xlim()
	y_range = ax.get_ylim()

	thextitle = x_range[0] + 0.05*(x_range[1]-x_range[0])
	theytitle = y_range[0] + 0.05*(y_range[1]-y_range[0])
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	ax2 = ax.twiny()

	#plt.text(thextitle,theytitle, figure_title)


	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)

	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(T_wanted)
	ax2.set_xticklabels(z_wanted_str)
	ax2.set_xlabel('z', fontsize=20)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ax2.xaxis.set_tick_params(width=2.5, length=8)
	ax2.yaxis.set_tick_params(width=2.5, length=8)    
	figure_title =r' $M_h$ ='+massstring+'$M_\odot$ (at z='+"{0:.4g}".format(zend)+')'
	#plt.title(figure_title, y=1.08, x=0.83, fontsize=15)
	#plt.tight_layout()
	#plt.draw()
	plt.savefig(name+'/' + str(today) + 'outandinflow_skinnyE_' +ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')
	plt.clf()




	###################################################################
	#Figure 4
	#z_outandinflow.pdf
	###################################################################

	fig3 = plt.figure(figsize=(10,6))
	#plt.axis([min(Zs), max(Zs), -15, 15])
	#plt.axis([min(Zs)*0.95, max(Zs)*1.05, -2, 2])
	plt.plot(Zs, SFR1 , color='k', markevery=1, lw=3)
	plt.plot(Zs, SFRouter, color='g',markevery=1, lw=3)
	plt.plot(Zs, Flow_in , color='r', ls='--', markevery=1, lw=3)
	plt.plot(Zs, Flow_out , color='r', ls=':',  markevery=1, lw=3)
	plt.plot(Zs, Flow_in2 , color='b', ls='--', markevery=1, lw=3)
	plt.plot(Zs, Flow_out2 , color='b', ls=':',  markevery=1, lw=3)
	plt.plot(Zs, (DenseGasMass/1e7), color='#483d8b', ls='-.', lw=3)
	plt.gca().invert_xaxis()
	#if (uselims == 'y'):
	#plt.ylim(the_y_lims)
	plt.xlim([min(Zs)*0.95, max(Zs)*1.05])
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	plt.ylabel(r'SFR, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)',fontsize=20)
	plt.xlabel('z', fontsize=20)
	plt.legend(thar22b, loc = 'best', prop={'size':9})
	ax = plt.gca()
	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
		label.set_family('serif')
		#label.set_weight('bold')
	ax.xaxis.set_tick_params(width=2.5, length=8)
	ax.yaxis.set_tick_params(width=2.5, length=8)
	plt.draw()
	plt.savefig(name+'/' + str(today) + 'z_outandinflow'  +ShellUseInner+ ShellUseOuter+'N'+haloN+'.pdf')
	#plt.show()

	###################################################################
	#Figure 5
	#averages.pdf
	###################################################################

	fig2 = plt.figure(figsize=(10, 6))
	#plt.axis([0, max(max(AvgSFR_ar), 5), min(min(AvgInflowInner_ar), min(AvgInflowOuter_ar),-5), max(AvgOutflowInner_ar,AvgOutflowOuter_ar 5)])
	plt.xlabel('Avg SFR (Msun / yr)')
	plt.ylabel('Avg Inflow/Outflow rates (Msun / yr) ')
	plt.plot(AvgSFR_ar, AvgOutflowInner_ar ,'.r', markevery=1)
	#plt.plot(AvgSFR_ar, AvgOutflowOuter_ar ,'.g', markevery=1)
	#plt.plot(AvgSFR_ar, AvgInflowInner_ar ,'.b', markevery=1)
	#plt.plot(AvgSFR_ar, AvgInflowOuter_ar ,'.k', markevery=1)
	plt.legend(thar2, loc = 'best')
	plt.savefig(name+'/' + str(today) + 'averages'  +ShellUseInner+ ShellUseOuter+'N'+haloN+'.pdf')


	###################################################################
	#Figure 6
	#averages_delay.pdf
	###################################################################

	fig4 = plt.figure(figsize=(10, 6))
	#plt.axis([0, max(max(AvgSFR_ar), 5), min(min(AvgInflowInner_ar), min(AvgInflowOuter_ar),-5), max(AvgOutflowInner_ar,AvgOutflowOuter_ar 5)])
	plt.xlabel('Avg SFR during dleay (Msun / yr)')
	plt.ylabel('Avg Inflow/Outflow rates (Msun / yr) ')
	plt.plot(AvgSFR_d_ar, AvgOutflowInner_ar ,'.r', markevery=1)
	#plt.plot(AvgSFR_d_ar, AvgOutflowOuter_ar ,'.g', markevery=1)
	#plt.plot(AvgSFR_d_ar, AvgInflowInner_ar ,'.b', markevery=1)
	#plt.plot(AvgSFR_d_ar, AvgInflowOuter_ar ,'.k', markevery=1)
	plt.legend(thar2, loc = 'best')
	plt.savefig(name+'/' + str(today) + 'averages_delay'  +ShellUseInner+ ShellUseOuter+'N'+haloN+'.pdf')
	#plt.show()

	###################################################################
	#Figure 7
	#formed_expelled.pdf
	#moving some cutting to this section
	###################################################################

	fig5 = plt.figure(figsize=(10, 6))

	Mformed_ar = np.array(Mformed_ar)
	Mexpelled_ar = np.array(Mexpelled_ar)
	Mcrossed_ar = np.array(Mcrossed_ar)

	logMformed_ar = np.log10(np.array(Mformed_ar)+ 1e-10)
	logMexpelled_ar = np.log10(np.array(Mexpelled_ar)+ 1e-10)
	logMcrossed_ar = np.log10(np.array(Mcrossed_ar) + 1e-10)
	bigdog1 = max(Mexpelled_ar)
	bigdog2 = max(Mcrossed_ar)
	cut5 = logMformed_ar > -8
	cut6 = logMexpelled_ar > -8
	cut7 = logMcrossed_ar > -8
	bigdogcut_f = Mexpelled_ar > 0.01 * bigdog1
	bigdogcut_c = Mcrossed_ar > 0.01 * bigdog2

	cut_f = cut5*cut6*bigdogcut_f
	cut_c = cut5*cut7*bigdogcut_c

	print 'ctu f ', cut_f
	print 'cut c ', cut_c
	cut = cut_f * cut_c 
	print 'cut ',cut
	print Mformed_ar[cut]
	#plt.axis([0, max(max(AvgSFR_ar), 5), min(min(AvgInflowInner_ar), min(AvgInflowOuter_ar),-5), max(AvgOutflowInner_ar,AvgOutflowOuter_ar 5)])
	plt.xlabel('Stars formed during episode(1e6 Msun)')
	plt.ylabel('Mass expelled during episode (1e6 Msun) ')
	plt.plot(Mformed_ar[cut], Mexpelled_ar[cut] ,'.k', markevery=1)
	plt.plot(Mformed_ar[cut], Mcrossed_ar[cut], '.b', markevery=1)
	#plt.plot(AvgSFR_d_ar, AvgOutflowOuter_ar ,'.g', markevery=1)
	#plt.plot(AvgSFR_d_ar, AvgInflowInner_ar ,'.b', markevery=1)
	#plt.plot(AvgSFR_d_ar, AvgInflowOuter_ar ,'.k', markevery=1)
	plt.legend(['Ejected Flux', 'Ejected Crossed'], loc = 'best')
	plt.savefig(name+'/' + str(today) + 'formed_expelled'+ShellUseInner+ ShellUseOuter +'N'+haloN+'.pdf')
	#plt.show()



	



	###################################################################
	#Figure 8
	#logoutflow.pdf
	###################################################################

	thar3 =['SFR', '0.'+ShellUseInner+'5 Rvir outflow', '0.'+ShellUseOuter+'5 Rvir outflow', 'Mgas > 1cc (1e7 Msun)']
	fig6 =  plt.figure(figsize=(10, 6))
	#plt.axis([min(Ts)*0.95, max(Ts)*1.05, -2, 2])
	plt.plot(Ts, (SFR1+1e-9) , color='k', markevery=1)
	plt.plot(Ts, (Flow_in+1e-9) , color='r', ls='--', markevery=1)
	plt.plot(Ts, (Flow_out+1e-9) , color='r', ls=':',  markevery=1)
	plt.plot(Ts, (DenseGasMass/1e7+1e-9), color='#483d8b', ls='-.')
	plt.ylim(1e-4, 1e2)
	plt.yscale('log')
	#if (uselims == 'y'):
		#plt.ylim(-15, 15)
	#plt.plot(Ts, Flow1 , color='g', ls='dashed', markevery=1)
	plt.ylabel(r'SFR, Outflow, Inflow  ($M_\odot$${\rm yr}^{-1}$)')
	plt.xlabel('cosmic time (Myr)')
	plt.legend(thar3, loc = 'best', prop={'size':9})
	plt.savefig(name+'/' + str(today) + 'logoutflow'  +ShellUseInner+ ShellUseOuter+'N'+haloN+'.pdf')


	###################################################################
	#Figure 9
	#ShellOutflowsN.pdf
	###################################################################

	fig7 =  plt.figure(figsize=(10, 6))
	count = 0
	while (count < interval_count):
		plt.plot(Shells[0:], Flow_total_interval_ar[count][0:],  color=thecolors[count], marker='D')
		count+=1
	plt.yscale('linear')
	plt.xlabel(' R/ Rvir') 
	plt.ylabel('total mass expelled per interval (1e6 Msun)')
	plt.legend(['episode 1','2','3','4','5','6','7','8','9'], loc = 'best', prop={'size':9})
	plt.savefig(name+'/' + str(today) + 'ShellOutflowsN'+haloN+'.pdf')


	###################################################################
	#Figure 10
	#baryonfraction.pdf
	###################################################################

	fig8 = plt.figure(figsize=(10, 6))
	plt.xlabel('cosmic time (Myr)')
	plt.ylabel('fraction')
	fbaruni =  0.044/0.28

	plt.plot(Ts, fbar/fbaruni, color='k', lw=3)
	plt.plot(Ts, fgas/fbaruni, color='b',lw=3)
	plt.plot(Ts, fstar/fbaruni, color='r', lw=3)
	#plt.plot(Ts, (1.0 - cumtotalmasses[9]/cumtotalmasses[9][-1]), color = 'g', ls='-.', lw=3)
	barfolemew = (cumtotalmasses[9] * 1e6) / (TotalMass * fbaruni / little_h)
	plt.plot(Ts, barfolemew, color = 'g', ls=':', lw=2)
	#plt.legend([r'$f_{bar}/f_{bar,uni}$',r'$f_{gas}/f_{bar,uni}$',r'$f_{stars}/f_{bar,uni}$', r'1.0 - $M_{out}$(z)/$M_{out}$(z=2)',  r'$M_{out}$(z) / $f_{bar,uni}M_h$(z)'], loc = 'best', prop={'size':9})
	plt.legend([r'$f_{bar}/f_{bar,uni}$',r'$f_{gas}/f_{bar,uni}$',r'$f_{stars}/f_{bar,uni}$',  r'$M_{out}$(z) / $f_{bar,uni}M_h$(z)'], loc = 'best', prop={'size':9})

	plt.savefig(name+'/' + str(today) + 'baryonfraction'+haloN+'.pdf')

	###################################################################
	#Figure 11
	#outflowfraction.pdf
	###################################################################
	fig9 = plt.figure(figsize=(10, 6))
	plt.xlabel('cosmic time (Myr)')
	plt.ylabel(r'fraction of total mass with $v_{rad} > 0$')
	plt.plot(Ts, TotalOutflow_fraction, lw=3, color='k')
	plt.plot(Ts, Outflow_fraction[0], lw=3, ls=':', color ='#483d8b')
	plt.plot(Ts, Outflow_fraction[2], lw=3, ls='--', color='b')
	plt.plot(Ts, Outflow_fraction[9], lw=3, ls='--', color='r')
	plt.legend(['total','0.05 Rvir','0.25 Rvir', '0.95 Rvir'], loc = 'best', prop={'size':9})

	plt.savefig(name+'/' + str(today) + 'outflowfraction'+haloN+'.pdf')

	###################################################################
	#Figure 12
	#outflowrate.pdf
	###################################################################
	fig10 = plt.figure(figsize=(10, 6))
	plt.xlabel('cosmic time (Myr)')
	plt.ylabel(r'total outflowing mass / dt(snapshot)   ($M_\odot$${\rm yr}^{-1}$)')
	plt.plot(Ts, TheTotalOutflow, lw=3)
	plt.savefig(name+'/' + str(today) + 'total_outflow'+haloN+'.pdf')






Mformed_ar = np.array(Mformed_ar)
Mexpelled_ar = np.array(Mexpelled_ar)
Mcrossed_ar = np.array(Mcrossed_ar)

logMformed_ar = np.log10(np.array(Mformed_ar)+ 1e-10)
logMexpelled_ar = np.log10(np.array(Mexpelled_ar)+ 1e-10)
logMcrossed_ar = np.log10(np.array(Mcrossed_ar) + 1e-10)
bigdog1 = max(Mexpelled_ar)
bigdog2 = max(Mcrossed_ar)
cut5 = logMformed_ar > -8
cut6 = logMexpelled_ar > -8
cut7 = logMcrossed_ar > -8
bigdogcut_f = Mexpelled_ar > 0.01 * bigdog1
bigdogcut_c = Mcrossed_ar > 0.01 * bigdog2

cut_f = cut5*cut6*bigdogcut_f
cut_c = cut5*cut7*bigdogcut_c

print 'ctu f ', cut_f
print 'cut c ', cut_c
cut = cut_f * cut_c 
print 'cut ',cut
print Mformed_ar[cut]

###################################################################
#text writing!
#formed_expelled.txt
#calculate mass loading slope
###################################################################

foutname = (name+'/' + str(today) + 'formed_expelled'+ShellUseInner+ ShellUseOuter +'N'+haloN+'.txt')
g = open(foutname, "w")
count = 0
cut = cut_f * cut_c 

interval_start_ar = np.array(interval_start_ar)
t_ar = np.array(t_ar)
z_ar = np.array(z_ar)
Mvir_ar = np.array(Mvir_ar)
Mstar_ar = np.array(Mstar_ar)
vcirc_ar = np.array(vcirc_ar)
AvgSFR_d_ar = np.array(AvgSFR_d_ar)
AvgOutflowInner_ar = np.array(AvgOutflowInner_ar)
ff = len(Mformed_ar[cut])

while (count < ff):
	#theline = str(interval_start_ar[count])+"  "+str(Mformed_ar[count]) + "  "+str(Mexpelled_ar[count])+"  "+str(t_ar[count])+"   "+str(+z_ar[count])+"  "+str(Mvir_ar[count]) + "  " +str(Mstar_ar[count])+"  "+str(vcirc_ar[count])+"  "+str(AvgSFR_d_ar[count]) +"  "+str(AvgOutflowInner_ar[count])+"  "+ str(Mnew_ar[count])+" \n"[cut]
	theline = str(interval_start_ar[cut][count])+"  "+str(Mformed_ar[cut][count]) + "  "+str(Mexpelled_ar[cut][count])+"  "+str(t_ar[cut][count])+"   "+str(z_ar[cut][count])+"  "+str(Mvir_ar[cut][count]) + "  " +str(Mstar_ar[cut][count])+"  "+str(vcirc_ar[cut][count])+"  "+str(AvgSFR_d_ar[cut][count]) +"  "+str(AvgOutflowInner_ar[cut][count])+"  "+ str(Mcrossed_ar[cut][count])+" \n"

	g.write(theline)
	count += 1
g.close()


cut = cut_f*cut_c

def residuals(b , y, x):
	err = y - (x + b)
	return err
def peval(x, p):
	return x + p[0]
p0 = 5
if (len((logMexpelled_ar[cut]) > 0) and (len(logMformed_ar[cut]) > 0) ):
	plsqtotal = optimize.leastsq(residuals,p0, args=(logMexpelled_ar[cut], logMformed_ar[cut]))
	yintercept_total = plsqtotal[0]

	MLFTX = np.linspace(min(logMformed_ar[cut]), max(logMformed_ar[cut]),num=50)
	MLFTX_lin = pow(10, MLFTX)

	massloadfit_total = MLFTX + yintercept_total
	massloadfit_total = pow(10.0, massloadfit_total)

	MassLoading_slope = pow(10.0, yintercept_total)
else:
	MassLoading_slope = 0

cut = cut_c*cut_f
p0 = 5
if (len((logMcrossed_ar[cut]) > 0) and (len(logMformed_ar[cut]) > 0) ):
	plsqtotal = optimize.leastsq(residuals,p0, args=(logMcrossed_ar[cut], logMformed_ar[cut]))
	yintercept_total = plsqtotal[0]
	MLFTX = np.linspace(min(logMformed_ar[cut]), max(logMformed_ar[cut]),num=50)
	MLFTX_lin = pow(10, MLFTX)

	massloadfit_total = MLFTX + yintercept_total
	massloadfit_total = pow(10.0, massloadfit_total)

	CrossLoading_slope = pow(10.0, yintercept_total)
else:
	CrossLoading_slope = 0






###################################################################
#writing text
#ML_integrated
###################################################################

print 'mass loading '
print 'Shell 2 flux cumulative ', ML_from_total
print 'Shell 2 cross cumulative ',ML_from_total_cross
print 'slope of line ',MassLoading_slope[0]
print 'slope of cross line ',CrossLoading_slope[0]





MLline_list = [int(haloN), MassLoading_slope[0], mend, vend, mstarend, len(logMexpelled_ar[cut]), ML_from_total, ML_from_total_cross, CrossLoading_slope[0], ep_rat_flux, ep_rat_cross, ep_rat_form, paramfac, stagger, zstart, zend, total_expelled_cross, total_gained_cross, totalmasses[2], crossmasses[2], intotalmasses[2], incrosstotalmasses[2], totalmasses[9], crossmasses[9], intotalmasses[9], incrosstotalmasses[9], mstart, vstart, mstarstart]

#totalmasses, crossmasses, intotalmasses, incrosstotalmasses

print MLline_list

MLline = line_to_string(MLline_list)

g = open(name+'/ML_integrated.txt', 'a')
g.write(MLline)
g.close()
#plt.show()

magicmstar = float(getMstar(magicZ))
magicvc = float(getvcirc(magicZ))
magicM = float(getmass(magicZ))

medmstar = float(getMstar(zmedpoint))
medvc = float(getvcirc(zmedpoint))
medM = float(getmass(zmedpoint))

g = open(name+'/N_integrated.txt', 'a')
MLline2_list = [int(haloN),magicZ, magicmstar, magicvc, magicM, zmedpoint, medmstar, medvc, medM ] 
print MLline2_list
MLline2 = line_to_string(MLline2_list)
g.write(MLline2)
g.close()

