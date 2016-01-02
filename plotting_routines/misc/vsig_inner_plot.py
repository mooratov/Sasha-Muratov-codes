import numpy as np
import Sasha_functions as SF
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import date

today = date.today()

Nhalo = 0

H = SF.read_halo_history(Nhalo)

Zs = H['redshift']
IDs = H['ID']
Vsig = H['Vsig']
M = H['M']
Rvir = H['Rvir']
Vmax = H['Vmax']
newtonG = 6.67384e-8
little_h = 0.7
As = 1.0 / (1.0 + Zs)
Vcirc = np.sqrt( (newtonG * M * 2e33 /little_h )/(Rvir * As * 3.08e21 / little_h)) / 1e5
Vsig1D = Vsig / np.sqrt(3)

Ns = []
for theZ in Zs:
	theN = SF.N_for_z(theZ)
	Ns.append(theN)
Ns = np.array(Ns)

count = 0
insig_ar = []
halfsig_ar = []
thirdsig_ar = []

M_1_ar = []
M_quart_ar = []

Vc_1_ar = []
Vc_quart_ar = []

Pot_1_ar = []
Pot_quart_ar  = []

vmax_ar = []
rmax_ar = []

while (count < len(Zs)):
	if (int(Ns[count]) < 10):
		Nsnapstring = '00'+str(Ns[count])
	elif (int(Ns[count]) < 100):
		Nsnapstring = '0'+str(Ns[count])
	else:
		Nsnapstring = str(Ns[count])
	redshiftstring = "{0:.3f}".format(Zs[count])
	prof = SF.read_profile(Nsnapstring, redshiftstring, int(IDs[count]))
	#print prof
	Sig = prof[:,7]
	rs = prof[:,0] * As[count] / little_h
	kpc_Ms = prof[:,2] / little_h
	kpc_vcircs = prof[:,5] 
	
	print Ns[count], rs
	get_sig = interp1d(rs, Sig, kind='linear')
	get_1kM = interp1d(rs, kpc_Ms, kind='linear')
	get_1kvcirc = interp1d(rs, kpc_vcircs, kind='linear')
	
	vmax = np.max(kpc_vcircs)
	rmax = rs[np.where(kpc_vcircs == vmax)][0]
	
	vmax_ar.append(vmax)
	rmax_ar.append(rmax)
	
	
	minrad = max([1.0, rs[0]])
	thesig1kpc = get_sig(minrad)
	theM1kpc = get_1kM(minrad)
	theVc1kpc = get_1kvcirc(minrad)
	thePot1kpc = (theVc1kpc * 1e5)**2
	
	
	
	thesig_halfRvir = get_sig(0.5*Rvir[count]* As[count]/little_h)
	thesig_thirdRvir =  get_sig((1.0/3.0)*Rvir[count]* As[count]/little_h)
	
	quartR = 0.25*Rvir[count]* As[count]/little_h
	theM_quart = get_1kM(quartR)
	theVc_quart = get_1kvcirc(quartR)
	thePot_quart = (theVc_quart* 1e5)**2
	
	insig = Sig[0]
	insig_ar.append(thesig1kpc)
	halfsig_ar.append(thesig_halfRvir)
	thirdsig_ar.append(thesig_thirdRvir)
	
	M_1_ar.append(theM1kpc)
	Vc_1_ar.append(theVc1kpc)
	Pot_1_ar.append(thePot1kpc)
	
	M_quart_ar.append(theM_quart)
	Vc_quart_ar.append(theVc_quart)
	Pot_quart_ar.append(thePot_quart)
	
	
	
	
	
	print insig, thesig1kpc
	count +=1
	

insig_ar = np.array(insig_ar)

masscolumns = np.column_stack([Ns, Zs, Vsig, Vcirc, Vmax, halfsig_ar, thirdsig_ar])

extracols = np.column_stack([Ns, Zs, M_1_ar, Vc_1_ar, Pot_1_ar, M_quart_ar, Vc_quart_ar,Pot_quart_ar, vmax_ar, rmax_ar])


foutname = ('./' + str(today) +'Velocity_props.txt')
np.savetxt(foutname, masscolumns, fmt='%1.6g')

foutname = ('./' + str(today) +'one_and_quart.txt')
np.savetxt(foutname, extracols, fmt='%1.6g')



fig1 = plt.figure(figsize=(10,8))
plt.plot(Zs, Vsig, '.k')
plt.plot(Zs, Vcirc, '.b')
plt.plot(Zs, Vmax, 'ob')
plt.plot(Zs, halfsig_ar, '.r')
plt.plot(Zs, thirdsig_ar, 'or')
plt.gca().invert_xaxis()
plt.legend(['Vsig', 'Vcirc', 'Vmax', 'Vsig 0.5 Rvir',  'Vsig 0.33 Rvir', ], loc = 'best')
plt.xlabel('z', fontsize=20)
plt.ylabel('km/s', fontsize=20)
foutname = str(today)+'blahR.pdf'
plt.savefig(foutname)

plt.clf()
fig2 = plt.figure(figsize=(10,8))
plt.plot(Zs, M_1_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('M (Msun)', fontsize=20)
plt.yscale('log')
foutname = str(today)+'M1kpc.pdf'
plt.savefig(foutname)

plt.clf()
fig3 = plt.figure(figsize=(10,8))
plt.plot(Zs, M_quart_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('M (Msun)', fontsize=20)
plt.yscale('log')
foutname = str(today)+'Mquart.pdf'
plt.savefig(foutname)

plt.clf()
fig4 = plt.figure(figsize=(10,8))
plt.plot(Zs, Vc_1_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('km/s', fontsize=20)
foutname = str(today)+'vc1kpc.pdf'
plt.savefig(foutname)

plt.clf()
fig5 = plt.figure(figsize=(10,8))
plt.plot(Zs, Vc_quart_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('km/s', fontsize=20)
foutname = str(today)+'vcquart.pdf'
plt.savefig(foutname)

plt.clf()
fig6 = plt.figure(figsize=(10,8))
plt.plot(Zs, Pot_1_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('Pot', fontsize=20)
plt.yscale('log')
foutname = str(today)+'pot1kpc.pdf'
plt.savefig(foutname)

plt.clf()
fig7 = plt.figure(figsize=(10,8))
plt.plot(Zs, Pot_quart_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('Pot', fontsize=20)
plt.yscale('log')
foutname = str(today)+'potquart.pdf'
plt.savefig(foutname)

plt.clf()
fig8 = plt.figure(figsize=(10,8))
plt.plot(Zs, vmax_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('vmax (km/s)', fontsize=20)
plt.yscale('log')
foutname = str(today)+'vmax.pdf'
plt.savefig(foutname)

plt.clf()
fig9 = plt.figure(figsize=(10,8))
plt.plot(Zs, rmax_ar, '.k')
plt.gca().invert_xaxis()
plt.xlabel('z', fontsize=20)
plt.ylabel('rmax (kpc)', fontsize=20)
plt.yscale('log')
foutname = str(today)+'rmax.pdf'
plt.savefig(foutname)
