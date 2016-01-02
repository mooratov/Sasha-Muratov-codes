import numpy as np

def tasz(omeganot, littleh, interval=1000):
    amin = 1.0 / interval
    #amin = 0.001
    amax = 1.0 
    alist = []
    tlist = []
    for i in range (0, (interval+1)):
        a = amin + i*(amax-amin)/interval
        f = omeganot * a ** (-3.0) / (1.0 - omeganot)
        t = (2.0/3.0)* 9.779e9 * (1.0/littleh) * np.log((1.0 + (1.0+f)**0.5)/f**0.5) / (1.0-omeganot)**0.5
        print a,'  ',t/1e9
        alist.append(a)
        tlist.append(t/1e9)
    return alist, tlist

def tfora(a, omeganot, littleh):
    f = omeganot * a ** (-3.0) / (1.0 - omeganot)
    t = (2.0/3.0)* 9.779e9 * (1.0/littleh) * np.log((1.0 + (1.0+f)**0.5)/f**0.5) / (1.0-omeganot)**0.5
    print (t/1e9)
    return (t/1e9)
    
def hubble_param(a, omeganot, littleh):
	omega_lamda = 1.0 - omeganot
	z = 1.0 / a - 1.0
	H = np.sqrt(omeganot * (1.0 + z)**3.0 + omega_lamda) * 100.0 * littleh
	return H
	
def critical_density(H_in_cgs):
	#assumes H given in cm / s / Mpc
	#thats right , cm not km 
	newtonG = 6.67e-8
	rho_crit = 3.0 * H_in_cgs * H_in_cgs / (8.0 * 3.14159* newtonG)
	#for return in g / cm^3
	rho_crit /= (3.08e24)**2.0
	return rho_crit

def baryonic_matter_density(rho_crit, omega_matter, fbar):
	bar = rho_crit * omega_matter * (f_bar/omega_matter)
	return bar

def compute_OD(z, omega_m):
	omega_l = 1.0 - omega_m
	omega_z = omega_m * pow((1.0+z), 3.0) / (  omega_m * pow((1.0+z), 3.0) + omega_l)
	x = omega_z - 1.0
	delta_c = 18 * np.pi* np.pi + 82*x - 39*x*x
	return delta_c

def virial_temperature(Mvir, Rvir):
	#assumes M  in solar masses
	#asumes R in physical kpc
	theTemp = (2.0 / 3.0) * (6.67e-8 * Mvir*2e33) / (Rvir*3.08e21)
	typical_mu = 0.6329113924050632 #this is calculated assuming electron fraction of 1 (fully ionized mixture of hydrogen, helium, metals).
	protonmass = 1.6726e-24 #in g
	BOLTZMANN  = 1.3806e-16 #in cgs
	
	theTemp*= typical_mu * protonmass / BOLTZMANN
	return theTemp
	
	
def compute_Rvir_for_Mvir(M, rho_crit, z, omega_m):
	#assumes M is in Solar Masses
	#rho crit in g/cm^3
	Rcubed = (M * 2.e33) / (compute_OD(z, omega_m) * rho_crit * 4.0 * np.pi / 3.0)
	R = Rcubed ** (1.0/3.0)
	return (R/3.08e21) #return in kpc
	
def compute_vcirc_for_Mvir(M, rho_crit, z, omega_m):
	vcsquared = 6.67e-8 *  (M * 2.e33)**(2.0/3.0) * ((4.0/3.0) * np.pi * rho_crit * compute_OD(z, omega_m))**(1.0/3.0)
	vc = np.sqrt(vcsquared)
	return vc/1e5
	#return (R/3.08e21) #return in kpc	
