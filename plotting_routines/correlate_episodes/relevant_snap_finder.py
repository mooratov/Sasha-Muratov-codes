import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys
today = str(date.today())
print today
from N_for_T_approx_func import N_for_T
import Sasha_functions as SF

Mstarlowlim = 1e0
Mstarupperlim = 5e15

if (len(sys.argv) < 2):
	print 'syntax  blah.py outflow_ep_file'
	sys.exit()

finname1 = str(sys.argv[1])
f = open(finname1)
dars = np.loadtxt(f)
f.close()

SF_start = dars[:,13]
SF_end = dars[:,15]
Mstar = dars[:,8]
eta = dars[:,17]
epN = dars[:,0]

cut1 = Mstar > Mstarlowlim
cut2 = Mstar < Mstarupperlim
cut = cut1*cut2


foutname = today+'relevant_snaps.txt'
f = open(foutname, 'w')
count = 0 
while (count < len(SF_start[cut])):
	Nstart = N_for_T(SF_start[cut][count])
	Nend = N_for_T(SF_end[cut][count])
	print epN[cut][count], Mstar[cut][count]/1e10, Nstart, Nend, eta[cut][count]
	line = [epN[cut][count], Mstar[cut][count]/1e10, Nstart, Nend, eta[cut][count]]
	thestring = SF.line_to_string(line)
	f.write(thestring)
	count+=1
f.close()
