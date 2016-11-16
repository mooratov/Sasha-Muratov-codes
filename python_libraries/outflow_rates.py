import numpy as np
import Sasha_functions as SF

year = 3.15569e7
pc = 3.08567758e18
kpc = pc * 1e3

#prints total outflow flux in all shells
def outflow_flux(GasMass, v_rad, OutflowCut, ShellCuts, ShellThickness_ar, h):
	bincount = 0
	R_bins = len(ShellThickness_ar)
	while (bincount < R_bins):
		ShellThickness = ShellThickness_ar[bincount]
		BinCuts = ShellCuts[bincount]
		outflowtotal = ((GasMass[OutflowCut * BinCuts]*1e10/h) * (v_rad[OutflowCut*BinCuts]*1e5*year/kpc)) / ShellThickness
		sumoutflow = np.sum(outflowtotal)
		print 'outflow at bin ',bincount,' is ', sumoutflow, 'Mun/yr'
		bincount+=1

#returns single outflow rate value
def outflow_flux_atshell(GasMass, v_rad, OutflowCut, Bincut, ShellThickness, h):
	outflowtotal = ((GasMass[OutflowCut * Bincut]*1e10/h) * (v_rad[OutflowCut*Bincut]*1e5*year/kpc)) / ShellThickness
	sumoutflow = np.sum(outflowtotal)
	print 'outflow is ', sumoutflow, 'Mun/yr'
	return sumoutflow


def prepare_X_close(X, Rmax, haloX, haloY, haloZ):
	XcloseHx = abs(X['p'][:,0] - haloX) < Rmax
	XcloseHy = abs(X['p'][:,1] - haloY) < Rmax
	XcloseHz = abs(X['p'][:,2] - haloZ) < Rmax
	#combine all cuts
	#can be used for Stars and DM as well i think
	Xclose = XcloseHx * XcloseHy * XcloseHz  
	return Xclose