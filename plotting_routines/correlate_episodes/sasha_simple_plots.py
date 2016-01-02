import numpy as np
def create_plottable_hist(the_histresults):
	histlim1 = the_histresults[1][:-1]
	histlim2 = the_histresults[1][1:]

	hist_table = (np.array([histlim1, histlim2, the_histresults[0]])).T

	center = (histlim1 + histlim2)/2.0
	width = (center[1]-center[0])

	histx = np.ravel(zip(center,center+width))
	histy = np.ravel(zip(the_histresults[0],the_histresults[0]))
	return [histx, histy]
