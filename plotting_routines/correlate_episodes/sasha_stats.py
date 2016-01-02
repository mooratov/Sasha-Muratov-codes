import numpy as np
import scipy.stats as stats

def sasha_chisquared(expec, observ, degree_of_freedom=2):
    tempnum = (expec - observ)**2
    thevar = stats.tvar(observ)
    MeanSquaredError = np.sum(tempnum) / (len(expec)-2)
    RootMeanSquaredError = np.sqrt(MeanSquaredError)
    tempnum /= thevar
    tempdenom = (len(expec)-degree_of_freedom)
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
