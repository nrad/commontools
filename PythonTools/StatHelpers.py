from uncertainties import ufloat, ufloat_fromstr


########################################################################################################
###################     ufloat helpers
########################################################################################################


def calcWeightedMean( uf_list ):
    vals = np.array([getVal(x)   for x in uf_list])
    sigs = np.array([getSigma(x) for x in uf_list])
    weights = 1/sigs**2
    mean = np.average( vals, weights=weights)
    return mean



def calcWeightedMeanAndStandardDeviation( uf_list, true_mean=None ):
    """
        Returns the average of the sample which each value weighted by 1/sig^2
        The error is returned as the standard deviation
        
    """
    vals = np.array([getVal(x) for x in uf_list])
    n    = len(uf_list)
    if true_mean==None:
        mean = calcWeightedMean(uf_list)
        std = np.sqrt( sum( 1/(n-1 ) * ((vals-mean)**2) ) )
    else:
        mean = true_mean
        std = np.sqrt( sum( 1/(n ) * ((vals-mean)**2) ) )
        
    return ufloat(mean,std)



def calcWeightedMeanAndError( uf_list ):
    """
        Returns the average of the sample which each value weighted by 1/sig^2
        The error is the uncertainty of the weighted average
        https://physics.stackexchange.com/questions/15197/how-do-you-find-the-uncertainty-of-a-weighted-average
    """
    vals = np.array([ getVal(x) for x in uf_list])
    sigs = np.array([ getSigma(x) for x in uf_list])
    weights = 1/sigs**2
    mean  = np.average( vals, weights=weights)
    error = 1/sum(weights) * np.sqrt( sum( (weights*sigs)**2  ) )
    return ufloat(mean,error)



def getVal(uf, strict=True):
    if hasattr(uf, 'nominal_value'):
        return uf.nominal_value
    elif hasattr(uf, 'val'):
        return uf.val
    else:
        if strict:
            raise ValueError("Can't tell if this is a ufloat!")
        else:
            return uf

def getSigma(uf, strict=True, def_val=None):
    if hasattr(uf, 'std_dev'):
        return uf.std_dev
    elif hasattr(uf, 'sigma'):
        return uf.sigma
    else:
        if strict:
            raise ValueError("Can't tell if this is a ufloat!")
        else:
            return def_val


def calcWeightedMean( uf_list ):
    vals = np.array([getVal(x)   for x in uf_list])
    sigs = np.array([getSigma(x) for x in uf_list])
    weights = 1/sigs**2
    mean = np.average( vals, weights=weights)
    return mean



def calcWeightedMeanAndStandardDeviation( uf_list, true_mean=None ):
    """
        Returns the average of the sample which each value weighted by 1/sig^2
        The error is returned as the standard deviation
        
    """
    vals = np.array([getVal(x) for x in uf_list])
    n    = len(uf_list)
    if true_mean==None:
        mean = calcWeightedMean(uf_list)
        std = np.sqrt( sum( 1/(n-1 ) * ((vals-mean)**2) ) )
    else:
        mean = true_mean
        std = np.sqrt( sum( 1/(n ) * ((vals-mean)**2) ) )
        
    return ufloat(mean,std)



def calcWeightedMeanAndError( uf_list ):
    """
        Returns the average of the sample which each value weighted by 1/sig^2
        The error is the uncertainty of the weighted average
        https://physics.stackexchange.com/questions/15197/how-do-you-find-the-uncertainty-of-a-weighted-average
    """
    vals = np.array([ getVal(x) for x in uf_list])
    sigs = np.array([ getSigma(x) for x in uf_list])
    weights = 1/sigs**2
    mean  = np.average( vals, weights=weights)
    error = 1/sum(weights) * np.sqrt( sum( (weights*sigs)**2  ) )
    return ufloat(mean,error)




