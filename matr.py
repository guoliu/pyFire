from configure import*
from . import GIS
import pandas as pd

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
def reBin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return np.nanmean(a.reshape(sh),axis=-1).mean(1)

####################################################################################################
def cleaner(dataList, mask=None, nameList=None, nanCut=True):
    """
    Through away the NaNs in multiple dataset with the same size, and/or output to Panda DataFrame if *nameList* is provided.
    
    Args:
        dataList: list or tuple. Input data arrays.
        mask: matrix of bool, optional. Mask of invalid data.
        nameList: list of str, optional. List of columns names. Output to DataFrame if provided. Otherwise output is a list of arrays.
        nanCut: bool, optinal. Through away NaNs.
    
    Return:
        list (when *nameList* not provided) or DataFrame (when *nameList* is provided).
    """

    comMask = np.ones(dataList[0].shape, dtype=np.bool)
    if nanCut:
        for i in range(len(dataList)):
            comMask = comMask & (~np.isnan(dataList[i]))

    if mask is not None:
        comMask = comMask & (~mask)
    
    if nameList is not None:
        if len(dataList)!=len(nameList):
            print 'Cleaner: Wrong length for variable name.'
        else:
            return pd.DataFrame({nameList[i]: dataList[i][comMask] for i in range(len(dataList))})
    return [ dataList[i][comMask] for i in range(len(dataList)) ]

####################################################################################################
def arraySamp(array, size = 50):
    pool = np.array(np.nonzero(array)).T
    if np.size(pool,0) < size:
        print 'Not enough sample. Sample size: ', np.size(pool,0)
        return None
    ind = np.array(random.sample(pool, size))
    new_array = np.zeros(array.shape,dtype=bool)
    new_array[ind[:, 0], ind[:, 1]] = True
    return new_array

####################################################################################################
def maskSamp(mask, contr, minEdge = None, maxEdge = None, inter = 100, size = 50):
    #sample mask layer according to contr to avoid distribution error in data
    new_mask = np.zeros(mask.shape,dtype=bool)
    if minEdge is None:
        minEdge = contr[mask].min().astype(int)
    if maxEdge is None:
        maxEdge = (contr[mask].max()+0.999*inter).astype(int)
    
    for edge in range(minEdge, maxEdge, inter):
        temMa = arraySamp((contr<=edge+inter)&(contr>=edge)&mask, size = size)
        if temMa is None:
            print 'Reduce sample size at edge: ', edge
            return None
        new_mask = new_mask|temMa    
    return new_mask

####################################################################################################
def quaReg(x, y):
    import rpy2.robjects as robjects 
    rObj = robjects.r
    qreg = robjects.packages.importr(name = 'quantreg')
    
    rObj.assign('xr', robjects.FloatVector(x))
    rObj.assign('yr', robjects.FloatVector(y))
    rObj("dt <- data.frame(xr,yr)")
    fmla = robjects.Formula('')

####################################################################################################
def binPer(x, y , nbin = None, percen = 90, med = True):
    print 'Calculating upper boundary. Percentage:', percen
    from accumarray import accum
    from scipy import stats
    x = x.flatten()
    y = y.flatten()
    
    xOrd = stats.rankdata(x,method='dense')
    if nbin is None:
        nbin = np.max((xOrd.max()//100, 5))
        print 'Automatic bin number:', nbin
    else:
        print 'User assigned bin number:', nbin
    binWid = xOrd.max()*1.0/nbin 
    mapper = (xOrd//binWid).astype(np.int)    

    yNew = accum(mapper, y, func=lambda x: np.percentile(x[:], percen))
    xNew = accum(mapper, x, func='mean')
    
    if med:
        return xNew, yNew, accum(mapper, y, func=lambda x: np.percentile(x[:], 50))
    else:
        return xNew, yNew

####################################################################################################
def reMap(xUpp, yUpp, xBas, template = None, flnm = None):
    xMask = ~np.isnan(xBas)
    
    yPo = np.empty(xBas.shape)
    yPo[:] = np.nan
    yPo[xMask] = np.interp(xBas[xMask], xUpp, yUpp)

    if template is not None and flnm is not None:
        GIS.write(yPo, flnm, template)

    return yPo


