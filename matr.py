from conf import*
import pandas as pd

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
def reBin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return np.nanmean(a.reshape(sh),axis=-1).mean(1)

####################################################################################################
def cleaner(dataList, mask=None, nameList=None, NaNCut=True, scalingPoint=None):
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
##########Calculate common mask and standard scale    
    if NaNCut:
        comMask = ~np.isnan(dataList[0])    
        for i in range(1,len(dataList)):
            comMask = comMask & (~np.isnan(dataList[i]))
    
    if scalingPoint is not None:
        percenF = lambda x: np.percentile(x[~np.isnan(x)], scalingPoint) #calculate percentile
        percenList = [percenF(dataList[i]) for i in range(len(dataList))]
        standard = np.max(percenList)
        print 'Scaling standard:', scalingPoint, 'percentile =', standard
        for i in range(len(dataList)):
            dataList[i] = dataList[i]*standard/percenList[i]

##########Apply additional mask if needed
    if mask is not None:
        if NaNCut:
            comMask = comMask & (~mask)
        else:
            for i in range(len(dataList)):
                dataList[i][mask] = np.nan 

    if nameList is not None: #return pandas dataFrame
        if len(dataList)!=len(nameList):
            print 'Cleaner: Wrong length for variable name.'
        else:
            return pd.DataFrame({nameList[i]: dataList[i][comMask] for i in range(len(dataList))})
    elif NaNCut: #return list of 1-d arrays
        return [ dataList[i][comMask] for i in range(len(dataList)) ]
    else: #return list of origional shape arrays
        return dataList

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
            print str(edge)+' to '+str(edge+inter)
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
        nbin = np.max((xOrd.max()//150, 5))
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
def autoBin(data, nBin, log=False):
    nanMask = ~np.isnan(data)
    nData = nanMask.flatten().sum()
    widBin = nData//nBin
    dataSort = np.sort(data[nanMask])
    if log:
        dataSort[dataSort==0] = dataSort.min()/2.0
        dataSort = np.log(dataSort)
    bound = dataSort[0::widBin]
    if len(bound)<nBin+1:
        bound = np.append(bound, dataSort[-1])
    else:
        bound[-1] = dataSort[-1]
    if log:
        bound = np.exp(bound)
    return bound
