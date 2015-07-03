from .conf import*
from . import GIS
from . import matr
import seaborn as sns

####################################################################################################
def mapDraw(flnm, titlTxt=None, barTxt = None, fignm = None, projMap = 'cea', scale = 1.0, log = False, maskNum = None, vMin=None, vMax=None, cmNm = 'terrain_r', lut = None, split = None):
    from mpl_toolkits.basemap import Basemap
    from numpy import ma
    from matplotlib.colors import LogNorm
        
    inproj = osr.SpatialReference()
    if type(flnm) is str: #GDAL recongnized GIS data
        print 'Input data: GIS dataset.'
        if fignm is None:
            fignm=flnm+'.png'
        if titlTxt is None:
            titlTxt=flnm
        print 'Plotting map from',flnm,'to', fignm
        # Read the data and metadata
        ds = gdal.Open(flnm)
        data = GIS.read(flnm)*scale
        gt = ds.GetGeoTransform()
        proj = ds.GetProjection()
        xres = gt[1]
        yres = gt[5]
        # get the edge coordinates and add half the resolution 
        # to go to center coordinates
        xmin = gt[0] + xres * 0.5
        xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
        ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
        ymax = gt[3] - yres * 0.5
        ds = None
        # create a grid of xy coordinates in the original projection
        inproj.ImportFromWkt(proj)
    elif type(flnm) is np.ndarray:
        print 'Input data: CLM dataset.'
        xmin, xres, xmax = [-180+1.25/2, 1.25, 180-1.25/2]
        ymin, yres, ymax = [-90+180.0/(193*2), -180.0/193, 90-180.0/(193*2)] 
        ##### CLM: 0~360 =====> -180~180
        inproj.ImportFromEPSG(4326)
        data = np.hstack([flnm[:,288/2:],flnm[:,:288/2]])
        del flnm
    else:
        print 'Input data type',type(flnm), 'not recongised.'
        return False
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    
    m = Basemap(projection=projMap,llcrnrlat=ymin,urcrnrlat=min(ymax,90),llcrnrlon=xmin,urcrnrlon=xmax,resolution='c') #extent
    m.drawcountries(linewidth=.1)
    m.drawcoastlines(linewidth=.2)
    m.drawparallels(np.arange(-60,60,30),labels=[1,0,0,0],linewidth=.3,color='grey')
    m.drawmeridians(np.arange(-180,180,30),labels=[0,0,1,0],linewidth=.3,color='grey')
    
    xy_source = np.mgrid[xmin:xmax+xres:xres,ymax+yres:ymin:yres] 
    outproj = osr.SpatialReference()
    outproj.ImportFromProj4(m.proj4string)
    xx, yy = GIS.cordConv(xy_source, inproj, outproj)

    # plot the data (first layer)
    if split is not None:
        lut = len(split)-1
        dump = np.zeros(data.shape) * np.NaN
        for j in range(len(split)):
            dump[data >= split[j]] = j
        data = dump
    colmap = plt.cm.get_cmap(cmNm, lut)
    colmap.set_bad('w')
    if maskNum is not None:
        data[data==maskNum] = np.nan
    
    if log is True:
        if vMin is None:
            vMin = np.max([np.nanmin(data),0.01])
        if vMax is None:
            vMax = 10**np.ceil(np.log10(np.nanmax(data))) 
        im = m.pcolormesh(xx, yy, ma.array(data.T,mask=np.isnan(data.T)), cmap=colmap, norm=LogNorm(vmin=vMin, vmax=vMax))
    else:
        if vMin is None:
            vMin = np.min(data[np.isfinite(data)])
        if vMax is None:
            vMax = np.max(data[np.isfinite(data)])
        im = m.pcolormesh(xx, yy, ma.array(data.T,mask=np.isnan(data.T)), cmap=colmap, vmin=vMin, vmax=vMax)
    
    ax = plt.gca()
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size="5%", pad=0.05)
    cbar = plt.colorbar(im,cax=cax) 
    if split is not None:
        tick_temp = np.arange(len(split)+1)
        cbar.set_ticks(tick_temp)
        string = []
        for l in range(len(split)):
            string.append(str(split[l]))            
        cbar.set_ticklabels(string)

    if barTxt is not None:
        cbar.set_label(barTxt)
    #cbar.update_ticks()
    plt.savefig(fignm, dpi = 500)
    plt.close()

####################################################################################################
def axAdj(ax):
# Save a nice dark grey as a variable
    almost_black = '#262626'
    
    # Remove top and right axes lines ("spines")
    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
    ax.get_yaxis().tick_left() 

    # For remaining spines, thin out their line and change the black to a slightly off-black dark grey
    almost_black = '#262626'
    spines_to_keep = ['bottom', 'left']
    for spine in spines_to_keep:
        ax.spines[spine].set_linewidth(0.5)
        ax.spines[spine].set_color(almost_black)
    
    # Change the labels to the off-black
    ax.xaxis.label.set_color(almost_black)
    ax.yaxis.label.set_color(almost_black)
    
    # Change the axis title to off-black
    ax.title.set_color(almost_black)
    
    # Remove the line around the legend box, and instead fill it with a light grey
    # Also only use one point for the scatterplot legend because the user will 
    # get the idea after just one, they don't need three.
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    light_grey = np.array([float(248)/float(255)]*3)
    legend = ax.legend(frameon=True, scatterpoints=1, fontsize=10, loc='best')#, bbox_to_anchor=(1, 0.5))
    if legend:
        rect = legend.get_frame()
        rect.set_facecolor(light_grey)
        rect.set_alpha(0.9)
        rect.set_linewidth(0.0)
    
    # Change the legend label colors to almost black, too
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)

####################################################################################################
####################################################################################################
def cover_hisPlt():
    his = np.bincount(cover[~(cover==17)])
    freq = np.argsort(his)[::-1]
    
    mpl.use('PDF')
    import matplotlib.pyplot as plt, prettyplotlib as ppl
        
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11, 6))
    ppl.bar(ax, range(len(his)), his, width = 1, grid='y')
    plt.xticks(range(len(his)), sorted(veg_dict , key=veg_dict.get), rotation=30,  fontsize=14)
    plt.ylabel('Pixel Number',  fontsize=14)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title('Histogram of Vegetation Type', fontsize=16)
    plt.tight_layout()
    plt.savefig('afr_cover_his.pdf', dpi=300)
    plt.close()

####################################################################################################
def hist(data, xtext, ttext, figname):
    #mpl.use('PDF')
    import matplotlib.pyplot as plt, prettyplotlib as ppl
    
    fig, ax = plt.subplots(1)
    ppl.hist(ax, data[~np.isnan(data)], bins=20, grid='y')
    plt.ylabel('Pixel Number',  fontsize=14)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel(xtext,  fontsize=14)
    plt.title('Histogram of '+ttext, fontsize=16)
    fig.savefig('afr_'+figname+'_his.pdf', dpi = 300)
    plt.close()

####################################################################################################
def aver2D(x,y,z,xlabel,ylabel,zlabel,flnm):
    binN = 40
    mask = ~(np.isnan(x)|np.isnan(y)|np.isnan(z)|(cover==veg_dict['urban'])|(cover==veg_dict['crop'])|(cover==veg_dict['cryo']))
    H1, _, _ = np.histogram2d(x[mask], y[mask], weights = z[mask], bins = binN)
    H2, xedges, yedges = np.histogram2d(x[mask], y[mask], bins = binN)
    H2[H2==0] = np.nan
    averTree = ma.masked_where(np.isnan(H2), H1/H2)
    
    averTree = np.rot90(averTree)
    averTree = np.flipud(averTree)
    
    mpl.use('PDF')
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    plt.pcolormesh(xedges,yedges,averTree)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(zlabel)   
    plt.tight_layout()
    plt.savefig(flnm, dpi=300)
    plt.close()

####################################################################################################
def scatter(xRaw, yRaw,  nmList=None, figNm=None, nbin=None, divider=None, text=None, percen=95, alpha=0.4, upp = True, med = True):    
    """
    Draw scatter plot in groups, and calculated upper percentile and median. Accept and ignore NaNs.

    Args:
        xRaw
    
    """
    
    if text is None:
        titText, xText, yText = ''
    else: 
        titText, xText, yText = text
    
    if divider is None:
        divider=np.zeros(xRaw.shape, dtype=np.int)
    if not np.nanmax(divider) == len(nmList)-1:
        print 'You sure you have the right length of name list?'
        print 'Maximum in divider:', np.nanmax(divider)
        print 'Length of name list:', len(nmList)
        return

    xRaw, yRaw, divider = matr.cleaner([xRaw, yRaw, divider])
    xNew, yUpp, yMed = matr.binPer(xRaw, yRaw, nbin = nbin, percen = percen)

    if figNm is not None:
        almost_black = '#262626'

        import matplotlib.pyplot as plt, brewer2mpl
    
        clm = brewer2mpl.get_map('Set3', 'qualitative', 8, reverse=True).mpl_colors
    
        fig, ax = plt.subplots(1)
        for i in range(np.nanmax(divider)+1):
            print 'Plotting', nmList[i]
            x = xRaw[divider==i]
            y = yRaw[divider==i]
            color = clm[i%8]
            ax.scatter(x, y, label=nmList[i], alpha=alpha, edgecolor=almost_black, facecolor=color, linewidth=0.15)
            alpha = alpha*.8
        ax.plot(xNew, yUpp, color='black', linewidth=0.7, label='Upper '+ str(percen) +' percentage', alpha=0.6)
        ax.plot(xNew, yMed, color='blue', linewidth=0.7, label='Median', alpha=0.6)
        axAdj(ax)
        ax.set_title(titText+': bin number='+str(xNew.size))
        ax.set_xlabel(xText, fontsize=14)
        ax.set_ylabel(yText, fontsize=14)
        ax.set_xlim([np.nanmin(xRaw)*0.85,np.nanmax(xRaw)*1.1])
        ax.set_ylim([np.nanmin(yRaw)*0.85,np.nanmax(yRaw)*1.1])  
        fig.savefig(figNm, dpi=300)
        fig.close()

    return xNew, yUpp

####################################################################################################
####################################################################################################
def resamPDF(dataList, nmList, valiMask=None, cutList = None, contrEdge = [1000,2000], suff='.png', nbin = 15, sampSize = 1, sampTime = 600):
    #noFireMaskRaw, indeValMaskRaw,
    indeNm, depenNm, contrNm = nmList
    indeVal, depenVal, contrVal = dataList
    del dataList

    if valiMask is None:
        valiMask=np.ones(depenVal.shape,dtype=np.bool)
    valiMask = valiMask&(~np.isnan(depenVal))

    if cutList is None:
        cutList = matr.autoBin(indeVal[valiMask], 4)
        print 'resamPDF: Auto bin boundaries: ',cutList

    SingSamp = lambda x: matr.maskSamp(x, contrVal, minEdge = contrEdge[0], maxEdge = contrEdge[1],  inter = 100, size=sampSize) #do sample for once 
    MultiSamp = lambda mask: np.hstack(depenVal[SingSamp(mask)] for _ in range(sampTime))

    if suff[-4:]=='.pdf':
        mpl.use('PDF')

    import matplotlib.pyplot as plt    
    fig, ax = plt.subplots(1)
    for i in range(len(cutList)-1):
        labelTxt = indeNm+': {:.3f}'.format(cutList[i])+' - {:.3f}'.format(cutList[i+1])
        print labelTxt
        mask = (indeVal>cutList[i])&(indeVal<=cutList[i+1])&valiMask
        print 'Sample size: '+str(mask.sum())
        plt.hist(MultiSamp(mask), histtype='stepfilled', color=plt.cm.jet(i/(len(cutList)-1.0)), alpha=0.3, edgecolor='k', bins=nbin, label=labelTxt, normed=True)
    plt.ylabel('Probability',  fontsize=14)
    plt.xlabel(depenNm,  fontsize=14)
    plt.title(indeNm+' vs '+depenNm+', controlled by '+contrNm, fontsize=16)
    plt.legend(loc='upper center')

    fig.savefig('resamPDF'+indeNm+depenNm+str(contrEdge[0])+'to'+str(contrEdge[1])+contrNm+suff, dpi = 300)
    plt.close()

####################################################################################################
def treeCompare():
    tempPath = outPath+'LandTypeScore/'
    treeDict = {'treeMOD44': outPath + 'MOD44B.tif', 'treeMCD12PFT': outPath+'LandTypeScore/MCD12PFTForestScore.tif', 'treeMCD12IGBP': outPath+'LandTypeScore/MCD12IGBPForestScore.tif', 'treeGLC': outPath+'LandTypeScore/GLCForestScore.tif'}
    dataList = matr.cleaner([GIS.read(treeDict[dataName]) for dataName in treeDict], NaNCut=False, standard=80, scalingPoint=100)
    
    i = 0
    for dataName in treeDict:
        flnm = tempPath+dataName+'Rescale.tif'
        GIS.write(dataList[i], flnm, treeDict[dataName])
        mapDraw(flnm, tempPath+dataName+'Rescale.png', dataName+' (rescaled)')
        i += 1

####################################################################################################
def explorer(data, name, hue=None, trel=True, corr=True):
    """        
    Draw and save Trellis plots including scatter plots (upper triangle) and kernal density (lower triangle and lower triangle), correlation map with person R and p value. Takes long time with big data.
    
    Args:
        data: dataFrame. Input data arrays.
        name: str. Name of output figure file. 
        hue: str, optional. Name of variable used as hue. 
    Return:
        PairGrid
    """

    if name[-4:]=='.pdf': 
        mpl.use('PDF')
    import matplotlib.pyplot as plt
    #sns.set_context("talk", font_scale=1.3)
    if trel:
        print 'Plotting Trellis plots.'
        #sns.set(style="white")
        #f, ax = plt.subplots(figsize=(7, 7))
        #ax.set(xscale="log", yscale="log")
        g = sns.PairGrid(data, hue=hue)
        g.map_lower(sns.kdeplot, cmap="Purples",shade=True)
        g.map_diag(plt.hist)
        g.map_upper(plt.scatter, s=10, alpha=.05)
    
        g.savefig('trel_'+name, dpi = 300)
        plt.close()

    if corr:
        print 'Plotting correlation map.'
        #sns.set_context(rc={"figure.figsize": (16, 16)})
        plt.figure()
        ax = sns.corrplot(data)
        ax.figure.savefig('corr_'+name, dpi = 300)
        plt.close()    
    
    
####################################################################################################
def bivar(data, name, var1=None, var2=None):
    """        
    Draw and save bivariant data, with linear regression and distribution.

    Args:
        data: dataFrame. Input data arrays.
        var1: str or list. Varibale name in data, plotted on x-axis. Use all varibale in data if not set.
        var2: str or list. Varibale name in data, plotted on y-axis. Use all varibale in data if not set.
        name: str. Name (suffix) of output figure file. 
    Return:
        PairGrid
    """

    from scipy import stats

    if name[-4:]=='.pdf': 
        mpl.use('PDF')
    import matplotlib.pyplot as plt
    
    def organize(var):
        if type(var) is str:
            var = [var]
        if var is None:
            var = data.columns
        return var
    var1 = organize(var1)
    var2 = organize(var2)
 
    for i in range(len(var1)):
        var1s = var1[i]
        for var2s in var2[i:]:
            if var1s==var2s:
                continue
            print 'Plotting relationship between', var1s, 'and', var2s 
            g = sns.JointGrid(var1s, var2s, data)
            g.plot_marginals(sns.distplot, color="seagreen")
            g.plot_joint(sns.regplot, color=".5", scatter_kws={'edgecolor':'white', 's':8, 'alpha':.4})
            g.annotate(stats.pearsonr)
            plt.savefig(var1s+'_'+var2s+'_'+name, dpi=300)
            plt.close()

####################################################################################################
def LarsPlot(X, y, figNm, alpha=1, titleTxt=''):    
    """
    Draw scatter plot and regression line computed by Lasso model fit with Least Angle Regression.

    Args:
        X
        y
        alpha
    
    """
    from sklearn import linear_model
    clf = linear_model.LassoLars(fit_path=False, alpha=alpha)
    clf.fit(X,y)
    coef = clf.coef_
    inter = clf.intercept_
    Rsq = clf.score(X, y)
    pred = clf.predict(X)
    print 'Coefficient:', coef
    print 'Interception:', inter
    print 'R square:', Rsq

    almost_black = '#262626'
    import matplotlib.pyplot as plt, brewer2mpl
    
    clm = brewer2mpl.get_map('Set3', 'qualitative', 8, reverse=True).mpl_colors
    
    fig, ax = plt.subplots(1)
    print 'Plotting', figNm
    color = clm[1]
    ax.scatter(pred, y, alpha=0.1, edgecolor=almost_black, facecolor=color, linewidth=0.15)
    
    limMin = np.min([pred.min(),y.min()])
    limMax = np.max([pred.max(),y.max()])
    ax.plot([limMin,limMax], [limMin,limMax], color='black', linewidth=0.7, alpha=0.6)
    axAdj(ax)
    ax.set_title(titleTxt+r'. $\alpha=$'+str(alpha)+r'$, R^2=$ %0.2f' % (Rsq,))
    ax.set_xlabel('Prediction', fontsize=14)
    ax.set_ylabel('Observation', fontsize=14)
    ax.set_xlim([limMin,limMax])
    ax.set_ylim([limMin,limMax])  
    fig.savefig(figNm, dpi=300)
