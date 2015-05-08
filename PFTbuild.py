from pyEarth.conf import *
from osgeo import gdal

####################################################################################################
def stepScale(dataRaw, target, sumAxis=0):
    '''
    data = dataRaw.astype(np.float64)
    space = target.astype(np.float64)
    for i in range(data.shape[sumAxis]):
        ratio = space/data[i:,...].sum(axis=sumAxis)
        ratio[np.isinf(ratio)] = 0
        if (ratio<0).sum()>0:
            print 'Something is wrong. Negative ratio in stepScale.'
        data[i,...] = data[i,...]*ratio
        space -= data[i,...]
    return data
    '''
    data = np.empty(dataRaw.shape)
    total = dataRaw.sum(axis=sumAxis)
    #ratio[np.isinf(ratio)] = target[np.isinf(ratio)]/data.shape[sumAxis]
    for i in range(data.shape[sumAxis]):
        ratio = dataRaw[i,...]/total
        data[i,...] = target*ratio
        data[i,...][total==0] = target[total==0]/data.shape[sumAxis]
    check = (np.abs(data.sum(axis=sumAxis)-target)>1e-8).sum()
    if check==0:
        print 'stepScale self-check passed.',
    else:
        print 'stepScale self-check not passed. Check number:',check

    return data

####################################################################################################
def getPreci():
    ##### Get precipitation data for CLM
    ###gdalwarp -te -180 -90 180 90 -tr 1.25 0.9375 -r average -overwrite 3B43.tif 3B43_lowRe.tif

    ds = gdal.Open(dataPath+'TRMM/3B43/3B43_lowRe.tif')
    preciRAW = ds.ReadAsArray()
    preci = np.empty(preciRAW.shape)
    preci[:,:288/2]=preciRAW[:,288/2:]
    preci[:,288/2:]=preciRAW[:,:288/2]
    preci[preci<0] = np.nan

    #CLMmap(preci, titlTxt = 'Precip in CLM', fignm = 'preciCLM.png')
    del preciRAW

    return preci

####################################################################################################
class PFT(object):
    rowNum = 192
    colNum = 288
    pixH = 180.0/rowNum
    pixW = 360.0/colNum

    lon, lat = np.meshgrid(np.arange(pixW/2, 360, pixW), np.arange(90-pixH/2, -90, -pixH))
    lon += 1.25

    mask = (lat>-23.25)&(lat<23.25) #tropical mask

    ################################################################################################  
    def __init__(self,flnm='surfdata_0.9x1.25_simyr1850_c130415.nc'):
        ds = gdal.Open('NETCDF:"'+flnm+'":PCT_PFT')
        self.source = flnm #filename of surface data
        self.pft = ds.ReadAsArray() #original pft data
    
        ds = gdal.Open('NETCDF:"'+flnm+'":PFTDATA_MASK')
        self.pftMask = ds.ReadAsArray() #pft mask

        ds = gdal.Open('NETCDF:"'+flnm+'":AREA')
        self.area = ds.ReadAsArray()

        self.landFrac = self.pft.sum(axis=0) #land fraction
        self.tree = self.pft[1:9,:,:].sum(axis=0) #total tree cover

    ################################################################################################  
    def boost(self, flnm='/data8/data/guol3/cutboard/MOD44B.boost.lowRe.tif'):
        from pyEarth import plot, GIS        
        
        '''
        self.preci = getPreci()
        
        self.mask = (self.mask&(self.preci>0))&(self.tree>0)
        preUpp, treeUpp = plot.scatter(self.preci[self.mask], self.tree[self.mask], ['PFT data in CLM'], 'preciCLMscatter.png', text=['Precipitation vs Tree Cover','Precipitation (mm/year)','Tree Cover (%)'])
        treeUpp = np.insert(treeUpp, 0, 0)
        preUpp = np.insert(preUpp, 0, 0)
        treeUpp[preUpp>2100] = treeUpp.max()
        treeUpp = treeUpp*100.0/treeUpp.max()
        '''
        
        treePoRaw = GIS.read(flnm)/0.8

        treePoHaf = np.empty(treePoRaw.shape)
        treePoHaf[:,:288/2]=treePoRaw[:,288/2:]
        treePoHaf[:,288/2:]=treePoRaw[:,:288/2]
       
        treePo = np.maximum(treePoHaf,self.tree)
        del treePoRaw, treePoHaf

        self.treePo = self.tree.copy()
        self.treePo[self.mask] = np.minimum(treePo,self.landFrac-self.pft[15,...])[self.mask]
        self.treePo[self.pftMask==0] = 0
        self.treePo[self.tree==0] = 0

        #####Initialize potential PFTs
        self.pftPo = self.pft.copy()
        #####Boost individual tree PFTs
        self.pftPo[1:9,...] = stepScale(self.pft[1:9,...], self.treePo)
        #####Maintain crops
        self.pftPo[15,...] = self.pft[15,...]
        #####Shrink non-tree PFTs
        self.nonTreePo = self.landFrac - self.pft[15,...] - self.treePo
        if (self.nonTreePo<-1e-8).sum()>0:
            print 'Oppppps, nonTreePo got negative.'
        li = range(9,15); li.extend([0])
        self.pftPo[li,...] = stepScale(self.pft[li,...], self.nonTreePo)
        
        print 'Sum of potential PFTs smaller than 0:', (self.pftPo<-1e-8).sum(axis=2).sum(axis=1)
        #test = (self.treePo + self.nonTreePo + self.pftPo[15,...])
        #plot.mapDraw(test, 'Sum of all PFTs', 'test.png')

        return self.pftPo
    
    ################################################################################################ 
    def layerMap(self):
        try:
            if self.pftPo is None:
                self.boost()
        except AttributeError:
            self.boost()

        from pyEarth import plot
        nameList = ['BareLand','NeedleEverTemp','NeedleEverBoreal','NeedleDeciduBoreal','BroadEverTrop','BroadEverTemp','BroadDeciduTrop','BroadDeciduTemp','BroadDeciduBoreal','BroadEverShrubTemp','BroadDeciduShrubTemp','BroadDeciduShrubBoreal','C3Arctic','C3','C4','C3Crop']        
        for i in range(len(nameList)):
            plot.mapDraw(self.pft[i,...], 'Original '+nameList[i], nameList[i]+'Old.png')
            plot.mapDraw(self.pftPo[i,...], 'New '+nameList[i], nameList[i]+'New.png')
        plot.mapDraw(self.pft[:len(nameList),...].sum(axis=0), 'Original PFT sum', 'PFTsumOld.png')
        plot.mapDraw(self.pftPo[:len(nameList),...].sum(axis=0), 'New PFT sum', 'PFTsumNew.png')
        plot.mapDraw(self.tree, 'Original Tree Cover', 'treeOld.png')
        plot.mapDraw(self.treePo, 'New Tree Cover', 'treeNew.png')

    ################################################################################################ 
    def flush(self):      
        try:
            if self.pftPo is None:
                self.boost()
        except AttributeError:
            self.boost()
        
        dst = self.source[:-3]+'.boosted.nc'
        from shutil import copyfile
        copyfile(self.source, dst)
        from netCDF4 import Dataset
        nc = Dataset(dst, mode="a")
        nc.variables['PCT_PFT'][:] = self.pftPo[:,::-1,:]
        nc.sync()
        nc.close()

####################################################################################################

Ipft=PFT('surfdata_0.9x1.25_simyr1850_c130415.nc')
#Ipft.boost()
#Ipft.flush()
#Ipft.layerMap()
