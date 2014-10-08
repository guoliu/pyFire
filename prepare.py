from analy import*
from . import GIS

####################################################################################################
### return iteration struture of datetime with given interval
def timerange(startY, endY, inter):
    import datetime
    if inter is '16day':
        return (datetime.date(y,1,1)+datetime.timedelta(d-1) for y in range(startY, endY) for d in range(1,366,16)) 
    if inter is 'month':
        return (datetime.date(y, m, 1) for y in range(startY, endY) for m in range(1, 13))
    if inter is 'year':
        return (datetime.date(y, 1, 1) for y in range(startY, endY))
    print 'Unsupported interval type'


def timeAver(oldPref, oldSuff, newFile, inter, method='average', start = startY, end = endY, dtype=gdal.GDT_Float32):
       
    nmlist = [oldPref + dp.strftime('%Y') + dp.strftime('%j') + oldSuff for dp in timerange(start, end, inter) if os.path.isfile(oldPref + dp.strftime('%Y') + dp.strftime('%j') + oldSuff)]
    timeSpan = len(nmlist) #total time series length
    
    l = 0
    for oldFile in nmlist:
        print oldFile
        try:
            g = gdal.Open(oldFile,gdal.GA_ReadOnly)
            dataName = g.GetSubDatasets()[0][0]
        except IndexError:
            dataName = oldFile                      
        dataset = gdal.Open(dataName, gdal.GA_ReadOnly)
        band = dataset.GetRasterBand(1)
        dataArr = band.ReadAsArray()
        if l==0:
            dataTotal = np.empty((dataArr.shape[0],dataArr.shape[1],timeSpan))
            dataTotal[:] = np.nan
        dataTotal[:,:,l] = dataArr.astype(np.float)
        l += 1
        
    noData = band.GetNoDataValue()
    print 'Fill value: ', noData
    print 'Calculating and writing data to ', newFile 
    if method is 'average':
        dataTotal[dataTotal==noData] = np.nan
        result = np.nanmean(dataTotal[:,:,:l], axis = 2)
        result[np.isnan(result)] = noData
    if method is 'mode':        
        result,_ = stats.mode(dataTotal[:,:,:l], axis = 2)
    if method is 'sum':
        dataTotal[dataTotal==noData] = np.nan
        result = np.nansum(dataTotal[:,:,:l], axis = 2)
        result[np.isnan(result)] = noData
    
    print 'Maxima: ', np.nanmax(result)
    print 'Minima: ', np.nanmin(result)

    GIS.write(result, newFile, dataName, nodata = noData, dtype=dtype)

####################################################################################################
def firePeak(tile, start=startY+1, end=endY):
    ## Data Transfer
    print 'Calculating peak month...'
    src = os.environ['DATA'] + '/MODIS/MCD64A1/raw/' + tile + '/'
    dst = os.environ['DATA'] + '/MODIS/MCD64A1/annual/'
    
    flnm = dst + 'MCD64A1.peakMonth.' + tile + '.tif'
    if not os.path.exists(src):
        return None
    if not os.path.exists(dst):
        os.makedirs(dst)

    fireFreq = np.zeros((2400,2400,12), dtype=np.int)
    for dp in timerange(start, end, 'month'):
        flRaw = src + 'MCD64A1.A' + dp.strftime('%Y') + dp.strftime('%j') + '.' + tile + '.hdf'
        if os.path.isfile(flRaw):                        
            dataset = gdal.Open(flRaw) #open data                    
            subdataset = gdal.Open(dataset.GetSubDatasets()[0][0], gdal.GA_ReadOnly)
            fire = subdataset.ReadAsArray()
            fireFreq[:,:,dp.month-1] += (fire>0)
            print 'Peak month: ', tile, dp
        else:
            print 'Data MCD64A1 raw ', tile, dp, ' not found.'
    print 'Peak month: writing data for ', tile 
    season = np.argmax(fireFreq, axis=2)+1 
    season[~fireFreq.any(axis=2)] = gNoData
    GIS.write(season, flnm, dataset.GetSubDatasets()[0][0], dtype=gdal.GDT_Int16)

####################################################################################################
def burnFracTile(tile, winLen = 9):
    ## Count the burnt frequency on annual-resolution
    src = os.environ['DATA'] + '/MODIS/MCD64A1/raw/' + tile + '/'
    dst = os.environ['DATA'] + '/MODIS/MCD64A1/annual/'
    peakFileNm = dst + 'MCD64A1.peakMonth.' + tile + '.tif'
    #if os.path.isfile(dst+'MCD64A1.burnFracWin'+str(winLen).zfill(2)+'.'+tile+'.tif'):
    #    return

    peak = GIS.read(peakFileNm) 
    staMon = (peak-np.ceil(winLen/2.0))%12+1 #when winLen=12 peak=5, staMon = 12, endMon = 11
    endMon = (peak+np.floor(winLen/2.0)-1)%12+1
    del peak
    
    annuPre, annuAcu, yeaCoun, outCoun, rebCoun = \
        [np.zeros((2400, 2400), dtype=np.int) for _ in xrange(5)] #present year annual burnt
    #present year fire, accumulate annual fire, number of years, fire outside fire season, reburnt

    nodMask = np.ones((2400, 2400), dtype=np.bool) #nodata mask
    recMask = np.zeros((2400, 2400), dtype=np.bool) #recording mask
    outMask = np.zeros((2400, 2400), dtype=np.bool) #recording mask of outCount

    for dp in timerange(startY, endY, 'month'):         
        print 'burnFracTile, window length', str(winLen)+';', tile, dp
        flnm = src + 'MCD64A1.A' + dp.strftime('%Y') + dp.strftime('%j') + '.' + tile + '.hdf'
        if not os.path.isfile(flnm):
            continue
        recMask[staMon==dp.month] = True #recording starts
        outMask[staMon==dp.month] = True #recording starts, no turning off

        fire = GIS.read(flnm)
        annuPre[recMask] += fire[recMask]>0
        nodMask = nodMask&np.isnan(fire)
        
        outCoun[(~recMask)&outMask] += fire[(~recMask)&outMask]>0 #burnt outside of fire season
        endMask = (endMon==dp.month)&recMask #Mask to end and re-initialize.
        annuAcu[endMask&(annuPre>0)] += 1 #Burnt pixels
        rebCoun[endMask&(annuPre>1)] += 1 #re-Burnt pixels
        annuPre[endMask] = 0 #initialize present count to zero
        yeaCoun[endMask] += 1 #add one to year number
        recMask[endMask] = False #recording ends
    
    noBurnMask = yeaCoun==0
    yeaCoun[noBurnMask] = 1
    
    def aver(dataRaw):
        data = dataRaw*100.0/yeaCoun
        data[noBurnMask] = 0
        data[nodMask] = gNoData
        return data

    frac = aver(annuAcu)
    out = aver(outCoun)
    reb = aver(rebCoun)

    GIS.write(frac, dst+'MCD64A1.burnFracWin'+str(winLen).zfill(2)+'.'+tile+'.tif', peakFileNm)
    GIS.write(reb, dst+'MCD64A1.reBurnWin'+str(winLen).zfill(2)+'.'+tile+'.tif', peakFileNm)
    GIS.write(out, dst+'MCD64A1.outCounWin'+str(winLen).zfill(2)+'.'+tile+'.tif', peakFileNm)

####################################################################################################
def burnFrac(winLen = 9):
    [burnFracTile('h'+str(tileList[num][0]).zfill(2)+'v'+str(tileList[num][1]).zfill(2), winLen = winLen) for num in range(len(tileList))]
    GIS.mosa(os.environ['DATA'] + '/MODIS/MCD64A1/annual/MCD64A1.reBurnWin'+str(winLen).zfill(2)+'.', outPath+'MCD64A1.reBurnWin'+str(winLen).zfill(2)+'.tif')
    GIS.mosa(os.environ['DATA'] + '/MODIS/MCD64A1/annual/MCD64A1.outCounWin'+str(winLen).zfill(2)+'.', outPath+'MCD64A1.outCounWin'+str(winLen).zfill(2)+'.tif')
    GIS.mosa(os.environ['DATA'] + '/MODIS/MCD64A1/annual/MCD64A1.burnFracWin'+str(winLen).zfill(2)+'.', outPath+'MCD64A1.burnFracWin'+str(winLen).zfill(2)+'.tif')

####################################################################################################
def burnArea():
    ## Count the burnt frequency on annual-resolution
    dst = os.environ['DATA'] + '/MODIS/MCD64A1/annual/'
    afri = os.environ['DATA']+'/Africa/'
    #water = GIS.read(afri+'MCD12Q1.waterMask.tif')
    #water[water==1] = np.nan
    
    area = np.zeros((2400,2400), dtype=np.int)
    for year in range(startY+1, endY):
        for num in range(len(tileList)):
            tile = 'h'+str(tileList[num][0]).zfill(2)+'v'+str(tileList[num][1]).zfill(2)
            src = os.environ['DATA'] + '/MODIS/MCD64A1/raw/' + tile + '/'
            outfl = dst+'MCD64A1.burnArea'+str(year)+'.'+tile+'.tif'
            if os.path.isfile(outfl):
                continue

            annuBurn = np.zeros((2400,2400), dtype=np.int)
            nodMask = np.ones((2400,2400), dtype=np.bool)
            for dp in timerange(year, year+1, 'month'):         
                flnm = src + 'MCD64A1.A' + dp.strftime('%Y') + dp.strftime('%j') + '.' + tile + '.hdf'
                fire = GIS.read(flnm)
                annuBurn[fire>0] = 1
                nodMask = np.isnan(fire)&nodMask
            annuBurn[nodMask] = gNoData
            GIS.write(annuBurn, outfl, flnm)
    arrStac = np.dstack((GIS.mosa(dst+'MCD64A1.burnArea'+str(year)+'.', dst+'MCD64A1.burnArea'+str(year)+'.tif') for year in range(startY+1, endY)))
    aveArea = np.nanmean(arrStac, axis=2) 
    aveArea[np.isnan(aveArea)] = gNoData 
    GIS.write(aveArea, afri+'MCD64A1.burnAreaNatrYear.tif', afri+'MCD64A1.freq.tif')  

####################################################################################################
def fireMosa():
    src = os.environ['DATA'] + '/MODIS/MCD64A1/annual/'
    dst = os.environ['DATA'] + '/Africa/'

    [fireCount('h'+str(tileList[num][0]).zfill(2)+'v'+str(tileList[num][1]).zfill(2)) for num in range(len(tileList))]
    GIS.mosa(src+'MCD64A1.freq.', dst+'MCD64A1.freq.africa.tif')
    GIS.mosa(src+'MCD64A1.reBurn.', dst+'MCD64A1.reBurn.africa.tif')

    freq = GIS.read(dst+'MCD64A1.freq.africa.tif')
    freq[freq==0] = np.nan
    reIn = 1/freq
    reIn[np.isnan(reIn)] = 17

    GIS.write(reIn, 'temp.tif', dst+'MCD64A1.freq.africa.tif')
    GIS.resamp('temp.tif', dst+'MCD64A1.interval.africa.tif', method='near')

####################################################################################################
##########pre-Process MOD44B: merge two fill values
def treePrePro():
    for num in range(len(tileList)):
        tile = 'h'+ str(tileList[num][0]).zfill(2) + 'v' + str(tileList[num][1]).zfill(2)
        src = os.environ['DATA']+'/MODIS/MOD44B/raw/'+tile+'/'
        dst = os.environ['DATA']+'/MODIS/MOD44B/'+tile+'/'
        if not os.path.exists(dst):
            os.makedirs(dst)
        
        for date in timerange(startY, endY, 'year'):
            flnm = src + 'MOD44B.A' + date.strftime('%Y') + date.strftime('%j') + '.' + tile + '.hdf'
            if os.path.isfile(flnm):
                namestr = gdal.Open(flnm,gdal.GA_ReadOnly).GetSubDatasets()[0][0]
                outfile = dst + 'MOD44B.A' + date.strftime('%Y') + date.strftime('%j') + '.' + tile + '.hdf'
                os.system('gdal_calc.py -A %s --outfile=%s --calc="255*((A==200)|(A==253))+A*((A!=200)|(A!=253))" --NoDataValue=255 --overwrite' %(namestr, outfile))

####################################################################################################
def fire_16day(tile):
    from dateutil.relativedelta import relativedelta
    ## Data Transfer
    src = os.environ['DATA'] + '/MODIS/MCD64A1/raw/' + tile + '/'
    dst = os.environ['DATA'] + '/MODIS/MCD64A1/' + tile + '/'
    if not os.path.exists(src):
        return

    if not os.path.exists(dst):
        os.makedirs(dst)

    for dp in timerange(startY, endY, '16day'):
        f = 0 #flag of fire data existence
        flnm_n = dst + 'MCD64A1.A' + dp.strftime('%Y') + dp.strftime('%j') + '.' + tile + '.hdf' 
        if not os.path.isfile(flnm_n): #data already transformed
            d = [datetime.date(dp.year,dp.month,1), datetime.date(dp.year,dp.month,1) - relativedelta(months=1)]
            fire_p = np.empty((2400, 2400, 2), dtype=np.int) # 0: present; 1: previous
            for p in range(2):
                flnm = src + 'MCD64A1.A' + d[p].strftime('%Y') + d[p].strftime('%j') + '.' + tile + '.hdf'
                if os.path.isfile(flnm):                        
                    dataset = gdal.Open(flnm) #open data                    
                    subdataset = gdal.Open(dataset.GetSubDatasets()[0][0], gdal.GA_ReadOnly)
                    geo_transform = subdataset.GetGeoTransform() #geotransform 
                    x_size = subdataset.RasterXSize # Raster xsize
                    y_size = subdataset.RasterYSize # Raster ysize
                    srs = subdataset.GetProjectionRef() # Projection
                    
                    temp = subdataset.ReadAsArray() #open subdataset  
                    fire_p[:, :, p] = (datetime.date(d[p].year, 1, 1).toordinal() + temp - 1).astype(int)
                    fire_p[temp <= 0, p] = 0
                else:
                    f += 1
            if f != 2: #data not transformed
                fire_mask = np.logical_or((fire_p[...,0] <= dp.toordinal()) & (fire_p[...,0] > dp.toordinal()-16), fire_p[...,1] > dp.toordinal()-16)
                driver = gdal.GetDriverByName ('HDF4Image')
                dataset_out = driver.Create (flnm_n, x_size, y_size, 1, gdal.GDT_Byte)
                dataset_out.SetGeoTransform ( geo_transform )
                dataset_out.SetProjection ( srs )
                dataset_out.GetRasterBand ( 1 ).WriteArray ( fire_mask.astype(np.bool) )
                dataset_out.FlushCache()
                dataset_out = None
                print tile, dp

####################################################################################################
def trmmRot(): #fix rotation problem
    import string
    src = os.environ['DATA']+'/TRMM/'
    if not os.path.exists(dst):
        os.makedirs(dst)
    
    for month in range(1, 13):
        for year in range(2000,2014): #fix rotation problem
            flnm = src + '3B43.' + str(year) + str(month).zfill(2) + '.hdf'
            monPreci = src + '3B43.' + str(year) + str(month).zfill(2) + '.tif' #monthly precip
            monPreci_afr = dst + '3B43.' + str(month).zfill(2) + '.afr.tif'

            namestr = gdal.Open(flnm,gdal.GA_ReadOnly).GetSubDatasets()[0][0]        
            subdataset = gdal.Open(namestr, gdal.GA_ReadOnly) #open data                    
            geo_transform = subdataset.GetGeoTransform() #geotransform 
            y_size = subdataset.RasterXSize # Raster xsize
            x_size = subdataset.RasterYSize # Raster ysize
    
            driver = gdal.GetDriverByName ('GTiff')
            dataset_out = driver.Create (monPreci, x_size, y_size, 1, gdal.GDT_Float32)
            dataset_out.SetGeoTransform ((-180, 0.25, 0, 50, 0, -0.25)) #0,1,0,0,0,1. upper left x, pixel width, rotation,upper left y, rotation, pixel height
            
            outRasterSRS = osr.SpatialReference() #set reference to lon/lat system
            outRasterSRS.ImportFromEPSG(4326)
            dataset_out.SetProjection (outRasterSRS.ExportToWkt())
            
            dataset_out.GetRasterBand ( 1 ).WriteArray (np.rot90(subdataset.ReadAsArray()*24*365) ) #in unit of mm per year
            dataset_out.FlushCache()
            dataset_out = None

            GIS.resamp(monPreci, monPreci_afr) 

####################################################################################################
##########Process TRMM data
def trmmProc():
    import string
    src = os.environ['DATA']+'/TRMM/'
    dst = os.environ['DATA']+'/Africa/TRMM/'
    if not os.path.exists(dst):
        os.makedirs(dst)
    
    mMonPreci = src + '3B43.tif' #multi-year MAP
    mMonPreci_afr = dst + '3B43.africa.tif'
    
    monlist = []
    for month in range(1, 13):
        monPreci = src + '3B43.' + str(month).zfill(2) + '.tif' #mean annual precipitation
        monPreci_afr = dst + '3B43.' + str(month).zfill(2) + '.africa.tif' #mean annual precipitation 
        
        imgstr = ' '.join('-'+string.ascii_uppercase[year-2000]+' '+src+'3B43.'+str(year)+str(month).zfill(2)+'.tif' for year in range(2000,2014)) #list to string
        methstr = 'mean(array([' + ','.join(string.ascii_uppercase[l] for l in range(14)) + ']), axis=0)'
        
        monlist += ['-'+string.ascii_uppercase[month-1]+' '+monPreci] #list of all MAP for gdalcal
        os.system('gdal_calc.py %s --outfile=%s --calc="%s"' %(imgstr, monPreci, methstr)) #MAP of one year
        GIS.resamp(monPreci, monPreci_afr)
    imgstr = ' '.join(monlist)
    methstr = 'mean(array([' + ','.join(string.ascii_uppercase[l] for l in range(len(monlist))) + ']),axis=0)'
    os.system('gdal_calc.py %s --outfile=%s --calc="%s" --overwrite' %(imgstr, mMonPreci, methstr)) #multiyear MAP    
    GIS.resamp(mMonPreci, mMonPreci_afr)

####################################################################################################
def coverMask(thre=.1, maskNum=[0,16], northBound=None, southBound=None, westBound=None, eastBound=None):
    src = dataPath + 'MODIS/MCD12Q1/'
    dst = outPath+'MCD12Q1/tile/'
    #water mask
    for num in range(len(tileList)):
        tile = 'h'+str(tileList[num][0]).zfill(2)+'v'+str(tileList[num][1]).zfill(2)
        infl =  src + tile + '/MCD12Q1.A2001001.' + tile + '.hdf'
        outfl = dst + 'MCD12Q1.mask.' + tile+'.tif'
        
        cover = GIS.read(infl)
        coverMask = np.zeros(cover.shape, dtype=np.bool)
        for i in len(maskNum):
            coverMask = coverMask|(cover==maskNum[i])
        coverMask[np.isnan(cover)] = gNoData
        GIS.write(coverMask, outfl, infl, noData = gNoData, dtype=gdal.GDT_Int16)
    GIS.mosa(dst+'MCD12Q1.mask.', outPath + 'MCD12Q1.maskAve.tif')
    maskAve = GIS.read(outPath + 'MCD12Q1.maskAve.tif')
    mask = maskAve>thre


####################################################################################################
def maskCalcu():
    africaPath = dataPath +'/Africa/'
    
    #elev = GIS.resamp(dataPath + '/FAO/GloElev_30as.asc', africaPath+'/elev.africa.tif', resol = 'fine', sSrs = 'EPSG:4326')
    #elevMask = (elev>1200)|(elev<0)
    #elevMask = GIS.write(elevMask, africaPath+'elevMask.africa.tif', africaPath+'elev.africa.tif')
    
    cover = GIS.resamp(dataPath + '/GEM/GLC/glc2000_v1_1.tif', africaPath+'/GLC.africa.tif', resol = 'origin', sSrs = 'EPSG:4326')
    coverMask = cover==19
    #coverMask = (cover>=15)&(cover<=22)
    
    GIS.write(coverMask, africaPath+'GLCMask.orgin.tif', africaPath+'GLC.africa.tif')
    coverMask = GIS.resamp(africaPath+'GLCMask.orgin.tif', africaPath+'GLCMask.tif')

    totalMask = np.zeros(coverMask.shape, dtype=np.int)
    totalMask[coverMask>0.1] = 1
    #totalMask[elevMask>0.1] = 2
    
    totalMask = GIS.write(totalMask, africaPath + 'totalMask.tif', africaPath + 'GLCMask.tif', dtype=gdal.GDT_Int16)

####################################################################################################
def concealer()
from osgeo import ogr

# Create ring
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(1179091.1646903288, 712782.8838459781)
ring.AddPoint(1161053.0218226474, 667456.2684348812)
ring.AddPoint(1214704.933941905, 641092.8288590391)
ring.AddPoint(1228580.428455506, 682719.3123998424)
ring.AddPoint(1218405.0658121984, 721108.1805541387)
ring.AddPoint(1179091.1646903288, 712782.8838459781)

# Create polygon
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

print poly.ExportToWkt()
