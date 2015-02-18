from conf import*
from osgeo import ogr
from osgeo import osr

####################################################################################################
def read(flnm, num=1):
    print 'Reading', flnm
    g = gdal.Open(flnm, gdal.GA_ReadOnly)
    band = g.GetRasterBand(num)
    if band is None:
        subdataset = gdal.Open(g.GetSubDatasets()[num-1][0], gdal.GA_ReadOnly)
        band = subdataset.GetRasterBand(1)

    data = band.ReadAsArray()
    if band.GetNoDataValue():
        data = data.astype(np.float)
        data[data==band.GetNoDataValue()] = np.nan
    #if g.RasterCount>1:
    #    mask = g.GetRasterBand(2)
    #    maskData = mask.ReadAsArray()==255
    #    data[~maskData] = np.nan
    g = None
    return data

####################################################################################################
def resamp(inFile, outFile, region='africa', method = 'average', noData = gNoData, resol = 'coarse', sSrs = '', tSrs='$DATA/MODIS.prf'):
    """        
    Resample raster data, and crop to cutline. 

    Args:
        inFile: Input raster file name.
        outFile: Output raster file name.
        method: average (default), near, bilinear, cubic, cubicspline, lanczos, mode.
    """
    #EPSG:4326
    print 'Resampling', inFile
    if tSrs!='$DATA/MODIS.prf' and resol!='origin':
        print 'Only original resolution is supported for non-Sinusoidal projection. Falling back to "origin".'
        resol='origin'

    if sSrs!='':
        sSrs='-s_srs '+sSrs
    if tSrs!='':
        tSrs='-t_srs '+tSrs

    if resol=='coarse': 
        tRes='-tr 25630.65 27685.95'
    elif resol == 'fine':
        #tRes='-tr 5565.95 5565.95'
        tRes='-tr 2782.975 2782.9'
    elif resol == 'origin':
        tRes=''
    else: 
        print 'Unsupported resulotion:', resol
        return
    
    cutFile = '$DATA/vector/'+region+'.shp'

    #print 'gdalwarp -ot Float32 -wt Float32 -overwrite %s %s -cutline $DATA/Africa/Africa.shp -crop_to_cutline -dstnodata %s -r %s %s %s %s' %(sSrs, tSrs, str(noData), method, tRes, inFile, outFile)
    os.system('gdalwarp -ot Float32 -wt Float32 -overwrite %s %s -cutline %s -crop_to_cutline -dstnodata %s -r %s %s %s %s' %(sSrs, tSrs, cutFile, str(noData), method, tRes, inFile, outFile))

    return

####################################################################################################
def write(indata, outFile, template, noData = gNoData, drivernm = 'GTiff', dtype=gdal.GDT_Float32):
    print 'Writing data to ', outFile 
    G = gdal.Open(template, gdal.GA_ReadOnly) #open data                    
    try:
        g = gdal.Open(G.GetSubDatasets()[0][0], gdal.GA_ReadOnly)
    except:
        g = G 
    geo_transform = g.GetGeoTransform() #geotransform 
    x_size = g.RasterXSize # Raster xsize
    y_size = g.RasterYSize # Raster ysize
    srs = g.GetProjectionRef() # Projection

    driver = gdal.GetDriverByName (drivernm)
    dataset_out = driver.Create (outFile, x_size, y_size, 1, dtype)
    dataset_out.SetGeoTransform ( geo_transform )
    dataset_out.SetProjection ( srs )
    raster_out = dataset_out.GetRasterBand ( 1 )
    
    nanMask = np.isnan(indata)
    indata[nanMask] = noData
    raster_out.WriteArray (indata)
    raster_out.SetNoDataValue(noData)
    dataset_out.FlushCache()
    dataset_out = None
    
    try:
        indata[nanMask] = np.nan
    except: 
        pass

    return
    #resamp('temp.tif', outFile, method = method)
    #return read(outFile)

####################################################################################################
### combine and resample tiles to GeoTiff with a given product name and date
def mosa(oldPref, newFile, method='average', oldSuff='.tif'):
    nmlist = [oldPref+'h'+str(tileList[num][0]).zfill(2)+'v'+str(tileList[num][1]).zfill(2) +oldSuff for num in range(len(tileList))]
    tilestr = ' '.join(nmlist)
    
    resamp(tilestr, newFile, method=method)
    return
 
####################################################################################################
def create(extent, outShapefile): #(left, right, bottom, top)
    # Create a Polygon from the extent tuple
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(extent[0],extent[2])
    ring.AddPoint(extent[1], extent[2])
    ring.AddPoint(extent[1], extent[3])
    ring.AddPoint(extent[0], extent[3])
    ring.AddPoint(extent[0],extent[2])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    # Save extent to a new Shapefile
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
    
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer("extent", geom_type=ogr.wkbPolygon)
    
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    outLayer.CreateField(idField)
    
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(poly)
    feature.SetField("id", 1)
    outLayer.CreateFeature(feature)
    
    # Close DataSource
    outDataSource.Destroy

####################################################################################################
def cordConv(xy_source, inproj, outproj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy

####################################################################################################
def area(longiDeg, latiDeg, resoDeg):
    R = 6371.0087714 #km. Arithmetic mean radius of Earth in WGS84
    
    if len(resoDeg)==1:
        resoDeg = [resoDeg,resoDeg]

    reOrg = lambda x, lim: np.array([x[0],x[1]+lim]) if x[1]<x[0] else np.array(x)
    longiDeg = reOrg(longiDeg, 180) #circula longitude
    latiDeg = reOrg(latiDeg, 0) #normal latitude

    tran = lambda x: x.astype(np.float64)*np.pi/180.0 #degree to radius    
    longi = tran(longiDeg) #longitude extend in radius
    lati = tran(latiDeg) #latitude extend in radius
    reso = tran(np.array(resoDeg))

    d = lambda x, r: (np.arange(x[0],x[1]-r*1e-2,r), np.append(np.arange(x[0]+r,x[1],r),x[1]))

    xLefVe, xRitVe = (x%180 for x in d(longi, reso[0])) #x vectors; round earth
    yButVe, yTopVe = d(lati, reso[1]) #y vectors

    xLef, yBut = np.meshgrid(xLefVe, yButVe) #buttom-left grids
    xRit, yTop = np.meshgrid(xRitVe, yTopVe) #top-right grids

    A = R**2*(xRit-xLef)*(np.sin(yTop)-np.sin(yBut))
    return A

####################################################################################################
def areaDiv(inFile, outFile):
    ds = gdal.Open(inFile)
    gt = ds.GetGeoTransform()

    xres = gt[1]
    yres = gt[5]

    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0]
    xmax = gt[0] + (xres * ds.RasterXSize)
    ymin = gt[3] + (yres * ds.RasterYSize)
    ymax = gt[3]

    data=read(inFile)
    gridArea=area([xmin, xmax], [ymin, ymax], [xres,-yres])
    
    write(gridArea, outPath+'gridArea.tif', inFile)
    write(data/gridArea, outFile, inFile)

####################################################################################################
def reMap(xUpp, yUpp, xBas, template = None, flnm = None):
    xMask = ~np.isnan(xBas)
    
    yPo = np.empty(xBas.shape)
    yPo[:] = np.nan
    yPo[xMask] = np.interp(xBas[xMask], xUpp, yUpp)
    
    if template is not None and flnm is not None:
        write(yPo, flnm, template)

    return yPo
