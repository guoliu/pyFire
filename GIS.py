from conf import*

####################################################################################################
def read(flnm):
    g = gdal.Open(flnm, gdal.GA_ReadOnly)
    band = g.GetRasterBand(1)
    if band is None:
        subdataset = gdal.Open(g.GetSubDatasets()[0][0], gdal.GA_ReadOnly)
        band = subdataset.GetRasterBand(1)

    data = band.ReadAsArray()
    if band.GetNoDataValue():
        data = data.astype(np.float)
        data[data==band.GetNoDataValue()] = np.nan
    #if g.RasterCount>1:
    #    mask = g.GetRasterBand(2)
    #    maskData = mask.ReadAsArray()==255
    #    data[~maskData] = np.nan
    return data

####################################################################################################
def resamp(infile, outfile, method = 'average', noData = -9999, resol = 'coarse', sSrs = ''):
    if sSrs is not '':
        sSrs = ' -s_srs ' + sSrs

    if resol == 'coarse':
        os.system('gdalwarp -ot Float32 -wt Float32 -overwrite%s -t_srs $DATA/MODIS.prf -cutline $DATA/Africa/Africa.shp -crop_to_cutline -dstnodata %s -r %s -tr 27829.75 27829.75 %s %s' %(sSrs, str(noData), method, infile, outfile))
    elif resol == 'fine':
        os.system('gdalwarp -ot Float32 -wt Float32 -overwrite%s -t_srs $DATA/MODIS.prf -cutline $DATA/Africa/Africa.shp -crop_to_cutline -dstnodata %s -r %s -tr 5565.95 5565.95 %s %s' %(sSrs, str(noData), method, infile, outfile))
    elif resol == 'origin':
        os.system('gdalwarp -overwrite%s -t_srs $DATA/MODIS.prf -cutline $DATA/Africa/Africa.shp -crop_to_cutline -dstnodata %s -r %s %s %s' %(sSrs, str(noData), method, infile, outfile))
    else: 
        print 'Unsupported resulotion: ', resol
        return

    #return read(outfile)

####################################################################################################
def write(indata, outfile, template, noData = gNoData, drivernm = 'GTiff', dtype=gdal.GDT_Float32):
    print 'Writing data to ', outfile 
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
    dataset_out = driver.Create (outfile, x_size, y_size, 1, dtype)
    dataset_out.SetGeoTransform ( geo_transform )
    dataset_out.SetProjection ( srs )
    raster_out = dataset_out.GetRasterBand ( 1 )
    
    indata[np.isnan(indata)] = noData
    raster_out.WriteArray (indata)
    raster_out.SetNoDataValue(noData)
    dataset_out.FlushCache()
    dataset_out = None
    
    #resamp('temp.tif', outfile, method = method)
    #return read(outfile)

####################################################################################################
### combine and resample tiles to GeoTiff with a given product name and date
def mosa(oldPref, newFile, method='average', oldSuff='.tif'):
    nmlist = [oldPref+'h'+str(tileList[num][0]).zfill(2)+'v'+str(tileList[num][1]).zfill(2) +oldSuff for num in range(len(tileList))]
    tilestr = ' '.join(nmlist)
    
    resamp(tilestr, newFile, method=method)
    return read(newFile)

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
