import os.path
import numpy as np
from pyEarth import GIS

startY=1998
endY=2014

####################################################################################################
def dailyAve():
    from nco import Nco
    import datetime
    nco = Nco()
    for d in range(365):
        dp =  datetime.date(startY,1,1)+datetime.timedelta(d)
        print "Averaging TRMM 3B42 for day "+dp.strftime('%j')+"..."
        ifile = ' '.join("3B42_daily."+str(year)+"."+dp.strftime('%m')+"."+dp.strftime('%d')+".7.nc" for year in range(startY,endY))
        ofile = "3B42_aver."+dp.strftime('%j')+".nc"           
        if not os.path.isfile(ofile):
            nco.ncra(input=ifile, output=ofile)
    nco.ncrcat(input="3B42_aver.*.nc", output="3B42_cat.nc", options="-d time,1,365")
    nco.ncwa(input="3B42_cat.nc", output="3B42_MAP.nc", options='-N -a time')

    return None

####################################################################################################
def TRMMwrite(data, fileNm):
    import gdal
    import osr
    y_size = 400 # Raster xsize
    x_size = 1440 # Raster ysize

    driver = gdal.GetDriverByName ('GTiff')
    dataset_out = driver.Create (fileNm, x_size, y_size, 1, gdal.GDT_Float32)
    dataset_out.SetGeoTransform ((-180, 0.25, 0, 50, 0, -0.25)) #0,1,0,0,0,1. upper left x, pixel width, rotation, upper left y, rotation, pixel height
    
    outRasterSRS = osr.SpatialReference() #set reference to lon/lat system
    outRasterSRS.ImportFromEPSG(4326)
    dataset_out.SetProjection (outRasterSRS.ExportToWkt())
    
    dataset_out.GetRasterBand(1).WriteArray(np.hstack((data[:,720:],data[:,:720]))) #in unit of mm per year; change from 0-360 to -180 to 180
    dataset_out.FlushCache()
    dataset_out = None

####################################################################################################
def SIcal(data):
    xSum = np.sum([data[:,:,day]*np.cos(day*2*np.pi/365) for day in range(365)],axis=0) #x dimensional sum
    ySum = np.sum([data[:,:,day]*np.sin(day*2*np.pi/365) for day in range(365)],axis=0) #y dimensional sum
    
    angle = np.arctan(ySum/xSum)+(xSum<0)*np.pi #x,y -> angle
    
    peak = np.round(angle*365/(2*np.pi)-0.5)%366 #retrieve angle -> day of year
    
    SI = np.sqrt(xSum**2+ySum**2)/np.sum(data,axis=2) #retrieve vector length -> seasonality index
    length = 365.0*(1-SI) #wet season length
    del SI
    staDay = np.round(peak-(length-1)/2.0).astype(np.int)%365
    endDay = np.round(peak+(length-1)/2.0).astype(np.int)%365
    
    seaSum = np.zeros(data.shape[:2]) #summary within season
    occSum = np.zeros(data.shape[:2]) #time of occurance
    sumMask = endDay<=staDay #mask of pixel in season
    for z in range(data.shape[2]):
        sumMask[z==staDay] = True
        sumMask[z==endDay] = False
        seaSum[sumMask]+=data[:,:,z][sumMask]
        occSum[sumMask]+=(data[:,:,z][sumMask]>0).astype(np.int)
    occRat = occSum*1.0/np.ceil(length)
    seaAve = seaSum/occSum
    
    TRMMwrite(occRat, 'wetRatio.tif')
    TRMMwrite(seaAve, 'wetMean.tif')
    TRMMwrite(peak, 'wetPeak.tif')
    TRMMwrite(length, 'wetLength.tif')
    TRMMwrite(seaSum, 'wetSum.tif')
    
    return None

####################################################################################################
def filter(data, weight=[0.1,0.1,0.1,0.1,0.2,0.1,0.1,0.1,0.1]):
    return np.sum(np.roll(data,i-4,axis=2)*weight[i] for i in range(9))

####################################################################################################

#dailyAve()
#data=np.dstack(GIS.read('3B42_cat.nc',num) for num in range(1,366))
SIcal(filter(data))
#data=GIS.read('3B42_MAP.nc')
#TRMMwrite(data,'3B42_MAP.tif')

from pyEarth import plot
plot.mapDraw('wetPeak.tif', 'Maxiumum Rainfall Day', cmNm='hsv',vMin=1, vMax=365, lut=366)
plot.mapDraw('wetSum.tif', r'$P_w$ (mm), Mean Wet season Precipitation', cmNm='rainbow_r',vMin=0, vMax=2000)
plot.mapDraw('wetMean.tif', r'$\alpha_w$ (mm), Mean Wet season Rainfall Depth', cmNm='rainbow_r',  vMin=0, vMax=18)
plot.mapDraw('wetLength.tif', 'Rainfall Season Length', cmNm='rainbow_r')
