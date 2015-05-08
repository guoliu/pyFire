from pyEarth import plot
import os
'''
plot.mapDraw('MCD12Q1.2009.cerealCrop.tropics.tif', scale=100, vMin=0, vMax=100, titlTxt='MCD12 PFT: Cereal Crop')
plot.mapDraw('MCD12Q1.2009.broadCrop.tropics.tif', scale=100, vMin=0, vMax=100, titlTxt='MCD12 PFT: Broad-leaf Crops')

os.system('gdal_calc.py -A GLOBCOVER_L4_200901_200912_V2.3.rainCrop.tropics.tif -B GLOBCOVER_L4_200901_200912_V2.3.mosaVege.tropics.tif -C GLOBCOVER_L4_200901_200912_V2.3.mosaCrop.tropics.tif -D GLOBCOVER_L4_200901_200912_V2.3.irrigCrop.tropics.tif --outfile=%s --calc="A+B*0.35+C*0.6+D"' %('GLOBCOVER_L4_200901_200912_V2.3.totaCrop.tropics.tif'))
'''
os.system('gdal_calc.py -A MCD12Q1.2009.cerealCrop.TRMM.tif -B MCD12Q1.2009.broadCrop.TRMM.tif --outfile=%s --calc="A+B"' %('MCD12Q1.2009.totaCrop.TRMM.tif'))

plot.mapDraw('MCD12Q1.2009.totaCrop.TRMM.tif', scale=100, vMin=0, vMax=100, titlTxt='MCD12 PFT: Total Crops')
plot.mapDraw('GLOBCOVER_L4_200901_200912_V2.3.totaCrop.world.tif', scale=100, vMin=0, vMax=100, titlTxt='GLC2009: Total Crops')
