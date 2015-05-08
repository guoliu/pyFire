import os    
#'africa','asiaPaci','america','tropic_africa','tropic_asiaPaci','tropic_america','tropics','world'

def knife(source, region, target=None, method='average', draw=True, scale=1, vMin=0, vMax=100, titlTxt=None, areaWeight=False, bound = [-180,90,180,-90]):
    if target is None:
        target = source.rsplit('/',1)[-1].rsplit('.',1)[0]+'.'+region+'.tif'
    cutFile = 'vector/'+region+'.shp'
    if not os.path.isfile(target):        
        os.system('gdalwarp -ot Float32 -wt Float32 -dstnodata 255 -overwrite -t_srs WGS84 -cutline %s -crop_to_cutline -r %s -tr 0.25 0.25 %s %s' %(cutFile, method, source, target))
    if draw:
        from pyEarth import plot
        plot.mapDraw(target, scale=scale, vMin=vMin, vMax=vMax, titlTxt=titlTxt)
    if areaWeight:
        from pyEarth import GIS
        data = GIS.read(target)
        area = GIS.areaGrid(res=[0.25, 0.25], bound=bound)
        GIS.write(data*area/100.0*scale, target.rsplit('/',1)[-1].rsplit('.',1)[0]+'.weighted.tif', target)

####################################################################################################
def knifeSet(source, target=None):
    for region in ['tropic_africa','tropic_asiaPaci','tropic_america']:
        knife(source, region, target)

####################################################################################################
def axe(source, region, target=None, method='mode'):
    if target is None:
        target = source.rsplit('/',1)[-1].rsplit('.',1)[0]+'.'+region+'.vrt'
    cutFile = 'vector/'+region+'.shp' 
    os.system('gdalwarp -ot Float32 -wt Float32 -dstnodata 255 -overwrite -cutline %s -crop_to_cutline -r %s %s %s' %(cutFile, method, source, target))

####################################################################################################
def axeSet(source, target=None):
    for region in ['tropic_africa','tropic_asiaPaci','tropic_america']:
        axe(source, region, target)

####################################################################################################


'''
for name in ['wetSum.tif','wetMean.tif']:
    knifeSet('/data8/data/guol3/TRMM/3B42/'+name)


knife('/data8/data/guol3/MODIS/MOD44B/MOD44B.tif', 'tropics')

os.system('gdalwarp -ot Float32 -wt Float32 -srcnodata 255 -dstnodata 255 -overwrite -t_srs WGS84 -r %s -tr 0.25 0.25 %s %s' %('average', '/data8/data/guol3/MODIS/MOD44B/MOD44B.tif', 'MOD44B.tif'))
knife('/data8/data/guol3/TRMM/3B42/3B42_MAP.tif', 'tropics')
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2010.deciBroad.vrt','tropics')
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2010.greenBroad.vrt','tropics')
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2010.grass.vrt','tropics')
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2010.shrub.vrt','tropics')
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2010.broadCrop.vrt','tropics')
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2010.cerealCrop.vrt','tropics')

knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.irrigCrop.tif',region,scale=100,titlTxt='GLC: Post-flooding or irrigated croplands')
knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.mosaCrop.tif',region,scale=100,titlTxt='GLC: Mosaic cropland/vegetation')
knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.mosaVege.tif',region,scale=100,titlTxt='GLC: Mosaic vegetation/cropland')
knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.rainCrop.tif',region,scale=100,titlTxt='GLC: Rainfed croplands')





plot.mapDraw('MCD12Q1.2009.broadCrop.world.tif', scale=100, vMin=0, vMax=100, titlTxt='MCD12 PFT: Broadleaf Crops')
plot.mapDraw('MCD12Q1.2009.cerealCrop.world.tif', scale=100, vMin=0, vMax=100, titlTxt='MCD12 PFT: Cereal Crops')
plot.mapDraw('MCD12Q1.2009.totaCrop.world.tif', scale=100, vMin=0, vMax=100, titlTxt='MCD12 PFT: Total Crops')

region = 'world'

knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.irrigCrop.tif',region,scale=100,titlTxt='GLC: Post-flooding or irrigated croplands', areaWeight=True)
knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.mosaCrop.tif',region,scale=100,titlTxt='GLC: Mosaic cropland/vegetation', areaWeight=True)
knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.mosaVege.tif',region,scale=100,titlTxt='GLC: Mosaic vegetation/cropland', areaWeight=True)
knife('/data8/data/guol3/GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.rainCrop.tif',region,scale=100,titlTxt='GLC: Rainfed croplands', areaWeight=True)

os.system('gdal_calc.py -A GLOBCOVER_L4_200901_200912_V2.3.rainCrop.world.weighted.tif -B GLOBCOVER_L4_200901_200912_V2.3.mosaVege.world.weighted.tif -C GLOBCOVER_L4_200901_200912_V2.3.mosaCrop.world.weighted.tif -D GLOBCOVER_L4_200901_200912_V2.3.irrigCrop.world.weighted.tif --outfile=%s --calc="A+B*0.35+C*0.6+D"' %('GLOBCOVER_L4_200901_200912_V2.3.totaCrop.world.weighted.tif'))

region = 'TRMM'

knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2009.broadCrop.vrt',region,scale=100,titlTxt='MCD12 PFT: Broadleaf Crop', areaWeight=True, bound = [-180,50,180,-50])
knife('/data8/data/guol3/MODIS/MCD12Q1/MCD12Q1.2009.cerealCrop.vrt',region,scale=100,titlTxt='MCD12 PFT: Cereal Crop', areaWeight=True, bound = [-180,50,180,-50])

os.system('gdal_calc.py -A MCD12Q1.2009.cerealCrop.TRMM.weighted.tif -B MCD12Q1.2009.broadCrop.TRMM.weighted.tif --outfile=%s --calc="A+B"' %('MCD12Q1.2009.totaCrop.TRMM.weighted.tif'))
'''

for dataset in ['GLOBCOVER_L4_200901_200912_V2.3.rainCrop.world.weighted.tif','GLOBCOVER_L4_200901_200912_V2.3.totaCrop.world.weighted.tif','GLOBCOVER_L4_200901_200912_V2.3.mosaVege.world.weighted.tif','GLOBCOVER_L4_200901_200912_V2.3.mosaCrop.world.weighted.tif','GLOBCOVER_L4_200901_200912_V2.3.irrigCrop.world.weighted.tif','MCD12Q1.2009.cerealCrop.TRMM.tif','MCD12Q1.2009.totaCrop.TRMM.weighted.tif']:
    knife(dataset, 'africa', target=dataset.rsplit('.',3)[0]+'.weighted.africa.tif',  draw=True)
