from pyEarth import GIS, plot

####################################################################################################
def LandTypeScore(DataNameList = ['GLC2009','GLC2000'], TypeNameList = ['CropScore'], decompose=False):
    PathDict = {'GLC2009':dataPath + 'GEM/GLC/GLC2009/GLOBCOVER_L4_200901_200912_V2.3.tif','MCD12IGBP':dataPath+'/MODIS/MCD12Q1/MCD12IGBP2009.vrt', 'MCD12PFT':dataPath+'/MODIS/MCD12Q1/MCD12PFT2009.vrt','GLC2000':dataPath + 'GEM/GLC/GLC2000/glc2000_v1_1.tif'}
    
    def cover2score(FilePath, DataName, TypeNameList = TypeNameList, decompose=False):
        import pandas as pd
        origin = FilePath
        if not os.path.isfile(origin):
            GIS.resamp(FilePath, origin, method='near',resol='origin')
        if not os.path.isfile(modeCover):
            GIS.resamp(FilePath, modeCover, method='mode')

        raw = GIS.read(origin)
        ref = pd.read_excel(outPath+'LandTypeScore/LandScoreRule.xls', DataName, index_col=None, na_values=['NaN'])
        if decompose:
            TypeNameList = list(ref['ShortName']) #decompose to seperate PFTs
        for TypeName in TypeNameList:
            if TypeName=='Ignore':
                continue

            RAWName = outPath+'LandTypeScore/'+DataName+TypeName+'RAW.tif'
            outName = outPath+'LandTypeScore/'+DataName+TypeName+'.tif'
            if True:#not os.path.isfile(RAWName):
                ScoreMat = np.zeros(raw.shape)
                if decompose:
                    ScoreTable = ref[['Label','Value']]
                    ScoreTable[TypeName] = ref['ShortName'].map(lambda x: 100 if x==TypeName else np.nan if x=='Ignore' else 0)
                else:
                    ScoreTable = ref[['Label','Value',TypeName]]
                for ValSco in ScoreTable.iterrows():
                    score = ValSco[1][TypeName]
                    if score!=0:
                        print 'cover2score:', DataName+',', TypeName+'.', ValSco[1]['Label']+':', ValSco[1][TypeName]
                        ScoreMat[raw==ValSco[1]['Value']]=ValSco[1][TypeName]
                GIS.write(ScoreMat, RAWName, origin)
                del ScoreMat
            if not os.path.isfile(outName):
                 os.system('gdalwarp -ot Float32 -wt Float32 -overwrite -r average -tr 0.25 0.25 %s %s' %(inFile, outFile))
                
                GIS.resamp(RAWName, outName)
            plot.mapDraw(outName, outName+'.png',DataName+', '+TypeName, vMin=0, vMax=100)
    #####
    for DataName in DataNameList:
        print 'LandTypeScore:', DataName
        cover2score(PathDict[DataName], DataName, decompose=decompose)
