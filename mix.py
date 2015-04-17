from conf import*
from . import GIS
from . import plot
from . import matr

####################################################################################################
def treeBoost(TreePath=outPath + 'MOD44B.tif', PreciPath=outPath + '3B43.tif', MaskPath=outPath + 'MCD12Q1.mask.tif', NamePre='MOD44'):
    print 'Boosting treecover, from', TreePath, 'to', PreciPath
    tree = GIS.read(TreePath)
    preci = GIS.read(PreciPath)
    fire = GIS.read(outPath + 'MCD64A1.burnFracWin07.tif')
    coverMask = GIS.read(MaskPath)>10

    treePur, preciPur, firePur = matr.cleaner([tree, preci, fire], mask=coverMask)
    
    preUpp, treeUpp = plot.scatter(preciPur, treePur, ['burt','unburnt'], NamePre+'PrecTree.png', divider=firePur==0, text=['Precipitation vs Tree Cover','Precipitation (mm/year)','Tree Cover (%)'])
    
    treeUpp[preUpp>2700] = treeUpp.max()
    treeUpp = np.insert(treeUpp, 0, 0)
    preUpp = np.insert(preUpp, 0, 0)
    
    treePo = GIS.reMap(preUpp, treeUpp, preci, template = outPath + '3B43.tif', flnm = outPath + NamePre + 'TreePo.tif')
     
    treeDif = treePo-tree
    treeDif[coverMask] = np.nan
    GIS.write(treeDif, outPath+NamePre+'TreeDif.tif', TreePath)

    treeDifRa = treeDif/treePo*100
    GIS.write(treeDifRa, outPath+NamePre+'TreeDifRa.tif', TreePath)

    plot.mapDraw(outPath+NamePre+'TreeDif.tif', NamePre+'TreeDif.png', 'Tree Cover Difference (%)', cmNm = 'gist_earth_r')
    plot.mapDraw(outPath+NamePre+'TreePo.tif', NamePre+'TreePo.png', 'Tree Cover Percentage (%)',cmNm = 'gist_earth_r')
    plot.mapDraw(outPath+NamePre+'TreeDifRa.tif', NamePre+'TreeDifRa.png', 'Tree Cover Difference Ratio (%)',cmNm = 'gist_earth_r')
