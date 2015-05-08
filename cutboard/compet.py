from pyEarth import matr, GIS, plot
import numpy as np

####################################################################################################
def uppBound(xRaw, yRawList, nmList, figNm, text, nbin=None, percen=95, alpha=0.4, xlim=None):    
    """
    Draw scatter plot in groups, and calculated upper percentile and median. Accept and ignore NaNs.

    Args:
        xRaw
    
    """
    
    titText, xText, yText = text

    almost_black = '#262626'
    import matplotlib.pyplot as plt, brewer2mpl
    clm = brewer2mpl.get_map('Set3', 'qualitative', 8, reverse=True).mpl_colors

    fig, ax = plt.subplots(1)
    for i in range(len(yRawList)):
        print 'Plotting', nmList[i]
        x = xRaw
        y = yRawList[i]
        color = clm[i%8]
        ax.scatter(x, y, alpha=alpha, s=10, facecolors='none', edgecolor=color, linewidth=0.5)
        xNew, yUpp = matr.binPer(xRaw, yRawList[i], nbin = nbin, percen = percen, med = False)
        ax.plot(xNew, yUpp, color=color, linewidth=1.4, label=nmList[i], alpha=1)
        ax.plot(xNew, yUpp, color='black', linewidth=0.2, alpha=1)
        alpha = alpha*0.8
    
    legend = ax.legend(frameon=True, fontsize=10, loc='best')#, bbox_to_anchor=(1, 0.5))
    rect = legend.get_frame()
    rect.set_alpha(0.7)
    rect.set_linewidth(0.0)
    
    ax.set_title(titText)
    ax.set_xlabel(xText, fontsize=14)
    ax.set_ylabel(yText, fontsize=14)
    if xlim is None:
        ax.set_xlim([np.nanmin(xRaw)*0.85,np.nanmax(xRaw)*1.1])
    else:
        ax.set_xlim(xlim)
    ax.set_ylim([0,110])  
    fig.savefig(figNm, dpi=300)
    plt.close()

####################################################################################################

treeRaw = GIS.read('MOD44B.tropics.tif')
greenBroadRaw = GIS.read('MCD12Q1.2010.greenBroad.tropics.tif')*100
deciBroadRaw = GIS.read('MCD12Q1.2010.deciBroad.tropics.tif')*100
grassRaw = GIS.read('MCD12Q1.2010.grass.tropics.tif')*100
shrubRaw = GIS.read('MCD12Q1.2010.shrub.tropics.tif')*100
preciRaw = GIS.read('3B42_MAP.tropics.tif')
PwRaw = GIS.read('wetSum.tropics.tif')
AwRaw = GIS.read('wetMean.tropics.tif')

'''
####################################################################################################
#####Compite under precipitation change
tree, greenBroad, deciBroad, grass, shrub, preci = matr.cleaner([treeRaw,greenBroadRaw,deciBroadRaw,grassRaw,shrubRaw,preciRaw])

uppBound(preci, [tree, greenBroad, deciBroad, shrub, grass], ['TreeCover','EvergreenBroadLeaf','DecidiousBroadLeaf','shrub','Grass'], 'preci_compet.png',['PFTs Competition','Precipitation (MM/year)','Percentage'], alpha=0.1,xlim=[0,6000])

####################################################################################################
#####Compite under 2D constrain

for a in [1, 3, 5, 7, 9, 11, 13, 15]:
    mask=(AwRaw>a-1.0)&(AwRaw<a+1.0)
    greenBroad, deciBroad, grass, shrub, Pw = matr.cleaner([greenBroadRaw,deciBroadRaw,grassRaw,shrubRaw,PwRaw], mask=(~mask))
    uppBound(Pw, [greenBroad, deciBroad, shrub, grass], ['EvergreenBroadLeaf','DecidiousBroadLeaf','shrub','Grass'], '2D_compet_a'+str(a)+'.png',['PFTs Competition, '+r'$\alpha_w$='+str(a)+'mm','Mean Wet Season Precipitation (mm)','Percentage'],alpha=1,xlim=[0,2000])

for p in [200,600,1000,1400,1800]:
    mask=(PwRaw>p-200)&(PwRaw<p+200)
    greenBroad, deciBroad, grass, shrub, Aw = matr.cleaner([greenBroadRaw,deciBroadRaw,grassRaw,shrubRaw,AwRaw], mask=(~mask))
    uppBound(Aw, [greenBroad, deciBroad, shrub, grass], ['EvergreenBroadLeaf','DecidiousBroadLeaf','shrub','Grass'], '2D_compet_p'+str(p)+'.png',['PFTs Competition, '+r'$P_w$='+str(p)+'mm','Mean Wet Season Rainfall Depth (mm)','Percentage'],alpha=1,xlim=[0,20])
'''
####################################################################################################
#####shrub vs grass
for p in [200,400,800]:
    mask=(PwRaw>p-200)&(PwRaw<p+200)
    grass, shrub, Aw = matr.cleaner([grassRaw,shrubRaw,AwRaw], mask=(~mask))
    uppBound(Aw, [shrub, grass], ['shrub','Grass'], '2D_ShrubGrassCompet_p'+str(p)+'.png',['Grass-Shrub Competition, '+r'$P_w$='+str(p)+'mm','Mean Wet Season Rainfall Depth (mm)','Percentage'],alpha=1,xlim=[0,8])

####################################################################################################
#####greenBroad vs deciBroad
for a in [7, 9, 11, 13]:
    mask=(AwRaw>a-1.0)&(AwRaw<a+1.0)
    greenBroad, deciBroad, Pw = matr.cleaner([greenBroadRaw,deciBroadRaw,PwRaw], mask=(~mask))
    uppBound(Pw, [greenBroad, deciBroad], ['EvergreenBroadLeaf','DecidiousBroadLeaf'], '2D_GreenDecidCompet_a'+str(a)+'.png',['Evergreen-Decidious Competition, '+r'$\alpha_w$='+str(a)+'mm','Mean Wet Season Precipitation (mm)','Percentage'],alpha=1,xlim=[800,1800])


