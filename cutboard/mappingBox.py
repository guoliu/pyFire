from pyEarth import plot

plot.mapDraw('MOD44B.tropics.boost.tif', 'Potential Tree Cover (%)', vMin=0, vMax=80, cmNm='terrain_r')
plot.mapDraw('MOD44B.tropics.tif', 'MOD44 Tree Cover (%)', vMin=0, vMax=80, cmNm='terrain_r')

plot.mapDraw('/data8/data/guol3/TRMM/3B42/wetPeak.tif', 'Maxiumum Rainfall Day', cmNm='hsv',vMin=1, vMax=365, lut=366)
plot.mapDraw('/data8/data/guol3/TRMM/3B42/wetLength.tif', 'Rainfall Season Length', cmNm='rainbow_r')
plot.mapDraw('/data8/data/guol3/TRMM/3B42/3B42_MAP.tif', 'Mean Annual Precipitation', cmNm='rainbow_r', vMin=0, vMax=4000)
plot.mapDraw('wetSum.tropics.tif', r'$P_w$ (mm), Mean Wet season Precipitation', cmNm='rainbow_r',vMin=0, vMax=2000)
plot.mapDraw('wetMean.tropics.tif', r'$\alpha_w$ (mm), Mean Wet season Rainfall Depth', cmNm='rainbow_r',  vMin=0, vMax=18)

plot.mapDraw('/data8/data/guol3/TRMM/3B42/wetSum.tif', r'$P_w$ (mm), Mean Wet season Precipitation', cmNm='rainbow_r',vMin=0, vMax=2000)
plot.mapDraw('/data8/data/guol3/TRMM/3B42/wetMean.tif', r'$\alpha_w$ (mm), Mean Wet season Rainfall Depth', cmNm='rainbow_r',  vMin=0, vMax=18)
