import os, random, numpy as np, matplotlib as mpl, numpy.ma as ma
from osgeo import gdal, osr

veg_dict = {'water': 0,
          'greenNeedle': 1,
          'greenBroad': 2,
          'deciNeedle': 3,
          'deciBroad': 4,
          'mixed': 5,
          'closedShrub': 6,
          'openShrub': 7,
          'woodySavanna': 8,
          'savanna': 9,
          'grass': 10,
          'wetland': 11,
          'crop': 12,
          'urban': 13,
          'mosaic': 14,
          'cryo': 15,
          'barren': 16,
          'unclass': 17}

tilelist = [[16,6],[16,7],[16,8],[17,5],[17,6],[17,7],[17,8],[18,5],[18,6],[18,7],[18,8],[18,9],[19,5],[19,6],[19,7],[19,8],[19,9],[19,10],[19,11],[19,12],[20,5],[20,6],[20,7],[20,8],[20,9],[20,10],[20,11],[20,12],[21,6],[21,7],[21,8],[21,9],[21,10],[21,11],[22,7],[22,8],[22,9],[23,7],[23,8]] 
startY, endY = 2000, 2014
dataPath = os.environ['DATA']
gNoData = -9999
