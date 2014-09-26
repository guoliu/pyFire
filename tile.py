import os, glob, random, datetime, sys, math, socket, collections, numpy as np, numpy.ma as ma,  matplotlib as mpl, pandas as pd
from osgeo import gdal
from dateutil.relativedelta import relativedelta
from datetime import date
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

qua_dic = {0: 'Good',
            1: 'Marginal',
            2: 'Snow/Ice',
            3: 'Cloudy'}
    
samp_dic = {0: 'Constant',
            1: 'Variable'}

startY = 2000
endY = 2014
annuN = 365//16+1
preLen_min = np.floor(1.5*annuN)

workPath = '/gdata/randerson/guol3/'

fireRaw = '/export/gdata/randerson2/giglio/MCD64A1/'
firePath = workPath + 'MCD64A1/'

EVIRaw = '/gdata/randerson2/group/MODIS/MOD13A1-by-tile/'
EVIPath = workPath + 'MOD13A1/'

treeRaw = '/gdata/randerson3/group/MODIS/MOD44B-by-tile/V005/'
treePath = workPath + 'MOD44B/'

coverRaw = '/gdata/randerson2/group/MODIS/MCD12Q1/V051/'
coverPath = workPath + 'MCD12Q1/'

if socket.gethostname() == 'Lycopodium':
    localPath = '/Volumes/LycopodiumHDrive/Research/Data/'
    firePath2 = localPath + 'MCD64A1/'
    EVIPath2 = localPath + 'MOD13A1/'
    coverPath2 = localPath + 'MCD12Q1/'
    os.system('scp guol3@gplogin1.ps.uci.edu:' + firePath + tile + '/* ' + firePath2 + tile)
    os.system('scp guol3@gplogin1.ps.uci.edu:' + EVIPath + tile + '/* ' + EVIPath2 + tile)
    os.system('scp guol3@gplogin1.ps.uci.edu:' + coverPath + tile + '/* ' + coverPath2 + tile)
    firePath = firePath2
    EVIPath = EVIPath2
    coverPath = coverPath2

####################################################################################################
def daterange16(start, end):
    return (datetime.date(y,1,1)+datetime.timedelta(d-1) for y in range(start, end) for d in range(1,366,16))

####################################################################################################
def daterangem(start, end):
    return (datetime.date(y, m, 1) for y in range(start, end) for m in range(1, 13))

####################################################################################################
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return np.nanmean(a.reshape(sh),axis=-1).mean(1)

####################################################################################################
def tileCord(tile):
    cPath = os.path.dirname(os.path.realpath(__file__)) 
    tile_co_list = np.loadtxt(cPath+'/tilelist.txt')

    lonDic = {-1: 'W',
               0: '',
               1: 'E'}

    latDic = {-1: 'S',
               0: '',
               1: 'N'}

    lonN = tile_co_list[(tile_co_list[:,0]==np.int(tile[1:3]))&(tile_co_list[:,1]==np.int(tile[4:6])),2][0]
    lon = np.str(np.absolute(lonN))+lonDic[np.sign(lonN)]
    
    latN = tile_co_list[(tile_co_list[:,0]==np.int(tile[1:3]))&(tile_co_list[:,1]==np.int(tile[4:6])),3][0]
    lat = np.str(np.absolute(latN))+latDic[np.sign(latN)]

    return lat+','+lon
####################################################################################################
def fire_read(tile):
## Data Transfer
    src = fireRaw + tile + '/'
    dst = firePath + tile  + '/'
    if not os.path.exists(dst):
        print 'Fire data is not ready. Trying to prepare data now...'
        os.makedirs(dst)
    open(tile + 'firemissing.txt', 'w').close()

    fire = np.ones((2400, 2400, (endY - startY) * annuN), dtype=np.bool)
    l = 0
    for dp in daterange16(startY, endY):
        f = 0 #flag of fire data existence
        flnm_n = dst + 'MCD64A1.' + dp.strftime('%Y') + dp.strftime('%j') + tile + '.npy' 
        if os.path.isfile(flnm_n): #data already transformed 
            fire[..., l] = np.load(flnm_n)
        else: #data not transformed or missing
            d = [datetime.date(dp.year,dp.month,1), datetime.date(dp.year,dp.month,1) - relativedelta(months=1)]
            fire_p = np.empty((2400, 2400, 2), dtype=np.int) # 0: present; 1: previous
            for p in range(2):
                flnm = glob.glob(src + 'MCD64A1.A' + d[p].strftime('%Y') + d[p].strftime('%j') + '.' + tile + '.051.*.hdf')
                if flnm:                        
                    temp = gdal.Open('HDF4_EOS:EOS_GRID:"' + flnm[0] + '":MOD_Grid_Monthly_500m_DB_BA:Burn Date', gdal.GA_ReadOnly).ReadAsArray()
                    fire_p[:, :, p] = (datetime.date(d[p].year, 1, 1).toordinal() + temp - 1).astype(int)
                    fire_p[temp <= 0, p] = 0
                else:
                    f += 1
            if f != 2: #data not transformed
                fire_mask = np.logical_or((fire_p[...,0] <= dp.toordinal()) & (fire_p[...,0] > dp.toordinal()-16), fire_p[...,1] > dp.toordinal()-16)
                fire[..., l] = fire_mask
                np.save(flnm_n, fire_mask)
                print dp
            else: #data missing
                myfile = open(tile + 'firemissing.txt', 'a')
                myfile.write('No fire data for ' + dp.strftime('%Y') + ',' + dp.strftime('%j') + ';\n')
                myfile.close()
        l += 1
    return fire

####################################################################################################
#0   Good Data          Use with confidence
#1   Marginal data      Useful, but look at other QA information
#2   Snow/Ice           Target covered with snow/ice
#3   Cloudy             Target not visible, covered with cloud
def EVI_read(tile, date, qua = 0, forced = False):     
    src = EVIRaw + tile + '/' #2000049-2013305
    dst = EVIPath + tile + '/'
    flnmS = src + 'MOD13A1.A' + date.strftime('%Y') + date.strftime('%j') + '.' + tile + '.hdf'
    flnmE = dst + 'MOD13A1.' + date.strftime('%Y') + date.strftime('%j') + tile + 'EVI_qua' + str(qua) + '.npy'
    flnmN = dst + 'MOD13A1.' + date.strftime('%Y') + date.strftime('%j') + tile + 'NBR_qua' + str(qua) + '.npy'  
    
    if not forced:
        if os.path.isfile(flnmE)&os.path.isfile(flnmN):
            EVI = np.load(flnmE)
            NBR = np.load(flnmN)
            return EVI, NBR

    if not os.path.exists(dst):
        os.makedirs(dst)
        print 'Transforming EVI', tile

    if os.path.isfile(flnmS):
        relia = gdal.Open('HDF4_EOS:EOS_GRID:"' + flnmS + '":MODIS_Grid_16DAY_500m_VI:500m 16 days pixel reliability', gdal.GA_ReadOnly).ReadAsArray().astype(np.int)
        reliaM = np.logical_or(relia==-1,relia>qua)
        EVI = gdal.Open('HDF4_EOS:EOS_GRID:"' + flnmS + '":MODIS_Grid_16DAY_500m_VI:500m 16 days EVI', gdal.GA_ReadOnly).ReadAsArray().astype(np.float)*.0001
        EVI[EVI<-.2] = np.nan 
        EVI[reliaM] = np.nan
        np.save(flnmE, EVI)
        temp = np.dstack((gdal.Open('HDF4_EOS:EOS_GRID:"' + flnmS + '":MODIS_Grid_16DAY_500m_VI:500m 16 days NIR reflectance', gdal.GA_ReadOnly).ReadAsArray(), gdal.Open('HDF4_EOS:EOS_GRID:"' + flnmS + '":MODIS_Grid_16DAY_500m_VI:500m 16 days MIR reflectance', gdal.GA_ReadOnly).ReadAsArray())).astype(np.float)*.0001
        temp[temp<0] = np.nan
        temp[reliaM,:] = np.nan
        temp_sum = np.sum(temp, axis=2)
        temp_sum[temp_sum==0] = np.nan
        NBR = (temp[..., 0] - temp[..., 1]) / temp_sum
        np.save(flnmN, NBR)
        return EVI, NBR    
    else:
        myfile = open(tile + 'EVImissing.txt', 'a')
        myfile.write(date.strftime('%Y') + date.strftime('%j') + tile + '\n')
        myfile.close()
        EVI = np.empty([2400, 2400])
        NBR = np.empty([2400, 2400])
        EVI[:] = np.nan
        NBR[:] = np.nan
        return EVI, NBR

####################################################################################################
def tree_read(tile, date = datetime.date(2001,1,1)):
    date = datetime.date(date.year,1,1)
    src = treeRaw + tile + '/'
    flnm = src + 'MOD44B.A' + date.strftime('%Y') + '065.' + tile + '.005.hdf'
    
    dst = treePath + tile + '/'
    flnm_n = dst + 'MOD44B.' +  date.strftime('%Y') + tile + '.npy' 
    if os.path.isfile(flnm_n):
        return np.load(flnm_n)
    
    if not os.path.exists(dst):
        os.makedirs(dst)
        print 'Transferring treecover data'
    if os.path.isfile(flnm):
        temp = gdal.Open('HDF4_EOS:EOS_GRID:"' + flnm + '":MOD44B_250m_GRID:Percent_Tree_Cover', gdal.GA_ReadOnly).ReadAsArray().astype(np.float)  
        temp[temp>100] = np.nan
        tree = rebin(temp, [2400,2400])
        #np.save(flnm_n, tree)
    else:
        myfile = open(tile + 'treemissing.txt', 'a')
        myfile.write(flnm + '; occurred for 0 times.\n')
        myfile.close
        tree = np.empty([2400,2400])
        tree[:] = np.nan
    return tree 

####################################################################################################
def cover_read(tile, date = datetime.date(2001,1,1), forced = False, detect = False):
    date = datetime.date(date.year,1,1)
    src = coverRaw + date.strftime('%Y') + '.01.01/'
    if not os.path.exists(src): #date out of bound
        return False
    
    dst = coverPath + tile + '/'
    flnm_n = dst + 'MCD12Q1.' + date.strftime('%Y') + tile + '.npy' #file name at destination
    
    if not forced and not detect: #if not foreced updated then return previous stored data
        if os.path.isfile(flnm_n):
            return np.load(flnm_n)    
    
    flnm = src + 'MCD12Q1.A' + date.strftime('%Y') + '001.' + tile + '.051.hdf'
    if os.path.isfile(flnm):
        if detect:
            return True
        cover = gdal.Open('HDF4_EOS:EOS_GRID:"' + flnm + '":MOD12Q1:Land_Cover_Type_1', gdal.GA_ReadOnly).ReadAsArray().astype(np.int)
        cover[cover>16] = len(veg_dict) #unclassified as 17        
        if not os.path.exists(dst): #new directory
            os.makedirs(dst)
        np.save(flnm_n, cover)
        return cover

    if cover_read(tile, date = date + relativedelta(years=1), detect = detect) is not False: #tile exis in next year
        myfile = open(tile + 'covermissing.txt', 'a')
        myfile.write(flnm + '\n')
        myfile.close
    return cover_read(tile, date = date + relativedelta(years=1), detect = detect)

####################################################################################################
def dTS(NBR, startM = 14):
    step = 23
    aNBR = NBR[..., range(startM, startM + step * 14, step)]
    butt = np.empty((ind[1] - ind[0], ind[3] - ind[2]))
    butt[:] = np.nan
    dNBR = aNBR - np.dstack((butt, aNBR[..., :-1]))
    x_v = []
    for dp in daterange16(2000, 2014):
        x_v = np.append(x_v, dp.toordinal())

    x_d = x_v[range(startM, step * 14, step)]
    return (dNBR, x_d)

####################################################################################################
def fireSamp(tile, veg_ty = 'savanna', fire=None):
    if os.path.isfile(tile + veg_ty+ 'Samp.npz'):
        with np.load(tile + veg_ty+ 'Samp.npz') as data:
            subfire = data['subfire']
            subEVI = data['subEVI']
            subNBR = data['subNBR']
            s = data['s']
        return subfire, subEVI, subNBR, s
    
    veg_num = veg_dict[veg_ty]
    direF = firePath + tile + '/'
    direV = EVIPath + tile + '/'
    if not os.path.exists(direV):
        print 'Vegetation data is not ready. Tring to prepare data now...'
        EVI_read(tile)
    if fire is None:
        fire = fire_read(tile)
    
    s = 0 # Start date number of data. Assuming fire data starts later than vegetation data.
    while fire[...,s].all(): 
        s += 1

    fireM = np.any(fire[...,s:], axis=2)
    cover = cover_read(tile)

    cover_m = cover == veg_num
    pool = np.array(np.nonzero(np.logical_and(fireM, cover_m))).T
    if np.size(pool,0) < 9:
        ind = pool
    else:
        ind = np.array(random.sample(pool, 8))
    print 'Available sample #:',np.size(pool,0)
    if np.size(pool,0) == 0:
        print 'Sample size too small.'
        print 'Percentage of', veg_ty, ':', cover_m.sum()*1.0/(cover<17).sum()
        print 'Total number of fire events:', fire.sum()
        print 'Number of pixles burnt on current vegetation type:', (fireM&cover_m).sum()

    subfire = fire[ind[:, 0], ind[:, 1], :]
    del fire
    print 'Reading vegetation data according to mask now...'
    open(tile + 'EVImissing.txt', 'w').close()
    subEVI = np.zeros((np.size(ind, 0), (endY - startY) * annuN))
    subNBR = np.zeros((np.size(ind, 0), (endY - startY) * annuN))
    subEVI[:] = np.nan
    subNBR[:] = np.nan
    l = 0
    for dp in daterange16(startY, endY):
        flnmE = direV + 'MOD13A1.' + dp.strftime('%Y') + dp.strftime('%j') + tile + 'EVI.npy'
        flnmN = direV + 'MOD13A1.' + dp.strftime('%Y') + dp.strftime('%j') + tile + 'NBR.npy'
        if os.path.isfile(flnmE):
            subEVI[..., l] = np.load(flnmE)[ind[:, 0], ind[:, 1]]
            subNBR[..., l] = np.load(flnmN)[ind[:, 0], ind[:, 1]]
        else:
            myfile = open(tile + 'EVImissing.txt', 'a')
            myfile.write(dp.strftime('%Y') + dp.strftime('%j') + tile + '; occurred for ' + '0 times.\n')
            myfile.close()
        l += 1

    np.savez(tile + veg_ty+ 'Samp.npz', subEVI = subEVI, subNBR = subNBR, subfire = subfire, s=s)
    return subfire, subEVI, subNBR, s
####################################################################################################
def fireSampPlt(tile, veg_ty = 'savanna', fire = None, out = ''):
    subfire, subEVI, subNBR, s = fireSamp(tile, veg_ty, fire)

    ### x axis of vegetation time series
    x_v = []
    for dp in daterange16(startY, endY):
        x_v = np.append(x_v, dp.toordinal())
    x_v = x_v[s:]
    
    ### ticks and labels on x axis
    xt = []
    xl = []
    k = 0
    for y in range(startY + 1, endY):
        xt = np.append(xt, datetime.date(y, 1, 1).toordinal())
        if k%((endY - startY) / 3) == 0:
            xl = np.append(xl, str(y))
        else:
            xl = np.append(xl, '')
        k += 1

### Plotting
    print 'Plotting samples now...'
    mpl.use('PDF')
    import matplotlib.pyplot as plt, prettyplotlib as ppl

    fig = plt.figure()
    f = []
    E = []
    N = []
    for p in range(np.size(subEVI,0)):
        fireSeri = subfire[p, :]
        plt.subplot(4, 2, p)
        k = 0
        for dp in daterange16(startY,endY):
            if fireSeri[k]&(k>=s):
                f = plt.axvline((dp-datetime.timedelta(days=8)).toordinal(), color='yellow', linewidth=1.5, zorder=0)
            k += 1
        E, = ppl.plot(x_v, subEVI[p, s:], color='blue', linewidth=0.9)
        N, = ppl.plot(x_v, subNBR[p, s:], color='black', linewidth=0.3)
        plt.xticks(xt, xl, size='small')
        plt.locator_params(axis='y', nbins=5)
        plt.xlim(x_v[0],x_v[-1])
    fig.suptitle('Samples from ' + veg_ty + ' @' + tileCord(tile), fontsize=14)
    leg = fig.legend((f, E, N), ('fire', 'EVI', 'NBR'), loc='lower center', ncol=3, shadow=True)
    leg.draw_frame(False)
    plt.savefig(out + tile + veg_ty + 'EVISamp' + '.pdf', dpi=300)
    plt.close()

####################################################################################################
def cover_freq(tile, cover=None, out=''):
    flnm = tile + '_cover_freq.npy'
    if os.path.isfile(flnm):
        return np.load(flnm)

    if cover is None:
        cover = cover_read(tile)
    
    his = np.bincount(cover.flat)
    freq = np.argsort(his)[::-1]
 
    mpl.use('PDF')
    import matplotlib.pyplot as plt, prettyplotlib as ppl
    
    #plt.figure(figsize=(11, 6), dpi=300) 
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11, 6))
    ppl.bar(ax, range(len(his)), his, width = 1)
    plt.xticks(range(len(his)), sorted(veg_dict , key=veg_dict.get), rotation=30,  fontsize=14)
    plt.ylabel('Pixel Number',  fontsize=14)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title('Vegetation Type ' + ' @' + tileCord(tile), fontsize=16)
    plt.tight_layout() 
    plt.savefig(out+tile+'_cover_freq.pdf', dpi=300)
    plt.close()

    np.save(flnm,freq)
    return freq

####################################################################################################
def tile_finder(veg_ty_list, h_v = None, out = ''):
    veg_ty_list.sort() #standardize veg_ty_list
    flnm = 'tile_finder'
    for i in range(len(veg_ty_list)):
        flnm += ('_'+veg_ty_list[i])
    flnm += '.npy'

    f = 0
    if h_v is None:
        if os.path.isfile(flnm): #return previous data if non-specified
            tile_list = list(np.load(flnm))
            if len(tile_list) is not 0:                 
                myfile = open(out + 'tile_finder.txt', 'a')
                myfile.write('-'.join(veg_ty_list + [': '] + tile_list + ['\n']))
                myfile.close
            return tile_list
        
        f = 1 #non-specified but no previous data
        #cPath = os.path.dirname(os.path.realpath(sys.argv[0]))
        cPath = os.path.dirname(os.path.realpath(__file__))
        tile_co_list = np.loadtxt(cPath+'/tilelist.txt')
        h_v = tile_co_list[:,0:2].astype(np.int)

    tile_list = []
    for i in range(h_v.shape[0]):
        tile = 'h' + str(h_v[i,0]).zfill(2) + 'v' + str(h_v[i,1]).zfill(2)
        cover = cover_read(tile)
        if cover is None:
            continue
        count = np.bincount(cover.flat)
        ind = np.argsort(count)[::-1]
        veg_num_list = np.array([veg_dict[veg_ty_list[p]] for p in range(len(veg_ty_list))])
        if (np.sort(ind[0:len(veg_ty_list)]) == np.sort(veg_num_list)).all():
            if count[ind[len(veg_ty_list)]]>200:
                tile_list.append(tile)

    if f: #save if h_v is non-specified
        np.save(flnm,tile_list)
    if len(tile_list) is not 0: 
        tile_list = list(tile_list)
        myfile = open(out + 'tile_finder.txt', 'a')
        myfile.write('-'.join(veg_ty_list + [': '] + tile_list + ['\n']))
        myfile.close
    
    return tile_list

####################################################################################################
def fire_freq(tile, veg_ty = 'savanna', fire = None, sea_s = 0, out=''):
    sea_s = np.min([sea_s,annuN-1])
    flnm = tile+veg_ty+'_fire_freq.npy'
    if os.path.isfile(flnm):
        fire_freq = np.load(flnm)
        seaso = fire_freq[sea_s]
        return seaso

    if fire is None:
        fire = fire_read(tile)
    veg_num = veg_dict[veg_ty]
    cover = cover_read(tile)
    cover_m = cover == veg_num
    cover_freq(tile, cover=cover)
    del cover

    fireMa = fire&np.tile(cover_m.reshape(2400,2400,1), (1,1,np.size(fire,2)))
    freq = np.sum(fireMa[...,annuN:].reshape((2400, 2400, -1, annuN)), axis=2)    
    freq_t = freq.sum(axis=1).sum(axis=0).flatten()#*1.0/np.sum(cover_m)/(endY-startY)
    if not freq_t.any():
        print 'Sample size too small.'
        print 'Number of fire events on', veg_ty, ':', fireMa.sum()
        return -1
    del freq

    mpl.use('PDF')
    import matplotlib.pyplot as plt, prettyplotlib as ppl
  
    x_a = []
    for dp in daterange16(startY, startY+1):
        x_a = np.append(x_a, (dp-datetime.timedelta(days=8)).toordinal())
    
    xt = []
    for m in range(1, 13):
        xt = np.append(xt, datetime.date(startY, m, 1).toordinal())
    
    plt.figure(figsize=(11, 6), dpi=300)
    ppl.bar(x_a, freq_t, width = 16)
    plt.xticks(xt, ('Jan', '', 'Mar', '', 'May', '', 'Jul', '', 'Sep', '', 'Nov', ''))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel('Burning Frequency')
    plt.suptitle('Fire Seasonality of ' + veg_ty + ' @' + tileCord(tile), fontsize=16)
    plt.savefig(out+tile + veg_ty + 'fire_freq' + '.pdf', dpi=300)
    plt.close()
    
    fire_freq = freq_t.argsort()[::-1]
    seaso = fire_freq[sea_s]
    np.save(flnm, fire_freq)
    return seaso

####################################################################################################
### fire window lenth: 3*16
def samp_cal(tile, veg_ty, fire = None, flag=[0,0,0]):
    winLen = 3
    qua, samp_po, sea_s = flag
    
    flnm = tile + veg_ty +'_samp_' + str(sea_s) + '.npz'
    if os.path.isfile(flnm):
        with np.load(flnm) as data:
            ds = sorted(data.iteritems())
        return [v[1] for v in ds]
    
    print 'Calculating sample mask of ', tile, veg_ty
    if fire is None:
        fire = fire_read(tile)
    veg_num = veg_dict[veg_ty]
    cover = cover_read(tile)
    cover_m = cover == veg_num #land cover mask
    del cover    
    s = 0 # Start date number of data. Assuming fire data starts later than vegetation data.
    while fire[...,s].all(): 
        s += 1
    seaso = fire_freq(tile, veg_ty=veg_ty, fire=fire, sea_s = sea_s) #fire season
    samp_m = (~fire[:,:,s:2*annuN+seaso-winLen//2]).prod(axis=2, dtype=np.bool)&(~(~fire[:,:,2*annuN+seaso-winLen//2:2*annuN+seaso+winLen//2+1]).prod(axis=2, dtype=np.bool))&(~fire[:,:,2*annuN+seaso+winLen//2+1:]).prod(axis=2, dtype=np.bool)&cover_m # only burnt in seaso at 2002
    con_m = (~samp_m)&cover_m
    np.savez(flnm, con_m = con_m, samp_m = samp_m)
    return con_m, samp_m

####################################################################################################
## Calculate post-fire trajectory of EVI and NBR. Will transform EVI, NBR, fire and landcover data if not already transformed. Also calls fire_freq to plot fire frequecy of each month.
## Input: 
###### tile number and vegetation type (default as 'savanna');
###### flag:
########## EVI Reliability: good data v.s. marginal data 
########## Sample Pool: constant sample pixels v.s. variable sample pixels
########## Fire Season Sequence: number of 16-period selected
## Output:
###### 1. histogram of sample size vs. time series length; 2. mean, upper 10% and lower 10% of EVI and NBR; 3. month of fire seaso. EVI and NBR include one year of pre-fire data in the beginning. Output is saved as npz file.
## Dependecy:
###### EVI_read, fire_read, fire_freq, cover_read, daterange16

def traj_cal(tile, veg_ty, fire = None, flag = [0, 0, 0]): 
    qua, samp_po, sea_s = flag

    flnm = tile + veg_ty +'_traj_' + str(qua) + str(samp_po) + str(sea_s) + '.npz'
    if os.path.isfile(flnm):
        with np.load(flnm) as data:
            ds = sorted(data.iteritems())
        return [v[1] for v in ds]
    
    print tile, veg_ty, 'post-fire trajectory.'  
### Prepare data
    if fire is None:
        fire = fire_read(tile)
      
    ## Start date of data. Assuming fire data starts later than vegetation data.
    #s = 0 # Start date number of data. Assuming fire data starts later than vegetation data.    
    #while fire[...,s].all(): 
    #    s += 1
    
    if samp_po==0: #constant sample
        con_m, samp_m = samp_cal(tile, veg_ty, fire = fire, flag = flag)
        rowN = samp_m.sum()
        print 'Sample size: ' + str(rowN)
        EVI = [np.empty((endY-startY)*annuN), np.empty([2,(endY-startY)*annuN]), np.empty((endY-startY)*annuN)]
        NBR = [np.empty((endY-startY)*annuN), np.empty([2,(endY-startY)*annuN]), np.empty((endY-startY)*annuN)]  
        for i in range(len(EVI)):
            EVI[i][:] = np.nan
            NBR[i][:] = np.nan
        his = np.empty((endY-startY)*annuN)
        his[:] = np.nan
        l = 0
        for dp in daterange16(startY, endY):
            print dp      
            EVIraw, NBRraw = EVI_read(tile, dp, qua = qua)
            EVItem = EVIraw[samp_m]
            NBRtem = NBRraw[samp_m]
            mask = ~np.isnan(EVItem)
            his[l] = mask.sum()
            EVIcon = EVIraw[con_m]
            EVI[1][1,l] = np.nanmean(EVIcon)
            NBRcon = NBRraw[con_m]
            NBR[1][1,l] = np.nanmean(NBRcon)
            if his[l]:
            ## EVI
                EVI[0][l] = np.percentile(EVItem[mask], 90)
                EVI[1][0,l] = np.mean(EVItem[mask])
                EVI[2][l] = np.percentile(EVItem[mask], 10)
            ## NBR
                NBR[0][l] = np.percentile(NBRtem[mask], 90)
                NBR[1][0,l] = np.mean(NBRtem[mask])               
                NBR[2][l] = np.percentile(NBRtem[mask], 10)
            l += 1
        colA = -np.zeros((2400,2400), dtype=np.int)
        colA[samp_m] = (endY - startY) * annuN

    else: #varibale sample
        seaso = fire_freq(tile, veg_ty, fire, sea_s) #fire season
        veg_num = veg_dict[veg_ty]
        cover = cover_read(tile)
        cover_m = cover == veg_num #land cover mask
        del cover
        preLen = np.int(((preLen_min-(seaso-s))//annuN+1)*annuN+(seaso-s))  #Length of pre-fire time-series
        print 'Length of pre-fire time-series:', preLen

        rowN = np.sum(np.any(fire[:, :, s+preLen::annuN],axis=2)&cover_m)
        print 'Estimated sample size: < ' + str(rowN)    
        if rowN==0:
            print 'Sample size too small.'
            return
    
        fireTotal = np.empty([2, rowN, (endY - startY) * annuN - s])
        fireTotal[:] = np.nan
        colA = -np.ones((2400, 2400), dtype=np.int)
        rowA = -np.ones((2400, 2400), dtype=np.int)
        samp_m = np.zeros((2400, 2400), dtype=np.bool) #Pixels alive.
        EVItem = np.empty([2400, 2400, preLen]) #EVI history data for preLen
        NBRtem = np.empty([2400, 2400, preLen]) #NBR history data for preLen 
        l = 0
    
        for dp in daterange16(startY, endY):
            print dp
            fire_c = fire[..., l]
            burnMa = fire_c & cover_m
            samp_m[burnMa] = False #Pixel killed 
            colA[samp_m] += 1 #Alive pixels
            if (l%annuN == seaso)&(l>=(s+preLen)): # Fire seaso
                preM = np.sum(fire[:, :, l-preLen:l], axis=2) == 0 #Didn't burn for preLen 16-days
                newPxl = burnMa & preM & (colA==-1) #New born pixels, never lived before
                if newPxl.any():
                    ## Creat new samples and update row and column
                    samp_m[newPxl] = True #Pixels alive
                    rowA[newPxl] = np.arange(np.sum(newPxl)) + np.amax(rowA) + 1 #Start new space for new born
                    colA[newPxl] = preLen #Start to record as the preLen-th 16-day 
                    ## Load pre-fire EVI and NBR
                    fireTotal[0, rowA[newPxl], 0:preLen] = np.roll(EVItem,preLen-(l-1)%preLen-1,axis=2)[newPxl, :]
                    fireTotal[1, rowA[newPxl], 0:preLen] = np.roll(NBRtem,preLen-(l-1)%preLen-1,axis=2)[newPxl, :]
            EVI, NBR = EVI_read(tile, dp, qua = qua)
            EVItem[..., l%preLen] = EVI
            NBRtem[..., l%preLen] = NBR
            fireTotal[0, rowA[samp_m], colA[samp_m]] = EVI[samp_m]
            fireTotal[1, rowA[samp_m], colA[samp_m]] = NBR[samp_m]
            l += 1    

    np.savez(flnm, EVI = EVI, NBR = NBR, colA = colA, his = his)

    return EVI, NBR, colA, his

####################################################################################################
def trajMap(tile, veg_ty = 'savanna', fire=None, out=''):
    preLen, s, _, _, _, _, colA = traj_cal(tile, veg_ty=veg_ty, fire = fire)
    cover = cover_read(tile)
    cover_m= ~(cover==veg_dict[veg_ty])
    tree = tree_read(tile)

    mpl.use('PDF')
    import matplotlib.pyplot as plt
    from matplotlib import colors

    fig, axes = plt.subplots(nrows=2, ncols=7, figsize=(10, 4.3), dpi=200, sharex=True, sharey=True)
        
    m_cmap = mpl.cm.YlGn
    m_cmap.set_bad(color='k')
    #m_cmap.set_under(color='k')
    
    c_cmap = plt.cm.Set3
    bounds = np.array(range(-1,17))+.5
    norm = colors.BoundaryNorm(bounds, c_cmap.N)
    #c_cmap.set_bad(color='k')
    #c_cmap.set_over(color='k', alpha=None)

    for y in range(0,13):
        samp = np.nonzero(colA>preLen+(y-1)*annuN)  
        #tree = tree_per
        #tree.data[tree.mask]=-100
        #tree.mask=colA>=preLen+(y-1)*annuN
        im1 = axes.flat[y].imshow(tree, origin='upper', interpolation='nearest', cmap=m_cmap, clim=(0,100))
        axes.flat[y].scatter(samp[1], samp[0] , s=1, c='r', alpha=.1, edgecolor = 'none', rasterized = True) #
        #axes.flat[y].set_title(str(y), fontsize=14)
        axes.flat[y].axis('off')
        #axes.flat[y].set_rasterized(True)

    im2 = axes.flat[-1].imshow(cover, origin='upper', interpolation='nearest', cmap=c_cmap, norm=norm)
    axes.flat[-1].axis('off')
    #axes.flat[-1].set_title('Vegetation Type', fontsize=14)
    
    cax1, kw1 = mpl.colorbar.make_axes([ax for ax in axes.flat[:7]], location='top', shrink=.7, aspect=50, pad=0.1)
    cbar1 = plt.colorbar(im1, cax=cax1, **kw1)
    cbar1.solids.set_edgecolor('face')

    cax2, kw2 = mpl.colorbar.make_axes([ax for ax in axes.flat[7:]], location='bottom', shrink=1.2, aspect=80, pad=0.1)
    cbar2 = plt.colorbar(im2, cax=cax2, norm=norm, ticks=range(17), **kw2)
    cbar2.solids.set_edgecolor('face')
    cbar2.ax.set_xticklabels(sorted(veg_dict, key=veg_dict.get), rotation=15, fontsize=10)

    plt.subplots_adjust(left = .01, bottom = .15, right = .99, top = .85, wspace = .015, hspace = .015)

    fig.suptitle('Sample change of ' + veg_ty + ' @' + tileCord(tile), fontsize=16)
    plt.savefig(out+tile + veg_ty + 'trajMap' + '.pdf', dpi=200)
    plt.close()

####################################################################################################
def tree_freq_cal(tile, fire=None, typ_n=8, bin_n=16):
    veg_num_list = cover_freq(tile)[:typ_n]
    veg_ty_list = sorted(veg_dict, key=veg_dict.get)

    flnm = tile + '_tree_freq_' + str(typ_n) + str(bin_n) + '.pkl'
    if os.path.isfile(flnm):
        return pd.read_pickle(flnm)

    print tile, 'fire frequency - tree cover.'  
    
    cover_total = np.dstack([cover_read(tile, datetime.date(y,1,1)) for y in range (2001,2010)])
    cover,_ = stats.mode(cover_total,axis=2)
    del cover_total
    cover = cover[...,0]

    tree_total = np.dstack([tree_read(tile, datetime.date(y,1,1)) for y in range (2001,2010)])
    tree = np.nanmean(tree_total,axis=2)    
    
    tree_ave_list = [np.nanmean(tree[cover==veg_num]) for veg_num in veg_num_list] #treecover average
    veg_num_sort = [veg_num for (tree_ave,veg_num) in sorted(zip(tree_ave_list,veg_num_list), reverse=True)] #sort according to treecover
    
    tree_min = np.nanmin(tree.flat)
    tree_max = np.nanmax(tree.flat)
    clas_mat = ((tree-tree_min)//((tree_max-tree_min)/bin_n+.01)).astype(np.int)
    
    nan_m = np.isnan(tree)
    tree[nan_m] = 0
    del tree_total

    if fire is None:
        fire = fire_read(tile)
    fire_num = np.sum(fire[..., annuN:], axis=2)
    del fire
    
    tree_freq = collections.OrderedDict({})
    for i in range(typ_n):
        cover_m = (cover == veg_num_sort[i])&(~nan_m)
        clas_tem = clas_mat[cover_m]
        fire_tem = fire_num[cover_m]

        count_num = np.bincount(clas_tem, minlength=bin_n).astype(np.float)
        count_num[count_num==0] = np.nan
        fire_sum = np.bincount(clas_tem, weights=fire_tem, minlength=bin_n)
        freq = np.append(fire_sum/count_num, np.mean(fire_tem))/(endY-startY-1)
        tree_freq[veg_ty_list[veg_num_sort[i]]] = freq[::-1]
    
    df = pd.DataFrame(tree_freq, index=['Total']+['{:.1f}'.format(tree_min+(0.5+i)*(tree_max-tree_min)/bin_n) for i in range(bin_n,0,-1)])
    df.to_pickle(flnm)
    
    return df
####################################################################################################
def tree_freq_plt(tile, fire=None, typ_n=8, bin_n=16, out=''):
    df = tree_freq_cal(tile, fire=None, typ_n=typ_n, bin_n=bin_n)
    
    mpl.use('PDF')
    import matplotlib.pyplot as plt
    
    c_cmap = plt.cm.binary
    c_cmap.set_gamma(.4)
    c_cmap.set_bad(color='yellow')
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 2.3), dpi=200, sharex=True, sharey=True)
    plt.imshow(ma.array(df.values, mask = np.isnan(df.values)).T, interpolation='nearest', cmap=c_cmap).get_axes()
    ax.set_xticks(np.linspace(0, bin_n, bin_n+1))
    ax.set_xticklabels(df.index, fontsize=8, rotation= 'vertical')
    ax.set_xlabel('Tree Cover (%)',fontsize=9) 
    
    ax.set_yticks(np.linspace(0, typ_n-1, typ_n))
    ax.set_yticklabels(df.columns,fontsize=8)
    ax.grid('off')
    
    ## set position for colorbar
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = plt.colorbar(shrink=0.8)
    tick_locator = mpl.ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    ## get rid of white lines
    cbar.solids.set_edgecolor('face')
    cbar.set_label('Burnt/year', fontsize=9)
        
    #for i in range(m):
    #    for j in range(n):
    #    ax.text(j, i, '{:.2f}'.format(df.iget_value(i, j)),
    #            size='medium', ha='center', va='center',
    #            path_effects=[patheffects.withSimplePatchShadow(shadow_rgbFace=(1,1,1))])
    #plt.figlegend(p, lab, loc='center right', prop={'size':9})
    plt.title('Fire frequency and Treecover @' + tileCord(tile), fontsize=10)
    #plt.subplots_adjust(left = .1, bottom = .25, right = .9, top = .85, wspace = .015, hspace = .015)
    plt.tight_layout()
    plt.savefig(out + tile + '_tree_freq_' + str(typ_n) + str(bin_n) + '.pdf', dpi=300)
    plt.close()

####################################################################################################

def traj_an_plt(tile, veg_ty = 'savanna', fire=None, flag = [0,0,0], out=''):
    preLen,s,_,mea,_,_,nonfireMea,_ = traj_cal(tile, veg_ty=veg_ty, fire = fire, flag = flag) 
    qua = flag[0]
    samp_po = flag[1]
    sea_s = flag[2]

    import matplotlib.pyplot as plt
    
    x_v = [dp.toordinal() for dp in daterange16(startY, endY)][s+preLen:s+preLen+annuN] #x-axis
    xt = [(date(date.fromordinal(x_v[0]).year, date.fromordinal(x_v[0]).month+1, 1)+relativedelta(months=m*2)).toordinal() for m in range(7)]
    xl = [(date(date.fromordinal(x_v[0]).year, date.fromordinal(x_v[0]).month+1, 1)+relativedelta(months=m*2)).strftime('%b') for m in range(7)]
    
    fire_d = date.fromordinal(np.int(x_v[0]-8)) #Center of 16-day period
    
    fig, ax = plt.subplots(1)
    plt.xticks(xt, xl)
    plt.ylim(mea[0,:].min(),mea[0,:].max())
    plt.xlim(x_v[0],x_v[-1])
    col_s = np.linspace(0, 1.0, num=(mea.shape[1]-preLen)//annuN)
    
    min_in = 1
    max_in = 0
    for i in range((mea.shape[1]-preLen)//annuN+1):
        x = np.array(x_v)
        y = mea[0,preLen+(i-1)*annuN:preLen+i*annuN]
        #mask = ~np.isnan(y)
        if i == 0:
            line_c = 'k'
        else:
            line_c = plt.cm.jet(col_s[i-1])
        plt.plot(x, y, label = str(i), color = line_c, linewidth=2.5, alpha = .5)
        min_in = np.min([min_in,np.nanmin(y)])
        max_in = np.max([max_in,np.nanmax(y)])
    plt.ylim(min_in,max_in)
    
    plt.legend(loc='best', framealpha = .3)
    fig.suptitle('Post-fire seasonal cycle of ' + veg_ty + ' burnt in ' + fire_d.strftime('%b') + fire_d.strftime('%d') + '$\pm$8 days @' + tileCord(tile) + '\nData reliability: ' + qua_dic[qua] + '; Sample Pool: ' + samp_dic[samp_po] + '; Fire Season: number ' + str(sea_s+1), fontsize=12)

    plt.savefig(out + tile + veg_ty + 'An_' + str(qua) + str(samp_po) + str(sea_s) + '.pdf', dpi=300)
    plt.close()
####################################################################################################
## Plot post-fire trajectory using result from function traj_cal, in PDF format. Will call traj_cal if required data is not found.
## Input:
###### tile number and vegetation type.
## Output: 
###### 1. Histogram of sample size; 2. Post-fire trajectory of EVI and NBR.
## Dependency:
###### traj_cal, daterange16

def traj_plt(tile, veg_ty, fire=None, flag = [0,0,0], out=''):
    EVI, NBR, _, his = traj_cal(tile, veg_ty=veg_ty, fire = fire, flag = flag)

    qua, samp_po, sea_s = flag
    seaso = fire_freq(tile, veg_ty, fire, sea_s)
    fire_d = list(daterange16(2002, 2002+1))[seaso]
########## EVI Reliability: good data v.s. marginal data 
########## Sample Pool: constant sample pixels v.s. variable sample pixels
########## Fire Season Sequence: number of 16-period selected
    qua_dic = {0: 'Good',
                1: 'Marginal',
                2: 'Snow/Ice',
                3: 'Cloudy'}
    
    samp_dic = {0: 'Constant',
                1: 'Variable'}

    ##Calculating x-tick and x-label
    x = np.array([dp.toordinal() for dp in daterange16(startY, endY)]) # x-axis
    xt = np.array([(fire_d+relativedelta(years=y)).toordinal() for y in range(startY-2002, endY-2002)]) # x tick
    xl = [str(y) for y in range(startY-2002, endY-2002)] # x label
    f_tick = fire_d-relativedelta(years=2)
    
    ##Plotting
    mpl.use('PDF')
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
    plt.xticks(xt, xl, size='small')
    axes[1].tick_params(axis='x', pad=12)
    #plt.xlim(x[~np.isnan(EVI[1][0,:])][0],x[~np.isnan(EVI[1][0,:])][-1])
    plt.xlabel('Time since fire (year)',fontsize=12)
    
    ### EVI
    A = axes[0].fill_between(x, EVI[0], EVI[2], linewidth=0, alpha=.2)
    B = axes[0].axvline(fire_d.toordinal(), color='yellow', linewidth=1.5, zorder=0)
    M, = axes[0].plot(x, EVI[1][0,:], color='black', linewidth=0.7)
    C, = axes[0].plot(x, EVI[1][1,:], color='green', linewidth=0.4)
    axes[0].set_ylabel('EVI')
    axes[0].yaxis.grid(b=True, alpha=.2)
    U = axes[0].axhline(y=np.nanmax(EVI[1][0, 0:2*annuN]), color='black', linewidth=.7, alpha=.3)
    L = axes[0].axhline(y=np.nanmin(EVI[1][0, 0:2*annuN]), color='black', linewidth=.7, alpha=.3)
    
    ### NBR
    A = axes[1].fill_between(x, NBR[0], NBR[2], linewidth=0, alpha=.2)
    B = axes[1].axvline(fire_d.toordinal(), color='yellow', linewidth=1.5, zorder=0)
    M, = axes[1].plot(x, NBR[1][0,:], color='black', linewidth=0.7)
    C, = axes[1].plot(x, NBR[1][1,:], color='green', linewidth=0.4)
    axes[1].set_ylabel('NBR')
    axes[1].yaxis.grid(b=True, alpha=.2)
    U = axes[1].axhline(y=np.nanmax(NBR[1][0, 0:2*annuN]), color='black', linewidth=.7, alpha=.3)
    L = axes[1].axhline(y=np.nanmin(NBR[1][0, 0:2*annuN]), color='black', linewidth=.7, alpha=.3)
    axes[1].xaxis.tick_top()
    ax2 = axes[1].twinx()
    S, = plt.plot(x, his, linewidth=1,  color='black', alpha = .3)
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.set_ylabel('Sample Size',fontsize=12)
        
    fig.suptitle(veg_ty + ' burnt ' + fire_d.strftime('%b') + fire_d.strftime('%d') + '$\pm$24 @' + tileCord(tile) + '\nData reliability: ' + qua_dic[qua] + '; Sample Pool: ' + samp_dic[samp_po] + '; Fire Season: number ' + str(sea_s+1), fontsize=14)
    leg = fig.legend((B, M, C, S), ('Burning Event', 'Mean and upper/lower 10 percentile','Control Mean','Sample Size'), loc='lower center', ncol=4, prop={'size':10})
    leg.draw_frame(False)

    plt.savefig(out + tile + veg_ty + '_traj_' + str(qua) + str(samp_po) + str(sea_s) + '.pdf', dpi=300)
    plt.close()

####################################################################################################
def cover_ch_cal(tile, veg_ty, fire = None, flag = [0,0,0]): 
    qua, samp_po, sea_s = flag
    flnm = tile + veg_ty +'_cover_ch_' + str(sea_s) + '.npy'
    if os.path.isfile(flnm):
        cover_ch = np.load(flnm)
        return cover_ch
     
    print tile, veg_ty, 'land cover change.' 
    con_m, samp_m = samp_cal(tile = tile, veg_ty = veg_ty, fire = fire, flag=flag)
    
    cover_ch = np.empty([len(veg_dict)+1,2012-2001+1,2])
    cover_ch[:] = np.nan
    
    for y in range(2001, 2013):
        cover_tem = cover_read(tile, datetime.date(y,1,1))
        rec_samp = cover_tem[samp_m]
        con_samp = cover_tem[con_m]
        for vegnum in range(len(veg_dict)):
            cover_ch[vegnum,y-2001,0] = (rec_samp==vegnum).sum()
            cover_ch[vegnum,y-2001,1] = (con_samp==vegnum).sum()
        tree_tem = tree_read(tile, datetime.date(y,1,1))
        cover_ch[len(veg_dict),y-2001,0] = np.nanmean(tree_tem[samp_m])
        cover_ch[len(veg_dict),y-2001,1] = np.nanmean(tree_tem[con_m])
    np.save(flnm, cover_ch)

    return cover_ch
####################################################################################################
def cover_ch_plt(tile, veg_ty, fire = None, flag = [0,0,0], out=''): 
    cover_ch = cover_ch_cal(tile, veg_ty, fire=fire, flag=flag)
    qua, samp_po, sea_s = flag    
    seaso = fire_freq(tile, veg_ty, fire, sea_s)
    fire_d = list(daterange16(2002, 2002+1))[seaso]

    mpl.use('PDF')
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex = True, )
    cover_cum = np.zeros([12,2])
    subtitle = ['Pixels Burnt only in 2002', 'All Other Pixels']
    p = []
    for lay in range(2):
        for vegnum in range(len(veg_dict)):
            bar = axes[lay].bar(range(12), cover_ch[vegnum,:,lay], 1, color = plt.cm.Set3(vegnum*1.0/len(veg_dict)), bottom=cover_cum[:,lay], edgecolor = 'None')
            if lay==0: p.append(bar[0])
            cover_cum[:,lay] = cover_cum[:,lay] + cover_ch[vegnum,:,lay]
        axes[lay].set_title(subtitle[lay], fontsize = 11)
        axes[lay].set_ylim(0,cover_cum[0,lay])
        axes[lay].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax = axes[lay].twinx() #axis for tree cover
        ll = ax.plot(np.arange(12)+.5,cover_ch[-1,:,lay]) #tree cover
        #p.append(ll[0])
        ax.set_ylim(np.nanmin(cover_ch[-1,:,:]),np.nanmax(cover_ch[-1,:,:]))
        ax.set_ylabel('Tree Cover (%)', fontsize=10,color='b')
        ax.yaxis.grid(b=True, color='b', alpha=.2)
        ax.tick_params(axis='y', colors='b', labelsize=10)
        ax.axvline(1.5, color='yellow', linewidth=1.5, zorder=0)

    axes[0].set_xticklabels([])
    axes[1].set_xlabel('Time since fire (year)',fontsize=10)
    axes[0].set_ylabel('Pixel Number',fontsize=10)
    #box = axes[0].get_position()
    #(axes[q].set_position([box.x0, box.y0, box.width * 0.7, box.height]) for q in range(2))
    plt.figlegend(p, sorted(veg_dict , key=veg_dict.get), loc='lower center', ncol=6, prop={'size':9})#, bbox_to_anchor=(1.18, 0.5),  prop={'size':11})
    plt.xticks(np.arange(12)+.5,np.arange(12)-1)
    plt.xlim(0,12)
    plt.subplots_adjust(left = .1, bottom = .2, right = .9, top = .88, wspace = .015, hspace = .2)
    #plt.tight_layout()

    #ax2.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    #plt.legend(sorted(veg_dict , key=veg_dict.get))

    fig.suptitle(veg_ty + ' burnt ' + fire_d.strftime('%b') + fire_d.strftime('%d') + '$\pm$24 (number ' + str(sea_s+1)+ ') @' + tileCord(tile), fontsize=14)
    plt.savefig(out + tile + veg_ty + '_cover_chan_' + str(sea_s) + '.pdf', dpi=300)
    plt.close()
