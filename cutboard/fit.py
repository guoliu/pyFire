from pyEarth import GIS, plot
import numpy as np, pandas as pd, os
from scipy.stats.mstats import mquantiles
kai = 80.0

####################################################################################################
def cleaner(data):
    mask = np.ones(data[0].shape, dtype=np.bool)
    for ind in data:
        mask[np.isnan(ind)] = False
    return [ind[mask] for ind in data]

####################################################################################################
def grow(xx,p50,s,d):
    return kai/(1 + np.exp(-4*s*(xx-p50)/kai))*(d/(d+xx))

####################################################################################################
def dataLoad(region):
    tree = GIS.read('MOD44B.'+region+'.tif')
    p_wet = GIS.read('wetSum.'+region+'.tif')
    alp_wet = GIS.read('wetMean.'+region+'.tif')
    
    return cleaner([tree, p_wet, alp_wet])

##################################################################################################
def Rvalue(para, xx, yy, a, tao=0.99): 
    from scipy.stats.mstats import mquantiles
    s = para[0]+para[1]*a
    p50 = para[2]+para[3]*a
    d = para[4]
    Q = mquantiles(yy,prob=[tao])
    ################################################## 
    def rho(x):
        return x*tao*(x>0)-x*(1-tao)*(x<=0)
    ################################################## 
    y_e = grow(xx, p50, s, d)
    R = 1-rho(yy-y_e).sum()/rho(yy-Q).sum()
    return (R,)

####################################################################################################
def evolve(a,xx,yy):
    import random
    from random import normalvariate
    from deap import base, creator, tools, algorithms 
    
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    
    toolbox = base.Toolbox()
    ######
    def seeding():
        para = [normalvariate(0.2,0.2), normalvariate(-0.004,0.04), normalvariate(300,200), normalvariate(22,30),normalvariate(10000,10000)]
        if ((para[0]+para[1]*a)>0)&((para[2]+para[3]*a)>0):
            return para
        else:
            return seeding()
    ######
    def tortilla(turkey, lettuce):
        return turkey(lettuce())
    ######
    toolbox.register("individual", tortilla, creator.Individual, seeding)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    
    toolbox.register("evaluate", Rvalue, a=a, xx=xx, yy=yy)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)
    
    pop = toolbox.population(n=800)
    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)
    
    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=40, stats=stats, halloffame=hof, verbose=True)
    
    return tools.selBest(pop,1)

####################################################################################################
def train(region):    
    yyTot, xxTot, alp = dataLoad(region) 
    
    ####################################################################################################
    d = 1.0
    stops = np.arange(d/2.0,25+d/2.0,d)
    data=np.empty((len(stops),9)) #a, para[0:5],s,p
    data[:] = np.nan
    i = 0
    for a in stops:
        mask=(alp>a-d/2.0)&(alp<a+d/2.0)
        data[i,0] = a
        if mask.sum()>10:
            xx = xxTot[mask]
            yy = yyTot[mask]
            res = evolve(a,xx,yy)[0]
            data[i,1:6] = res
            data[i,6] = res[0]+res[1]*a #s
            data[i,7] = res[2]+res[3]*a #P50
            data[i,8] = res.fitness.values[0]
        else:
            data[i,1:] = np.nan
        i += 1
    import pandas as pd
    table = pd.DataFrame(data[:,1:], columns=['c1','c2','c3','c4','d','s','P50','fitness'],index=data[:,0])
    table.to_pickle('para.' + region + '.pkl')
    ax = table[['s','P50','fitness']].plot(secondary_y=['P50'])
    ax.figure.savefig('para.' + region + '.png')
    ax.figure.clf()
    return table

#####################################################################################################

####################################################################################################
def paraLoad(region):
    table = pd.read_pickle('para.' + region + '.pkl')
    cleanTable = table[table['fitness']>0 & table['fitness'].apply(np.isnan)]
    para = [(cleanTable.iloc[:,c]*(cleanTable['fitness'])/cleanTable['fitness'].sum()).sum() for c in range(5)]
    return para

####################################################################################################
def boost():
    buf = 20
    tree = GIS.read('MOD44B.tif')
    p_wet = GIS.read('wetSum.tropics.tif')
    alp = GIS.read('wetMean.tropics.tif')
    tree_po = tree.copy()
    mar = ((90-23.25)*4,(90+23.25)*4)

    limits = {'america':(0,600),'africa':(600,980),'asiaPaci':(0,1440)}
    for block in limits:
        para = paraLoad(block)
        s = para[0]+para[1]*alp[:,limits[block][0]:limits[block][1]]
        p = para[2]+para[3]*alp[:,limits[block][0]:limits[block][1]]
        tree_regPo = grow(p_wet[:,limits[block][0]:limits[block][1]],p,s,para[4]) #regional potential tree cover

        #####Gradient patch
        bufNorth = tree_po[mar[0]:mar[0]+buf,limits[block][0]:limits[block][1]]
        bufSouth = tree_po[mar[1]-buf:mar[1],limits[block][0]:limits[block][1]]
        
        grad = np.repeat(np.linspace(0,1,buf+2)[1:-1].reshape((-1,1)),tree_regPo.shape[1],axis=1)
        tree_regPo[:buf,:] = grad*tree_regPo[:buf,:] + (1-grad)*bufNorth
        tree_regPo[-buf:,:] = (1-grad)*tree_regPo[-buf:,:] + grad*bufSouth
        tree_po[mar[0]:mar[1],limits[block][0]:limits[block][1]] = tree_regPo
    
    tree_po[np.isnan(tree)] = np.nan
    tree_po[:,0:240][tree[:,0:240]==0] = np.nan #fix missing tile
    GIS.write(tree_po, 'MOD44B.boost.tif','MOD44B.tif')
    plot.mapDraw('MOD44B.boost.tif', 'Potential Tree Cover (%)', vMin=0, vMax=80, cmNm='terrain_r')
    
    os.system('gdalwarp -te -180 -90 180 90 -tr 1.25 0.9375 -r average -overwrite MOD44B.boost.tif MOD44B.boost.lowRe.tif')

####################################################################################################
def fitPlot(region,local=False):
    d = 1.0
    yyTot, xxTot, alp = dataLoad(region) 

    import matplotlib.pyplot as plt
    almost_black = '#262626'
    fig, ax = plt.subplots(1)
    cList = ['maroon','darkgreen','navy']
    para = paraLoad(region)
    di = para[4]
    if local:
        table = pd.read_pickle('para.' + region + '.pkl')
        best = table['fitness'].order(ascending=False).head(3).sort_index(ascending=False) #3 highist fitness, sorted by a
        a = list(best.index[:])
        fits = list(best)
        label= [r'$\alpha_w$='+str(a[i])+', fitness='+"{0:.2f}".format(fits[i]) for i in range(3)]
        p50 =  [table['P50'][a[i]] for i in range(3)]
        s = [table['s'][a[i]] for i in range(3)]
        figNm = 'local/fit.' + region + '.local.png'
    else:
        a = [14.5, 9.5, 4.5]
        label = [r'$\alpha_w$='+str(a[i]) for i in range(3)]
        s = [para[0]+para[1]*a[i] for i in range(3)]
        p50 = [para[2]+para[3]*a[i] for i in range(3)]
        figNm = 'fit.' + region + '.png'
    for i in range(len(a)):
        mask=(alp>a[i]-d/2.0)&(alp<a[i]+d/2.0)
        xx = xxTot[mask] #p_wet
        yy = yyTot[mask] #tree
        ax.scatter(xx, yy, alpha=0.6, facecolors='none', edgecolor=cList[i])
        ax.plot(np.arange(0,2500), grow(np.arange(0,2500), p50[i], s[i], di), color=cList[i], linewidth=0.7, label=label[i], alpha=0.7)
    ax.set_xlabel('Mean Wet Season Precipitation (mm)', fontsize=14)
    ax.set_ylabel('Tree Cover (%)', fontsize=14)
    ax.set_xlim([0,2500])
    ax.set_ylim([0,90])  
    ax.legend(frameon=True, fontsize=10, loc='best')
    fig.savefig(figNm, dpi = 500)
    plt.close()

####################################################################################################
for region in ['africa','america','asiaPaci','tropics']:#,'tropic_africa','tropic_america','tropic_asiaPaci']:    
    train(region)
    fitPlot(region)
    fitPlot(region, local=True)
boost()
'''
for region in ['africa','america','asiaPaci','tropics','tropic_africa','tropic_america','tropic_asiaPaci']:
    train(region)
'''
