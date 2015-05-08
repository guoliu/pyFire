def ktrSm(source, ktrLst=None):
    import json, numpy as np
    from pyEarth import GIS

    with open('countryVal.json') as f:
        ktrDx = json.load(f)

    ktrCrpDx = {}
    if ktrLst is None:
       ktrLst = list(ktrDx)
    crpMp = GIS.read(source)
    ktrCd = GIS.read('/data8/data/guol3/cutboard/countries.tif')
    if crpMp.shape[0]<ktrCd.shape[0]:
        mar = np.empty(((ktrCd.shape[0]- crpMp.shape[0])/2, crpMp.shape[1],))*np.nan
        crpMp = np.concatenate((mar, crpMp, mar), axis=0)
    for ktr in ktrLst:
        mask = ktrCd==ktrDx[ktr]
        ktrCrpDx[ktr] = np.nansum(crpMp[mask])
    with open('/data8/data/guol3/cutboard/ktrCrpDx.json', 'w') as f:
        json.dump(ktrCrpDx, f)
    
    return ktrCrpDx

ktrCrpDx = ktrSm('GLOBCOVER_L4_200901_200912_V2.3.totaCrop.world.weighted.tif')
