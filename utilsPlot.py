import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors


def plot_gmfs_hist(gmfs):
    plt.hist(gmfs['gmv'])
    plt.ylabel('Frequency')
    plt.xlabel(gmfs['IMT'][0])

    return plt.show()


def plot_hist(data, title, xlabel, ylabel, bins):
    plt.hist(data,bins=bins)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    
    return plt.show()


def plot_gmfs_hist2(gmfs,bins):
    n, bins, patches = plt.hist(gmfs['gmv'],
                            bins=bins,
                            normed=1,
                            facecolor='green')


def plot_fitted_lognorm(values, xlabel):
    shape = np.std(np.log(values))
    scale = np.exp(np.mean(np.log(values)))
    x = np.linspace(min(values[0]), max(values[0]), num=400)
    plt.plot(x, stats.lognorm.pdf(x, shape, loc=0, scale=scale), 'r', linewidth=3) # Plot fitted curve
    plt.xlabel(xlabel)
    #plt.xticks(fontsize=15)
    plt.show()


def plot_ff(fragility_model):
    #for i in xrange(len(fragility_model['types'])):
    for i in xrange(1):
        if fragility_model['types'][i] == 'continuous':
            plot_continuous_ff(fragility_model['values'][i], fragility_model['IMT'][i])

        if fragility_model['types'][i] == 'discrete':
            plot_discrete_ff(fragility_model['values'][i], fragility_model['IMT'][i])


def plot_continuous_ff(values, IMT):
    fig = plt.figure(figsize=(10,6))
    fig = plt.figure()
    for i in xrange(len(values)):
        stddev = values[i][1]
        mean = values[i][0]
        median = mean**2/np.sqrt(stddev**2+mean**2)
        beta = np.sqrt(np.log(stddev**2/mean**2+1))
        x = np.linspace(0.01, 2.5, num=100)
        cdf = stats.lognorm.cdf(x, beta, loc=0, scale=median) 
        ds1 = plt.plot(x, cdf, label = values[i][2])
    ax = fig.gca()
    ax.set_xticks(np.arange(0,2.51,0.2))
    ax.set_yticks(np.arange(0,1.01,0.2))
    ax.set_xlabel(str(IMT)+' (g)', fontsize=12)
    ax.set_ylabel('Probability', fontsize=12)
    plt.legend(loc=4, fontsize=12)
    plt.grid()

    return plt.show()


def plot_discrete_ff(values, IMT):
	#needs revision and add legend
    fig = plt.figure()
    x = values[0]
    for i in xrange(1,len(values)):
        y = values[i]
        ds1 = plt.plot(x, y, label = values)
    ax = fig.gca()
    ax.set_xticks(np.arange(int(values[0][0]),int(values[0][-1]),1))
    ax.set_yticks(np.arange(0,1.,0.2))
    ax.set_xlabel(IMT)
    ax.set_ylabel('Probability')
    plt.legend(loc=4)
    plt.grid()
    
    return plt.show()


def plot_vulnFunction(vulnerability_model,ylabel):
    fig = plt.figure(figsize=(10,6))
    fig = plt.figure()
    x = map(float, vulnerability_model['values'][0])
    y = map(float, vulnerability_model['values'][1])
    ds1 = plt.plot(x, y)
    ax = fig.gca()

    ax.set_xticks(np.arange(0,np.max(x)+0.1*np.max(x),0.5))
    ax.set_yticks(np.arange(0,np.max(y)+0.1*np.max(y),20))
    ax.set_xlabel(vulnerability_model['IMT'][0]+' (g)', fontsize = 12)
    ax.set_ylabel(ylabel, fontsize = 12)
    plt.grid()
    
    return plt.show()


def plot_gmf_path(gmfs, measure, percentileValue, logScale, save_gmf):
    plt.figure(figsize=(10,6))
    lon1 = float(min(gmfs['lon']))
    lon2 = float(max(gmfs['lon'])) 
    minLon = np.min([lon1,lon2])-0.1
    maxLon = np.max([lon1,lon2])+0.1
    lat1 = float(min(gmfs['lat']))
    lat2 = float(max(gmfs['lat'])) 
    minLat = np.min([lat1,lat2])-0.1
    maxLat = np.max([lat1,lat2])+0.1
    map = Basemap(lat_0=(minLat+maxLat)/2, lon_0=(minLon+maxLon)/2,resolution = 'i', area_thresh = 1000.0,
            llcrnrlon=minLon, llcrnrlat=minLat, urcrnrlon=maxLon, urcrnrlat=maxLat)
    map.drawcountries(linewidth=1)
    map.fillcontinents(color='coral',alpha=0.3)
    map.drawmapboundary(linewidth=1)

    x,y = map(gmfs['lon'], gmfs['lat'])

    avegmfs=np.zeros(shape=(len(gmfs['lon']),3))
    for i in xrange(len(gmfs['lon'])):
        if measure == 'mean':
            avegmfs[i,0] = np.mean(gmfs['gmv'][i])
        elif measure == 'median':
            avegmfs[i,0] = np.median(gmfs['gmv'][i])
        elif measure == 'percentile':
            avegmfs[i,0] = np.percentile(gmfs['gmv'][i], percentileValue)
    mag = avegmfs[:,0]
    if logScale == True:
        plt.scatter(x,y,c=mag,norm=colors.LogNorm(),marker='o',lw = 0,s=70)
    else:
        map.scatter(x,y,c=mag,marker='o',s=70,lw = 0)

    c = plt.colorbar(orientation='vertical')
    if measure == 'percentile':    
        c.set_label("%d th %s gmf values - SA(1.0)" % (percentileValue, measure))
    else:
        c.set_label("%s gmf values - SA(1.0)" % measure)
    #c.set_clim(vmin=0,vmax=0.13)
    
    if save_gmf:
        avegmfs[:,1] = gmfs['lon']
        avegmfs[:,2] = gmfs['lat']
        np.savetxt('./gmf.txt', avegmfs, fmt = '%.6f' )


def filter_list(data, limitValue):
    dataFiltered = []
    dataFiltered = map(lambda x: limitValue if x<limitValue else x, data)
    
    return dataFiltered


def filter_nested_lists(data, limitValue):
    dataFiltered = []
    for i in xrange(len(data)):
        data1 = map(lambda x: limitValue if x<limitValue else x, data[i])
        dataFiltered.append(data1)
    
    return dataFiltered


def plot_country(data, exposure_model, measure, percentileValue, logScale):
    plt.figure(figsize=(10,6))
    lon1 = float(min(exposure_model['lon']))
    lon2 = float(max(exposure_model['lon'])) 
    minLon = np.min([lon1,lon2])-0.2
    maxLon = np.max([lon1,lon2])+0.2
    lat1 = float(min(exposure_model['lat']))
    lat2 = float(max(exposure_model['lat'])) 
    minLat = np.min([lat1,lat2])-0.2
    maxLat = np.max([lat1,lat2])+0.2
    map = Basemap(lat_0=(minLat+maxLat)/2, lon_0=(minLon+maxLon)/2,resolution = 'i', area_thresh = 1000.0,
            llcrnrlon=minLon, llcrnrlat=minLat, urcrnrlon=maxLon, urcrnrlat=maxLat)
    
    map.drawcountries(linewidth=1)
    map.fillcontinents(color='coral',alpha=0.3)
    map.drawmapboundary(linewidth=1)

    parallels = np.arange(37.,43.5,1.)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=10, linewidth=0.1)
    meridians = np.arange(-9,-6,1.)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10, linewidth=0.1)

    x,y = map(exposure_model['lon'], exposure_model['lat'])

    aveData = np.zeros(shape=(len(exposure_model['lat']),1))
    for i in xrange(len(exposure_model['lat'])):
        aveData[i,0] = np.median(data[i])
        #aveData[i,0] = np.median(np.log(data[i]))
    mag = aveData[:,0]

    if logScale == True:
        plt.scatter(x,y,c=mag,norm=colors.LogNorm(),marker='o',lw = 0,s=20)
    else:
        map.scatter(x,y,c=mag,marker='o',s=20,lw = 0)
    
    c = plt.colorbar(orientation='vertical')
    c.ax.tick_params(labelsize=10)
    #c.set_label("median log gmf values - Sa(1sec)")
    plt.savefig('./destination_path.tiff', format='tiff', dpi=200)

