from scipy import stats
from lxml import etree
import numpy as np
import copy
import matplotlib.pyplot as plt
import glob
import os
import scipy
from scipy.linalg import toeplitz


def read_gmfs_file(gmfs_file, exposure_model):

    lon = []
    lat = []
    IMT = []
    num_gmfs = []
    rupture_id = []
    nrml = '{http://openquake.org/xmlns/nrml/0.5}'
    lon = exposure_model['lon']
    lat = exposure_model['lat']
    coords = []
    for i in xrange(len(lon)):
        coords.append((lon[i],lat[i]))

    num_gmfs = 0
    for _, element in etree.iterparse(gmfs_file):
        if element.tag == '%sgmf' % nrml:
            num_gmfs += 1
            element.clear()
    gmv = [ ([0] * num_gmfs) for i in range(len(lon)) ]
    j=0
    for _, element in etree.iterparse(gmfs_file):
        if element.tag == '%sgmf' % nrml:
            imt = element.attrib['IMT']
            IMT.append(imt)
            ruptureId = element.attrib['ruptureId']
            rupture_id.append(ruptureId)
            for node in element:
                gmvi = float(node.attrib['gmv'])
                loni = format(float(node.attrib['lon']),'.6e')
                lati = format(float(node.attrib['lat']),'.6e')           
                gmv[coords.index((loni,lati))][j] = gmvi
                node.clear()
            j+=1
            element.clear()

    gmfs = {'gmv': None, 'IMT': None, 'lon': None, 'lat': None, 'num_gmfs': None, 'rupture_id': None}   
    gmfs['gmv'] = gmv
    gmfs['IMT'] = IMT
    gmfs['lon'] = lon
    gmfs['lat'] = lat
    gmfs['num_gmfs'] = num_gmfs
    gmfs['rupture_id'] = rupture_id

    return gmfs


def combine_gmfs(gmfs_file1, gmfs_file2, exposure_model, weight_1, weight_2):
    if weight_1+weight_2 !=1:
        print("The sum of the weights needs to be one")
    else:
        gmfs1 = read_gmfs_file(gmfs_file1, exposure_model)
        gmfs2 = read_gmfs_file(gmfs_file2, exposure_model)
        gmv = []
        weight_a = weight_1*np.asarray(gmfs1['gmv'])
        weight_b = weight_2*np.asarray(gmfs2['gmv'])
        gmv = (weight_a + weight_b).tolist()

        gmfs = {'gmv': None, 'IMT': None, 'lon': None, 'lat': None, 'num_gmfs': None, 'rupture_id': None}   
        gmfs['gmv'] = gmv
        gmfs['IMT'] = gmfs1['IMT']
        gmfs['lon'] = gmfs1['lon']
        gmfs['lat'] = gmfs1['lat']
        gmfs['num_gmfs'] = gmfs1['num_gmfs']
        gmfs['rupture_id'] = gmfs1['rupture_id']

    return gmfs


def read_exposure_model(exposure_file):
    ids = []
    number = []
    taxonomies = []
    lon = []
    lat = []
    types = []
    values = []
    nrml = '{http://openquake.org/xmlns/nrml/0.4}'
    
    tree = etree.iterparse(exposure_file, tag='%sasset' % nrml )
    for event, elem in tree:
        ids.append(elem.attrib['id'])
        number.append(float(elem.attrib['number']))
        taxonomies.append(elem.attrib['taxonomy'])
    tree = etree.iterparse(exposure_file, tag='%slocation' % nrml )
    for event, elem in tree: 
        loni = format(float(elem.attrib['lon']),'.5f')
        lon.append(format(float(loni),'.6e'))
        lati = format(float(elem.attrib['lat']),'.5f')
        lat.append(format(float(lati),'.6e'))
    tree = etree.iterparse(exposure_file, tag='%scost' % nrml )
    for event, elem in tree: 
        types.append(elem.attrib['type'])
        values.append(float(elem.attrib['value']))

    exposure_model = {'ids': None, 'number': None, 'taxonomies': None, 
        'lon': None, 'lat': None, 'types': None, 'values': None}   
    exposure_model['ids'] = ids
    exposure_model['number'] = number
    exposure_model['taxonomies'] = taxonomies
    exposure_model['lon'] = lon
    exposure_model['lat'] = lat
    exposure_model['types'] = types
    exposure_model['values'] = values
    
    return exposure_model


def read_fragility_model(fragility_file):
    damage_states = []
    taxonomies = []
    types = []
    IMT = []
    values = []
    nrml = '{http://openquake.org/xmlns/nrml/0.5}'
    
    tree = etree.iterparse(fragility_file, tag='%slimitStates' % nrml)
    for event, elem in tree:
        damage_states = str(elem.text).split()
        
    tree = etree.iterparse(fragility_file, tag='%sfragilityFunction' % nrml)
    for event, elem in tree:
        format = elem.attrib['format']
        types.append(format)
        taxonomies.append(elem.attrib['id'])
        if format == "continuous":
            iIMT,ivalues = parse_continuous_function(elem,nrml)
            IMT.append(iIMT)
            values.append(ivalues)
        if format == "discrete":
            iIMT,ivalues = parse_discrete_function(elem,nrml)
            IMT.append(iIMT)
            values.append(ivalues)
            
    fragility_model = {'damage_states': None, 'taxonomies': None, 'types': None, 
        'IMT': None, 'values': None}   
    fragility_model['damage_states'] = damage_states
    fragility_model['taxonomies'] = taxonomies
    fragility_model['types'] = types
    fragility_model['IMT'] = IMT
    fragility_model['values'] = values
    
    return fragility_model


def parse_continuous_function(elem,nrml):
    values = []
    for params in elem:
        if params.tag == '%simls' % nrml:
            IMT = params.attrib['imt']
        if params.tag == '%sparams' % nrml:
            mean = float(params.attrib['mean'])
            stddev = float(params.attrib['stddev'])
            ls = str(params.attrib['ls'])            
            values.append([mean,stddev,ls])

    return IMT, values 

###not working, ver como fazer o parse e por a legenda

def parse_discrete_function(elem,nrml):
    values = []
    for params in elem:
        if params.tag == '%simls' % nrml:
            IMT = params.attrib['imt']
            noDamageLimit = str(params.text).split()
            values.append(noDamageLimit)
            
        if params.tag == '%spoes' % nrml:
            poesValues = str(params.text).split()
            values.append(poesValues)            

    return IMT, values 


def read_vulnerability_model(vulnerability_file):
    taxonomies = []
    IMT = []
    values = []
    nrml = '{http://openquake.org/xmlns/nrml/0.5}'
    
    tree = etree.iterparse(vulnerability_file, tag='%svulnerabilityFunction' % nrml)
    for event, elem in tree:
        taxonomies.append(elem.attrib['id'])
        for params in elem:
            if params.tag == '%simls' % nrml:
                imt = params.attrib['imt']
                imls = str(params.text).split()
            if params.tag == '%smeanLRs' % nrml:
                meanLRs = str(params.text).split()
            if params.tag == '%scovLRs' % nrml:
                covLRs = str(params.text).split()     
        values.append(imls)
        values.append(meanLRs)
        values.append(covLRs)
        IMT.append(imt)    

    vulnerability_model = {'taxonomies': None, 'IMT': None, 'values': None}   
    vulnerability_model['taxonomies'] = taxonomies
    vulnerability_model['IMT'] = IMT
    vulnerability_model['values'] = values
    
    return vulnerability_model


def read_repair_model(repair_file):
#only works for discrete models
    from lxml import etree
    taxonomies = []
    model_units = []
    values = []
    damage_states = []
    types = []

    nrml = '{http://openquake.org/xmlns/nrml/0.5}'
    tree = etree.iterparse(repair_file, tag='%slimitStates' % nrml)
    for event, elem in tree:
        damage_states = str(elem.text).split()
        
    tree = etree.iterparse(repair_file, tag='%srepairFunction' % nrml )
    for event, elem in tree:
        format = elem.attrib['format']
        types.append(format)
        taxonomies.append(elem.attrib['id'])
        model_units.append(elem.attrib['model_units'])
        if format == "discrete":
            values_repair_model = []             
            for params in elem:
                if params.tag == '%sparams' % nrml:         
                    values_repair_model.append(params.attrib['time'])
            values.append(values_repair_model)

    repair_model = {'model_units': None, 'taxonomies': None, 
        'types': None, 'damage_states': None, 'values': None}   
    repair_model['model_units'] = model_units
    repair_model['taxonomies'] = taxonomies
    repair_model['types'] = types
    repair_model['damage_states'] = damage_states
    repair_model['values'] = values
        
    return repair_model


def calculate_damage(gmfs, exposure_model, fragility_model):
    fragility_model_median_beta = new_fragility_model_median_beta(fragility_model)
    lon = []
    lat = []
    taxonomies = []
    damageLevels = []     
    lon = exposure_model['lon']
    lat = exposure_model['lat']
    taxonomies = exposure_model['taxonomies']
    damage_states = fragility_model['damage_states']
    for i in xrange(len(lon)):
        if fragility_model['types'][fragility_model['taxonomies'].index(taxonomies[i])] == 'continuous':
            damageNode = calculate_damage_ff_continuous(i,lon[i],lat[i],taxonomies[i],gmfs[i],fragility_model_median_beta)
                    
        if fragility_model['types'][fragility_model['taxonomies'].index(taxonomies[i])] == 'discrete':
            damageNode = calculate_damage_ff_discrete(lon[i],lat[i],taxonomies[i],gmfs,fragility_model_median_beta)
        
        damageLevels.append(damageNode)
    
    damage_levels = {'lon': None, 'lat': None, 'taxonomies': None, 'damageLevels': None, 'damage_states':None}   
    damage_levels['lon'] = lon
    damage_levels['lat'] = lat
    damage_levels['damageLevels'] = damageLevels
    damage_levels['taxonomies'] = taxonomies
    damage_levels['damage_states'] = damage_states
        
    return damage_levels


def new_fragility_model_median_beta(fragility_model):
    indicesContinuousFF = [j for j, x in enumerate(fragility_model['types']) if x == 'continuous'] 
    fragility_model_median_beta = copy.deepcopy(fragility_model)
    for i in xrange(len(indicesContinuousFF)):
        valuesFF = fragility_model['values'][indicesContinuousFF[i]]
        newValues = []
        for j in xrange(len(fragility_model['damage_states'])):  
            median = valuesFF[j][0]**2/np.sqrt(valuesFF[j][1]**2+valuesFF[j][0]**2)
            beta = np.sqrt(np.log(valuesFF[j][1]**2/valuesFF[j][0]**2+1))
            newValues.append([round(median,4), round(beta,4)])
        
        fragility_model_median_beta['values'][indicesContinuousFF[i]] = newValues

    return fragility_model_median_beta


def calculate_damage_ff_continuous(i,lonNode,latNode,taxonomyNode,gmvNode,fragility_model_median_beta):
    damageNode = []
    ffindex = fragility_model_median_beta['taxonomies'].index(str(taxonomyNode))
    beta = fragility_model_median_beta['values'][ffindex][-1][1]
    median = fragility_model_median_beta['values'][ffindex][-1][0] 
    prob_Coll = stats.lognorm.cdf(gmvNode, beta, loc=0, scale=median)
    previous_prob = prob_Coll
    damage_probs_all = [prob_Coll]
        
    for j in xrange(len(fragility_model_median_beta['damage_states'])-2,-1,-1):
        median_j = fragility_model_median_beta['values'][ffindex][j][0]  
        beta_j = fragility_model_median_beta['values'][ffindex][j][1]  
        probDamage = stats.lognorm.cdf(gmvNode, beta_j, loc=0, scale=median_j) - previous_prob 
        previous_prob = stats.lognorm.cdf(gmvNode, beta_j, loc=0, scale=median_j)
        damage_probs_all.append(probDamage)   
      
    for i in xrange(len(gmvNode)):
        damageEachGmf = []
        for j in xrange(len(damage_probs_all)):
            damageEachGmf.append(damage_probs_all[j][i])
        damageNode.append(damageEachGmf)

    return damageNode


def calculate_damage_ff_discrete(lonNode,latNode,taxonomyNode,gmfs,fragility_model_median_beta):
####Not working yet####
    damageNode=[]   
    indicesLon = [j for j, x in enumerate(gmfs['lon']) if x == lonNode]
    indicesLat = [j for j, x in enumerate(gmfs['lat']) if x == latNode]
    indexGmfs = list(set(indicesLon) & set(indicesLat))
    gmvNode = gmfs['gmv'][int(str(indexGmfs[0]))]

    ffindex = fragility_model_median_beta['taxonomies'].index(str(taxonomyNode))
    
#ver o que da isto e como se calcula a partir daqui
    values_coll = fragility_model_median_beta['values'][ffindex][-1]

    for i in xrange(len(gmvNode)):
        x = gmvNode[i] 
        prob_Coll =  0

    return damageNode


def calculate_repair_time(repair_model, damage_levels):
    lon = []
    lat = []
    taxonomies = []
    repairTime = []     
    lon = damage_levels['lon']
    lat = damage_levels['lat']
    taxonomies = damage_levels['taxonomies']
    
    numDamageLevels = len(repair_model['values'][0])
    damageLevels = np.array(damage_levels['damageLevels'])
    repairTime = np.empty(damageLevels.shape)
    for i in xrange(len(lon)):
        repair_index = int(repair_model['taxonomies'].index(taxonomies[i]))   
        for x in xrange(numDamageLevels):
            repairTime[i,:,x] = damageLevels[i,:,x]*int(repair_model['values'][repair_index][numDamageLevels-1-x])

    repair_times = {'lon': None, 'lat': None, 'taxonomies': None, 'damageLevels': None}   
    repair_times['lon'] = lon
    repair_times['lat'] = lat
    repair_times['repairTime'] = repairTime
    repair_times['taxonomies'] = taxonomies
    
    return repair_times


def total_repair_time(repair_times):
    lon = []
    lat = []
    totalTimeByNode = []
    lon = repair_times['lon']
    lat = repair_times['lat']
    for i in xrange(len(lon)):
        jtotalTime = []
        for j in xrange(len(repair_times['repairTime'][0])):    
            total_repair_time_gmf = np.sum(repair_times['repairTime'][i][j])
            jtotalTime.append(total_repair_time_gmf)
        totalTimeByNode.append(jtotalTime)    

    return totalTimeByNode


def total_disruption_time(repair_times):
    lon = []
    lat = []
    totalDTimeByNode = []
    totalTime = []
    lon = repair_times['lon']
    lat = repair_times['lat']
    for i in xrange(len(lon)):
        jtotalTime = []
        for j in xrange(len(repair_times['repairTime'][0])):    
            total_disr_time_gmf = np.sum(repair_times['repairTime'][i][j][0]+repair_times['repairTime'][i][j][1])
            jtotalTime.append(total_disr_time_gmf)
        totalDTimeByNode.append(jtotalTime)
    
    return totalDTimeByNode


def factory_repair_time_vuln(vulnerability_model, gmfs):
    lon = []
    lat = []
    repairTime = []  
    taxonomies = []
    lon = gmfs['lon']
    lat = gmfs['lat']
    taxonomies = vulnerability_model['taxonomies']
    
    repairTime = np.interp(gmfs['gmv'][:], vulnerability_model['values'][0], vulnerability_model['values'][1])
 
    repair_times = {'lon': None, 'lat': None, 'taxonomies': None, 'repairTime': None}   
    repair_times['lon'] = lon
    repair_times['lat'] = lat
    repair_times['taxonomies'] = taxonomies
    repair_times['repairTime'] = repairTime

    return repair_times


def prob_level(damage_levels, prob_level):
    prob_coll = []
    prob_non_coll = []
    prob_coll_path = []
    prob_index = damage_levels['damage_states'].index(prob_level)
    for i in xrange(len(damage_levels['lon'])):
        prob_coll.append(list(zip(*damage_levels['damageLevels'][i])[prob_index]))
        prob_non_coll.append([1-x for x in prob_coll[i]]) 
    
    probNonCollprod = np.cumprod(prob_non_coll,axis=0)[-1]
    probCollPath = 1-probNonCollprod

    return probCollPath


def disr_prob_path(damage_levels):
    prob_disr = []
    prob_non_disr = []
    for i in xrange(len(damage_levels['lon'])):
        a = list(zip(*damage_levels['damageLevels'][i])[0])
        b = list(zip(*damage_levels['damageLevels'][i])[1])
        prob_disr.append([x + y for x, y in zip(a, b)])     
        prob_non_disr.append([1-x for x in prob_disr[i]])
       
    probNonDisrProd = np.cumprod(prob_non_disr,axis=0)[-1]
    prob_disr_path = 1-probNonDisrProd
  
    return prob_disr_path


def print_stats(value):
    mean = np.mean(value)
    median = np.median(value)
    stddev = np.std(value)
    print( 'Mean = %e' %mean)
    print( 'St Dev = %e' %stddev)
    print( 'Median = %e' %median)
    
    return mean, stddev, median



def apply_dam_RT_corr(damCorr, RTCorr, num_dam_samples, damage_levels, gmfs, fragility_model, exposure_model, timeShinozuka):
    
    #Build Damage correlation matrix
    damCorrUni = 2*np.sin(damCorr*np.pi/6)
    A = [1]
    first_row_dam = np.hstack((A,[damCorrUni]*(len(damage_levels['damageLevels'])-1)))
    corr_dam = toeplitz(first_row_dam, first_row_dam)
    chole_dam = np.linalg.cholesky(corr_dam)

    #Sample random numbers from a standard normal distribution
    mTot_dam = np.random.normal(0,1,num_dam_samples*len(gmfs['gmv'][0]))
    for i in xrange(len(damage_levels['damageLevels'])-1):
        mTot_dam = np.vstack([mTot_dam, np.random.normal(0,1,num_dam_samples*len(gmfs['gmv'][0]))])  
    #Convert the correlated samples to an uniform distribution
    mTotCorr_dam = np.dot(chole_dam,mTot_dam)
    mTotCorrUni_dam = 1./2*scipy.special.erf(-mTotCorr_dam/np.sqrt(2.))+0.5


    #Build Repair Time correlation matrix
    RTCorrUni = 2*np.sin(RTCorr*np.pi/6)
    first_row_RT = np.hstack((A,[RTCorrUni]*(len(damage_levels['damageLevels'])-1)))
    corr_RT = toeplitz(first_row_RT, first_row_RT)
    chole_RT = np.linalg.cholesky(corr_RT)

    #Sample random numbers from a standard normal distribution
    mTot_RT = np.random.normal(0,1,num_dam_samples*len(gmfs['gmv'][0]))
    for i in xrange(len(damage_levels['damageLevels'])-1):
        mTot_RT = np.vstack([mTot_RT, np.random.normal(0,1,num_dam_samples*len(gmfs['gmv'][0]))])
    #Convert the correlated samples to an uniform distribution   
    mTotCorr_RT = np.dot(chole_RT,mTot_RT)
    mTotCorrUni_RT = 1./2*scipy.special.erf(-mTotCorr_RT/np.sqrt(2.))+0.5   


    DamState_total = []
    RT_total = []
    DT_total = []
    
    for i in xrange(len(gmfs['lon'])):
        damDistrNode = []
        RTDistrNode = []
        DTDistrNode = []    
    
        for x in xrange(len(gmfs['gmv'][0])):
            pND = 1.-sum(damage_levels['damageLevels'][i][x])
            pS = damage_levels['damageLevels'][i][x][3]
            pM = damage_levels['damageLevels'][i][x][2]
            pE = damage_levels['damageLevels'][i][x][1]
            ffindex = fragility_model['taxonomies'].index(str(exposure_model['taxonomies'][i])) 
        
            for z in xrange(num_dam_samples):

                if mTotCorrUni_dam[i,x*num_dam_samples+z] < pND:
                    DSNode = 1
                elif mTotCorrUni_dam[i,x*num_dam_samples+z] < pND+pS:
                    DSNode = 2
                elif mTotCorrUni_dam[i,x*num_dam_samples+z] < pND+pS+pM:
                    DSNode = 3
                elif mTotCorrUni_dam[i,x*num_dam_samples+z] < pND+pS+pM+pE:
                    DSNode = 4
                else:
                    DSNode = 5
                damDistrNode.append(DSNode)

                RT = mTotCorrUni_RT[i][x*num_dam_samples+z]*(timeShinozuka[ffindex][DSNode-1][1]-timeShinozuka[ffindex][DSNode-1][0])+timeShinozuka[ffindex][DSNode-1][0]
                RTDistrNode.append(RT)
            
                if DSNode > 3:
                    DT = mTotCorrUni_RT[i][x*num_dam_samples+z]*(timeShinozuka[ffindex][DSNode-1][1]-timeShinozuka[ffindex][DSNode-1][0])+timeShinozuka[ffindex][DSNode-1][0]
                    DTDistrNode.append(DT)
                else:
                    DTDistrNode.append(0.0)
        
        DamState_total.append(damDistrNode)
        RT_total.append(RTDistrNode)
        DT_total.append(DTDistrNode)

    return DamState_total, RT_total, DT_total

def distance_on_sphere_numpy(coordinate_array):
    EARTH_RADIUS = 6371.0

    latitudes = coordinate_array[:, 0]
    longitudes = coordinate_array[:, 1]
    n_pts = coordinate_array.shape[0]

    # Convert latitude and longitude to spherical coordinates in radians.
    degrees_to_radians = np.pi/180.0
    phi_values = (90.0 - latitudes)*degrees_to_radians
    theta_values = longitudes*degrees_to_radians

    # Expand phi_values and theta_values into grids
    theta_1, theta_2 = np.meshgrid(theta_values, theta_values)
    theta_diff_mat = theta_1 - theta_2

    phi_1, phi_2 = np.meshgrid(phi_values, phi_values)

    # Compute spherical distance from spherical coordinates
    angle = (np.sin(phi_1) * np.sin(phi_2) * np.cos(theta_diff_mat) + 
           np.cos(phi_1) * np.cos(phi_2))
    arc = np.arccos(angle)

    # Multiply by earth's radius to obtain distance in km
    return arc * EARTH_RADIUS
  
