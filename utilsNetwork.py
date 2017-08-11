import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import shapely.geometry as shgeom
import fiona
from heapq import heappush, heappop
from itertools import count


def haversine(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = np.radians([lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    km = 6367* 2 * np.arcsin(np.sqrt(a)) 
    return km 


def shp_adj(networkShp):
    G = nx.read_shp(networkShp)
    mat = nx.adjacency_matrix(G)
    maxNumConn = max((mat).sum(1))
    connectionMat = np.zeros(shape=(mat.shape[0],maxNumConn))
    for i in xrange(mat.shape[0]):
        noZero = np.nonzero(mat[i])
        z = np.asarray(noZero[1])
        connectionMat[i, 0:z.shape[0]] = (z[0:z.shape[0]]+1)

    shpAdj = np.zeros(shape=(len(G.nodes()),maxNumConn*2+5))
    Gnodes = np.array(G.nodes())
    shpAdj[:,0] = range(1,len(G.nodes())+1)
    shpAdj[:,1] = Gnodes[:,0]
    shpAdj[:,2] = Gnodes[:,1]
    for i in xrange(maxNumConn):
        shpAdj[:,3+i] = connectionMat[:,i]
        nonZero = np.nonzero(shpAdj[:,3+i])
        shpAdj[nonZero[0],5+i+maxNumConn] = haversine(shpAdj[nonZero[0],1],shpAdj[nonZero[0],2],shpAdj[(shpAdj[nonZero[0],3+i]-1).astype(int),1],shpAdj[(shpAdj[nonZero[0],3+i]-1).astype(int),2])

    #np.savetxt(resultsFolder+'network-shpAdj.txt', shpAdj, fmt='%i %f %f %i %i %i %i %i %i %i %f %f %f %f %f', delimiter=',')
    return shpAdj,maxNumConn


def divide_edges(shpAdj,maxNumConn,maxDist):
    shpAdjMaxDist = np.zeros(shape=(shpAdj.shape[0],shpAdj.shape[1]))
    shpAdjMaxDist = shpAdj
    linesToDivide = np.where(shpAdjMaxDist[:,5+maxNumConn]>maxDist)
    while len(linesToDivide[0])!=0:
        for i in range(len(linesToDivide[0])):
            newNode = len(shpAdjMaxDist)
            newLine = np.zeros(shape=(1,len(shpAdjMaxDist[0])))
            newLine[0,0] = len(shpAdjMaxDist)+1
            newLine[0,1] = (shpAdjMaxDist[linesToDivide[0][i],1]+shpAdjMaxDist[shpAdjMaxDist[linesToDivide[0][i],3]-1,1])/2
            newLine[0,2] = (shpAdjMaxDist[linesToDivide[0][i],2]+shpAdjMaxDist[shpAdjMaxDist[linesToDivide[0][i],3]-1,2])/2
            newLine[0,3] = shpAdjMaxDist[linesToDivide[0][i],3]
            newLine[0,5+maxNumConn] = shpAdjMaxDist[linesToDivide[0][i],10]/2
            shpAdjMaxDist[linesToDivide[0][i],3] = len(shpAdjMaxDist)+1
            shpAdjMaxDist[linesToDivide[0][i],5+maxNumConn] = newLine[0,5+maxNumConn]
            shpAdjMaxDist = np.vstack((shpAdjMaxDist,newLine))
        linesToDivide = np.where(shpAdjMaxDist[:,5+maxNumConn]>maxDist)

    linesToDivide2 = np.where(shpAdjMaxDist[:,6+maxNumConn]>maxDist)
    while len(linesToDivide2[0])!=0:
        for i in range(len(linesToDivide2[0])):
            newNode = len(shpAdjMaxDist)
            newLine = np.zeros(shape=(1,len(shpAdjMaxDist[0])))
            newLine[0,0] = len(shpAdjMaxDist)+1
            newLine[0,1] = (shpAdjMaxDist[linesToDivide2[0][i],1]+shpAdjMaxDist[shpAdjMaxDist[linesToDivide2[0][i],4]-1,1])/2
            newLine[0,2] = (shpAdjMaxDist[linesToDivide2[0][i],2]+shpAdjMaxDist[shpAdjMaxDist[linesToDivide2[0][i],4]-1,2])/2
            newLine[0,3] = shpAdjMaxDist[linesToDivide2[0][i],4]
            newLine[0,5+maxNumConn] = shpAdjMaxDist[linesToDivide2[0][i],11]/2
            shpAdjMaxDist[linesToDivide2[0][i],4] = len(shpAdjMaxDist)+1
            shpAdjMaxDist[linesToDivide2[0][i],6+maxNumConn] = newLine[0,5+maxNumConn]
            shpAdjMaxDist = np.vstack((shpAdjMaxDist,newLine))
        linesToDivide2 = np.where(shpAdjMaxDist[:,6+maxNumConn]>maxDist)

    #np.savetxt(resultsFolder+'network-shpAdjMaxDist.txt', shpAdjMaxDist, fmt='%i %f %f %i %i %i %i %i %i %i %f %f %f %f %f', delimiter=',')
    return shpAdjMaxDist


def shp_adj_br(bridgesShp,maxDistBr,resultsFolder):
    G = nx.read_shp(bridgesShp)
    mat = nx.adjacency_matrix(G)
    maxNumConn = max((mat).sum(1))
    connectionMat = np.zeros(shape=(mat.shape[0],maxNumConn))
    for i in xrange(mat.shape[0]):
        noZero = np.nonzero(mat[i])
        z = np.asarray(noZero[1])
        connectionMat[i, 0:z.shape[0]] = (z[0:z.shape[0]]+1)

    shpAdjBr = np.zeros(shape=(len(G.nodes()),maxNumConn*2+5))
    Gnodes = np.array(G.nodes())
    shpAdjBr[:,0] = range(1,len(G.nodes())+1)
    shpAdjBr[:,1] = Gnodes[:,0]
    shpAdjBr[:,2] = Gnodes[:,1]
    for i in xrange(maxNumConn):
        shpAdjBr[:,3+i] = connectionMat[:,i]
        nonZero = np.nonzero(shpAdjBr[:,3+i])
        shpAdjBr[nonZero[0],5+i+maxNumConn] = haversine(shpAdjBr[nonZero[0],1],shpAdjBr[nonZero[0],2],shpAdjBr[(shpAdjBr[nonZero[0],3+i]-1).astype(int),1],shpAdjBr[(shpAdjBr[nonZero[0],3+i]-1).astype(int),2])

    numConn = np.array(shpAdjBr[:,3:3+maxNumConn]).tolist()
    numConnReshape = (np.reshape(numConn,maxNumConn*len(numConn))).astype(int)

    aab = np.zeros(shape=(len(shpAdjBr),10))
    aab[:,0] = np.bincount(numConnReshape)[1:len(numConnReshape-1)]-1 
    aab[:,0] = aab[:,0].clip(0,1) 
    aac = np.nonzero(aab[:,0])[0]+1

    nodeBifurcation = np.nonzero(shpAdjBr[:,4])[0]+1
    nodeExit = np.where(shpAdjBr[:,3]==0)[0]+1
    nodes3Cases = np.concatenate((aac,nodeBifurcation,nodeExit),axis=0)
    nodesNoDelete = np.unique(nodes3Cases)
    for i in xrange(len(nodesNoDelete)):
        la=np.where(shpAdjBr[:,3:3+maxNumConn] == nodesNoDelete[i])
        aab[i,1:1+len(la[0])]=la[0]+1

    aao = np.unique(aab[:,1:])-1
    shpAdjBr[aao[1:].astype(int),maxNumConn+3]=1
    for j in xrange(len(shpAdjBr)):
        sd = np.nonzero(shpAdjBr[j,3:3+maxNumConn+1+1])
        if len(sd[0]) == 1:
            while len(sd[0]) == 1:
                mm = shpAdjBr[j,3]-1
#Condition to limit distance between points
#(if condition not needed delete next line and the "else" part)
                if shpAdjBr[j,5+maxNumConn]+shpAdjBr[mm,5+maxNumConn]<maxDistBr:
                    shpAdjBr[j,3:3+maxNumConn+1+1] = shpAdjBr[mm,3:3+maxNumConn+1+1]
                    shpAdjBr[j,5+maxNumConn] = shpAdjBr[j,5+maxNumConn]+shpAdjBr[mm,5+maxNumConn]
                    shpAdjBr[j,6+maxNumConn:6+maxNumConn*2-1] = shpAdjBr[mm,6+maxNumConn:6+maxNumConn*2-1]
                    shpAdjBr[mm,:] = 0
                    sd = np.nonzero(shpAdjBr[j,3:3+maxNumConn+1+1])
                else: #This else serves only to stop the while if dist>maxDistBr
                    sd = [(1, 1)]

    bbb = np.nonzero(shpAdjBr[:,0])
    br_adj = shpAdjBr[bbb,0:3+maxNumConn*2+3]

    np.savetxt(resultsFolder+'network-br_adj.txt', br_adj[0], fmt='%i %f %f %i %i %i %i %i %i', delimiter=',')
    return br_adj[0]


def one_node_per_bridge(br_adj,shpAdjMaxDist,resultsFolder):
    p = np.where(br_adj[:,3]==0)
    b = np.delete(br_adj,p,0)
    maxNumConn = (len(br_adj[1])-5)/2
    br_XY = b[:,[1,2,5+maxNumConn]]

    brRows = []
    for i in xrange(len(br_XY)):
        ab = np.logical_and(shpAdjMaxDist[:,1]==br_XY[i,0], shpAdjMaxDist[:,2]==br_XY[i,1]).nonzero()[0]
        brRows.append([ab[0],br_XY[i,2]])
    
    np.savetxt(resultsFolder+'br_length.txt', br_XY, fmt='%f %f %f', delimiter=',')

    return np.asarray(brRows)


def find_br_rows_closest_coord(shpAdjMaxDist,brNodes):
    brRows = []
    NDIM = 2
    brNodes = np.loadtxt(brNodes, delimiter=",", skiprows=1, usecols=(0,1), unpack=False)
    nodes = shpAdjMaxDist[:,1:3]
    nodes.shape = nodes.size / NDIM, NDIM
    nodesrad = np.radians(nodes)
    for i in xrange(len(brNodes)):
        point =  brNodes[i,:]
        d = ((nodes-point)**2).sum(axis=1)
        ndx = d.argsort()
        brRows.append([ndx[0],0])

    return np.asarray(brRows)


def find_br_rows_exact_coord(shpAdjMaxDist,brNodes):
    brNodes = np.loadtxt(brNodes, delimiter=",", skiprows=1, usecols=(0,1), unpack=False)
    brNodes1 = np.zeros(shape=(len(brNodes),2))
    brNodes1[:,0:2] = np.round(brNodes[:,0:2],6)
    shpAdjMaxDist[:,1] = np.round(shpAdjMaxDist[:,1],6)
    shpAdjMaxDist[:,2] = np.round(shpAdjMaxDist[:,2],6)
    brRowsArray = np.zeros(shape=(len(brNodes1),1))
    for i in xrange(len(brNodes)):
        ab = np.logical_and(shpAdjMaxDist[:,1]==brNodes1[i,0], shpAdjMaxDist[:,2]==brNodes1[i,1]).nonzero()[0]
        if len(ab) !=0:
            brRowsArray[i,0] = ab[0]

    p1 = np.where(brRowsArray[:,0]==0)
    brRowsArray = np.delete(brRowsArray, p1, 0)
    brRows = []
    for i in xrange(len(brRowsArray)):
        brRows.append([brRowsArray[i][0],0])

    return np.asarray(brRows)


def sim_adj(shpAdjMaxDist,brRows,brPerpRows,maxDist,limit_length1,limit_length2):
    maxNumConn = (len(shpAdjMaxDist[1])-5)/2
    numConn = np.array(shpAdjMaxDist[:,3:3+maxNumConn]).tolist()
    numConnReshape = (np.reshape(numConn,maxNumConn*len(numConn))).astype(int)
    aab = np.zeros(shape=(len(shpAdjMaxDist),10))
    aab[:,0] = np.bincount(numConnReshape)[1:len(numConnReshape-1)]-1 
    aab[:,0] = aab[:,0].clip(0,1) 
    aac = np.nonzero(aab[:,0])[0]+1

    nodeBifurcation = np.nonzero(shpAdjMaxDist[:,4])[0]+1
    nodeExit = np.where(shpAdjMaxDist[:,3] == 0)[0]+1
    nodes3Cases = np.concatenate((aac,nodeBifurcation,nodeExit),axis=0)
    nodesNoDelete = np.unique(nodes3Cases)
    for i in xrange(len(nodesNoDelete)):
        la = np.where(shpAdjMaxDist[:,3:3+maxNumConn] == nodesNoDelete[i])
        aab[i,1:1+len(la[0])] = la[0]+1

    aao = np.unique(aab[:,1:])-1
    shpAdjMaxDist[aao[1:].astype(int),maxNumConn+3] = 1
    
    #calcular dist para cada ponte e por em maxNumConn+4 a FF correspondente
    for i in xrange(len(brRows)):
        if brRows[i,1] < limit_length1/1000.:
            shpAdjMaxDist[brRows[i,0],maxNumConn+4] = 1
        elif brRows[i,1] < limit_length2/1000.:
            shpAdjMaxDist[brRows[i,0],maxNumConn+4] = 2
        else:
			shpAdjMaxDist[brRows[i,0],maxNumConn+4] = 3
    if len(brPerpRows) != 0:
		shpAdjMaxDist[brPerpRows[:,0],maxNumConn+4] = 4
    aaa1 = shpAdjMaxDist
    for j in xrange(len(shpAdjMaxDist)):
        sd = np.nonzero(shpAdjMaxDist[j,3:3+maxNumConn])
        if (len(sd[0]) == 1 and shpAdjMaxDist[j,maxNumConn+3] == 0 and shpAdjMaxDist[j,maxNumConn+4] == 0):
            while len(sd[0]) == 1:
                mm = shpAdjMaxDist[j,3]-1
#Condition to limit distance between points 
#(if condition not needed delete next line and the "else" part)
                if shpAdjMaxDist[j,5+maxNumConn]+shpAdjMaxDist[mm,5+maxNumConn]<maxDist:
                    shpAdjMaxDist[j,3:5+maxNumConn] = shpAdjMaxDist[mm,3:5+maxNumConn]
                    shpAdjMaxDist[j,5+maxNumConn] = shpAdjMaxDist[j,5+maxNumConn]+shpAdjMaxDist[mm,5+maxNumConn]
                    shpAdjMaxDist[j,6+maxNumConn:6+maxNumConn*2-1] = shpAdjMaxDist[mm,6+maxNumConn:6+maxNumConn*2-1]
                    shpAdjMaxDist[mm,:] = 0
                    sd = np.nonzero(shpAdjMaxDist[j,3:3+maxNumConn+1+2])
                else: #This else serves only to stop the while if dist>maxDist
                    sd = [(1, 1)]

    bbb = np.nonzero(shpAdjMaxDist[:,0])
    adj = shpAdjMaxDist[bbb,0:3+maxNumConn*2+3]

    #np.savetxt(resultsFolder+'network-adj.txt', adj[0], fmt='%i %f %f %i %i %i %i %i %i %i %f %f %f %f %f', delimiter=',')
    return adj[0]


def save_files_networkx(adj,maxNumConn,resultsFolder):
    #adj = np.loadtxt(adj, delimiter=" ", unpack=False)
    nodes = adj[:,0:3]	
    np.savetxt(resultsFolder+'network-nodes.txt', adj[:,0:3], fmt='%i %f %f', delimiter=',')
    matadj = np.zeros(shape=(len(adj),maxNumConn*3))
    matadj[:,0] = adj[:,0]
    matadj[:,1] = adj[:,3]
    matadj[:,2] = adj[:,5+float(maxNumConn)]
    if maxNumConn >1:
        matadj[:,3] = adj[:,0]
        matadj[:,4] = adj[:,4]
        matadj[:,5] = adj[:,6+float(maxNumConn)]
    if maxNumConn >2:
        matadj[:,6] = adj[:,0]
        matadj[:,7] = adj[:,5]
        matadj[:,8] = adj[:,7+float(maxNumConn)]
    if maxNumConn >3:	
        matadj[:,9] = adj[:,0]
        matadj[:,10] = adj[:,6]
        matadj[:,11] = adj[:,8+float(maxNumConn)]
    if maxNumConn >4:
        matadj[:,12] = adj[:,0]
        matadj[:,13] = adj[:,7]
        matadj[:,14] = adj[:,9+float(maxNumConn)]

    matadjAll = np.reshape(matadj,(len(adj)*maxNumConn,3))
    matadjZeros = np.where(matadjAll[:,1]==0)
    matadj = np.delete(matadjAll, matadjZeros, 0)
    edges = matadj[:,0:2]
    weights = matadj[:,0:3]

    np.savetxt(resultsFolder+'network-edges.txt', matadj[:,0:2], fmt='%i %i', delimiter=',')
    np.savetxt(resultsFolder+'network-weights.txt', matadj[:,0:3], fmt='%i %i %f', delimiter=',')
    return nodes, edges, weights


def node_in_out(adj,lon1,lat1,lon2,lat2):
    relations = np.reshape(adj[:,3:8],(len(adj)*5,1))
    relations = np.unique(relations)

    nodesIN = np.zeros(shape=(len(adj),1))
    for i in xrange(len(adj)):
        if adj[i,0] not in relations:
            nodesIN[i,0] = adj[i,0]

    nodesIN = np.unique(nodesIN)[1:,]
    distIN = np.zeros(shape=(len(nodesIN),3))
    distIN[:,0] = nodesIN
    for i in xrange(len(distIN)):
        a = np.where(adj[:,0]==distIN[i,0])
        distIN[i,1] = haversine(lon1,lat1,adj[a[0][0],1],adj[a[0][0],2])

    rowNodeIN = np.where(min(distIN[:,1])==distIN[:,1])[0][0]
    nodeIN = distIN[rowNodeIN,0]

    rowNodesOut = np.where(adj[:,3]==0)
    distOUT = np.zeros(shape=(len(rowNodesOut[0]),3))
    distOUT[:,0] = adj[rowNodesOut[0],0]
    for i in xrange(len(distOUT)):
        distOUT[i,1] = haversine(lon2,lat2,adj[rowNodesOut[0][i],1],adj[rowNodesOut[0][i],2])

    rowNodeOUT = np.where(min(distOUT[:,1])==distOUT[:,1])[0][0]
    nodeOUT = distOUT[rowNodeOUT,0]
    return nodeIN,nodeOUT


def draw_network(nodes,edges,weights,G):
    #import nodes
    for i in xrange(len(nodes)):
        G.add_node(nodes[i,0], pos = (nodes[i,1], nodes[i,2]))
    pos = nx.get_node_attributes(G,'pos')
    #import edges
    for i in xrange(len(edges)):
        G.add_edge(edges[i,0], edges[i,1])
    #import weights (distances)
    for i in xrange(len(weights)):
        G.add_weighted_edges_from([(weights[i,0],weights[i,1],weights[i,2])])
    #draw the network
    plt.figure(figsize=(13,7))
    nx.draw(G,pos,node_size=20,marker='o',linewidths = 0.4)
    lon1 = min(nodes[:,1])
    lon2 = max(nodes[:,1])
    minLon = np.min([lon1,lon2])-0.2
    maxLon = np.max([lon1,lon2])+0.2
    lat1 = min(nodes[:,2])
    lat2 = max(nodes[:,2])
    minLat = np.min([lat1,lat2])-0.2
    maxLat = np.max([lat1,lat2])+0.2
    pngLimits = [minLon,minLat,maxLon,maxLat]
    map = Basemap(resolution = 'i', area_thresh = 1000.0,
            llcrnrlon=minLon, llcrnrlat=minLat, urcrnrlon=maxLon, urcrnrlat=maxLat)
    map.drawcountries(linewidth=1)
    map.drawmapboundary(linewidth=1)
    map.fillcontinents(color='coral',alpha=0.3)
    parallels = np.arange(int(minLat),int(maxLat)+1,1.)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    meridians = np.arange(int(minLon),int(maxLon)+1,1.)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    
    return plt.show(),G,pos,pngLimits


def draw_nodes(source,destination,G,pos,pngLimits):
    plt.figure(figsize=(13,7))
    map = Basemap(resolution = 'i', area_thresh = 1000.0,
        llcrnrlon=pngLimits[0], llcrnrlat=pngLimits[1], urcrnrlon=pngLimits[2], urcrnrlat=pngLimits[3])
    map.drawcountries(linewidth=1)
    map.fillcontinents(color='coral',alpha=0.3)
    map.drawmapboundary(linewidth=1)
    parallels = np.arange(int(pngLimits[1]),int(pngLimits[3])+1,1.)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    meridians = np.arange(int(pngLimits[0]),int(pngLimits[2])+1,1.)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    nx.draw(G,pos,node_size=20,marker='o',linewidths = 0.5)
    nx.draw_networkx_nodes(G, pos, nodelist=[source,destination],node_size=50, node_color='b') 
    return plt.show()


def shortest_path(G,source,destination,pos,filename,saveImage,resultsFolder,pngLimits):
    plt.figure(figsize=(13,7))
    map = Basemap(resolution = 'i', area_thresh = 1000.0,
        llcrnrlon=pngLimits[0], llcrnrlat=pngLimits[1], urcrnrlon=pngLimits[2], urcrnrlat=pngLimits[3])
    map.drawcountries(linewidth=1)
    map.fillcontinents(color='coral',alpha=0.3)
    map.drawmapboundary(linewidth=1)
    parallels = np.arange(int(pngLimits[1]),int(pngLimits[3])+1,1.)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    meridians = np.arange(int(pngLimits[0]),int(pngLimits[2])+1,1.)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)

    path1 = (nx.dijkstra_path(G,source,destination,weight='weights'))
    nx.draw(G,pos,node_size=20,marker='o',linewidths = 0.5)
    nx.draw_networkx_nodes(G, pos, nodelist=path1, node_size=20, node_color='b') 
    np.savetxt(resultsFolder+'path1.txt', path1, fmt='%i', delimiter=',')
    if saveImage:
        plt.savefig(filename, format='png')

    return plt.show()


def paths_1to5(paths2,G,pos,resultsFolder,pngLimits):
    paths2nodes = np.asarray(paths2[1]).transpose()
    paths2nodesl = paths2nodes.tolist()
    #np.savetxt(resultsFolder+'paths2dist.txt', paths2[0], fmt='%f', delimiter=',')
    np.savetxt(resultsFolder+'path1.txt', paths2nodes[0], fmt='%i', delimiter=',')
    np.savetxt(resultsFolder+'path2.txt', paths2nodes[1], fmt='%i', delimiter=',')
    np.savetxt(resultsFolder+'path3.txt', paths2nodes[2], fmt='%i', delimiter=',')
    np.savetxt(resultsFolder+'path4.txt', paths2nodes[3], fmt='%i', delimiter=',')
    np.savetxt(resultsFolder+'path5.txt', paths2nodes[4], fmt='%i', delimiter=',')

    plt.figure(figsize=(15,8))
    map = Basemap(resolution = 'i', area_thresh = 1000.0,
        llcrnrlon=pngLimits[0], llcrnrlat=pngLimits[1], urcrnrlon=pngLimits[2], urcrnrlat=pngLimits[3])
    map.drawcountries(linewidth=1)
    map.fillcontinents(color='coral',alpha=0.3)
    map.drawmapboundary(linewidth=1)
    parallels = np.arange(int(pngLimits[1]),int(pngLimits[3])+1,1.)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    meridians = np.arange(int(pngLimits[0]),int(pngLimits[2])+1,1.)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    
    nx.draw(G,pos,node_size=20,marker='o',linewidths = 0.5)
    #nx.draw_networkx_nodes(G, pos, nodelist=paths2nodes[0],node_size=20, 	node_color='g') 
    nx.draw_networkx_nodes(G, pos, nodelist=paths2nodes[1],node_size=20, node_color='y') 
    #nx.draw_networkx_nodes(G, pos, nodelist=paths2nodes[2],node_size=20, node_color='r') 
    #nx.draw_networkx_nodes(G, pos, nodelist=paths2nodes[3],node_size=20, node_color='b') 
    #nx.draw_networkx_nodes(G, pos, nodelist=paths2nodes[4],node_size=20, node_color='c') 
    return plt.show()


def RSA_zones(zoneA,zoneB,zoneC,zoneD,pathNodes,resultsFolder):
    zoneAopen = fiona.open(zoneA)
    polygonA = shgeom.shape(zoneAopen.next()['geometry'])
    zoneBopen = fiona.open(zoneB)
    polygonB = shgeom.shape(zoneBopen.next()['geometry'])
    zoneCopen = fiona.open(zoneC)
    polygonC = shgeom.shape(zoneCopen.next()['geometry'])
    zoneDopen = fiona.open(zoneD)
    polygonD = shgeom.shape(zoneDopen.next()['geometry'])
    #This aditional polygon is to include the bridges that cross Tejo,
    #because the shp does not include that part of the teritory
    polygonLx = shgeom.Polygon([(-9.2, 38.6), (-9.2, 38.8), (-9, 38.8), (-9, 38.6)])

    nodes = np.loadtxt(pathNodes, delimiter=" ", unpack=False)
    nodesZones = []
    taxonomy = []
    for i in xrange(len(nodes)):
        if nodes[i,2]!=0:
            X = nodes[i,1]
            Y = nodes[i,2]
            point = shgeom.Point(X, Y)
            if polygonA.contains(point):
                if nodes[i,3]==1:
                    taxonomy.append('1A')
                if nodes[i,3]==2:
                    taxonomy.append('2A')
            elif polygonLx.contains(point):
                if nodes[i,3]==1:
                    taxonomy.append('1A')
                if nodes[i,3]==2:
                    taxonomy.append('2A')
            elif polygonB.contains(point):
                if nodes[i,3]==1:
                    taxonomy.append('1B')
                if nodes[i,3]==2:
                    taxonomy.append('2B')
            elif polygonC.contains(point):
                if nodes[i,3]==1:
                    taxonomy.append('1C')
                if nodes[i,3]==2:
                    taxonomy.append('2C')
            elif polygonD.contains(point):
                if nodes[i,3]==1:
                    taxonomy.append('1D')
                if nodes[i,3]==2:
                    taxonomy.append('2D')
            else:
                print('Node not found')
                print(i)

    nodesZones = np.zeros(len(nodes), dtype=[('var1',float),('var2',float),('var3','a5')])
    nodesZones['var1'] = nodes[:,1]
    nodesZones['var2'] = nodes[:,2]
    nodesZones['var3'] = taxonomy
        
    np.savetxt(resultsFolder+'network-ffZones.txt', nodesZones, fmt="%f %f %s", delimiter=',')


def path_coord(path,nodes,path_number):
    #Extract path coordinates # Repeat for all the paths to calculate in OQ
    path1 = np.loadtxt(path, delimiter=" ", unpack=False)
    nodes = np.loadtxt(nodes, delimiter=" ", unpack=False)
    nodes3 = np.zeros(shape=(len(path1),3))
    for i in xrange(len(path1)):
        b = np.where(nodes[:,0]==path1[i])
        nodes3[i,:] = nodes[b[0][0],:]
    #np.savetxt(resultsFolder+'path2-'+str(path_number)+'-coord.txt', nodes3, fmt='%i %.6f %.6f', delimiter=',')


def path_br_ff(path,adj,path_number,resultsFolder):
    path1 = np.loadtxt(path, delimiter=" ", unpack=False)
    maxNumConn = (len(adj[1])-5)/2
    path1Br = np.zeros(shape=(len(path1),4))
    for i in xrange(len(path1)):
        b = np.where(adj[:,0]==path1[i])
        if adj[b[0][0],4+maxNumConn]!=0:
            path1Br[i,0:3] = adj[b[0][0],0:3]
            path1Br[i,3] = adj[b[0][0],4+maxNumConn]

    BrZeros = np.where(path1Br[:,1]==0)
    path1Br = np.delete(path1Br, BrZeros, 0)

    np.savetxt(resultsFolder+'path'+str(path_number)+'-exposure.txt', path1Br, fmt='%i %.6f %.6f %i', delimiter=',')
    return path1Br


def path_ff(path,adj,path_number,resultsFolder):
    path1 = np.loadtxt(path, delimiter=" ", unpack=False)
    #adj = np.loadtxt(adj, delimiter=" ", unpack=False)
    maxNumConn = (len(adj[1])-5)/2
    path1Br = np.zeros(shape=(len(path1),4))
    for i in xrange(len(path1)):
        b = np.where(adj[:,0]==path1[i])
        path1Br[i,0:3] = adj[b[0][0],0:3]
        path1Br[i,3] = adj[b[0][0],4+maxNumConn]

    BrZeros = np.where(path1Br[:,1]==0)
    path1Br = np.delete(path1Br, BrZeros, 0)

    np.savetxt(resultsFolder+'path'+str(path_number)+'-exposure.txt', path1Br, fmt='%i %.6f %.6f %i', delimiter=',')
    return path1Br


def exposure_1path(exposureTxt,resultsFolder):
    nodesFF = np.loadtxt(exposureTxt, delimiter=" ", unpack=False)
    np.savetxt(resultsFolder+'path-expCsv.csv', nodesFF[:,1:3], fmt='%f,%f', delimiter=',')


def exposure_n_paths(nodesAll,resultsFolder):
    nodesFF = unique(nodesAll)
    np.savetxt(resultsFolder+'paths-expCsv.csv', nodesFF[:,1:3], fmt='%f,%f', delimiter=',')


def exposure_only_br(nodes,adj,resultsFolder):
    #nodes = np.loadtxt(nodes, delimiter=" ", unpack=False)
    #adj = np.loadtxt(adj, delimiter=" ", unpack=False)
    maxNumConn = (len(adj[1])-5)/2
    nodesBrFF = np.zeros(shape=(len(nodes),4))
    for i in xrange(len(nodes)):
        b = np.where(adj[:,0]==nodes[i,0])
        if adj[b[0][0],4+maxNumConn]!=0:
            nodesBrFF[i,0:3] = adj[b[0][0],0:3]
            nodesBrFF[i,3] = adj[b[0][0],4+maxNumConn]

    BrZeros = np.where(nodesBrFF[:,1]==0)
    nodesBrFF = np.delete(nodesBrFF, BrZeros, 0)

    np.savetxt(resultsFolder+'network-ff.txt', nodesBrFF, fmt='%i %.6f %.6f %i', delimiter=',')
    np.savetxt(resultsFolder+'country-expCsv.csv', nodesBrFF[:,1:3], fmt='%f,%f', delimiter=',')


def exposure_br_pav(nodes,adj,resultsFolder):
    #nodes = np.loadtxt(nodes, delimiter=" ", unpack=False)
    #adj = np.loadtxt(adj, delimiter=" ", unpack=False)
    maxNumConn = (len(adj[1])-5)/2
    np.savetxt(resultsFolder+'locationscountry.csv', nodes[:,1:3], fmt='%f,%f', delimiter=',')
    nodesFF = np.zeros(shape=(len(adj),3))
    nodesFF[:,0] = adj[:,1]
    nodesFF[:,1] = adj[:,2]
    nodesFF[:,2] = adj[:,maxNumConn+4]
    np.savetxt(resultsFolder+'network-ff.txt', nodesFF, fmt='%.6f %.6f %i', delimiter=',')


def k_shortest_paths(G, source, destination, k, weight='weight'):
    if source == destination:
        return ([0], [[source]]) 
       
    length, path = nx.single_source_dijkstra(G, source, destination, weight=weight)
    if destination not in length:
        raise nx.NetworkXNoPath("node %s not reachable from %s" % (source, destination))
        
    lengths = [length[destination]]
    paths = [path[destination]]
    c = count()        
    B = []                        
    G_original = G.copy()    
    for i in range(1, k):
        for j in range(len(paths[-1]) - 1):            
            spur_node = paths[-1][j]
            root_path = paths[-1][:j + 1]
            
            edges_removed = []
            for c_path in paths:
                if len(c_path) > j and root_path == c_path[:j + 1]:
                    u = c_path[j]
                    v = c_path[j + 1]
                    if G.has_edge(u, v):
                        edge_attr = G.edge[u][v]
                        G.remove_edge(u, v)
                        edges_removed.append((u, v, edge_attr))
            
            for n in range(len(root_path) - 1):
                node = root_path[n]
                # out-edges
                for u, v, edge_attr in G.edges_iter(node, data=True):
                    G.remove_edge(u, v)
                    edges_removed.append((u, v, edge_attr))
                
                if G.is_directed():
                    # in-edges
                    for u, v, edge_attr in G.in_edges_iter(node, data=True):
                        G.remove_edge(u, v)
                        edges_removed.append((u, v, edge_attr))
            
            spur_path_length, spur_path = nx.single_source_dijkstra(G, spur_node, destination, weight=weight)            
            if destination in spur_path and spur_path[destination]:
                total_path = root_path[:-1] + spur_path[destination]
                total_path_length = get_path_length(G_original, root_path, weight) + spur_path_length[destination]                
                heappush(B, (total_path_length, next(c), total_path))
                
            for e in edges_removed:
                u, v, edge_attr = e
                G.add_edge(u, v, edge_attr)
                       
        if B:
            (l, _, p) = heappop(B)        
            lengths.append(l)
            paths.append(p)
        else:
            break
    
    return (lengths, paths)


def get_path_length(G, path, weight='weight'):
    length = 0
    if len(path) > 1:
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]
            length += G.edge[u][v].get(weight, 1)
    
    return length     
