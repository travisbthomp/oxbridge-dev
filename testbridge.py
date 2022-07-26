#=============================================================================
# ** Oxford Mathematical Brain Modeling Group **
#   Test script for the Oxbridge Graph Reader / Writer software for bridging 
#   various OxMBM utilities 
#
#   -----------------------------------------------------------
#   Authors
#   -----------------------------------------------------------
#
#       Travis Thompson             thompsont@maths.ox.ac.uk
#                      ----
#       Prama Putra                 putra@maths.ox.ac.uk
#                      ----
#       Pavanjit Chaggar            chaggar@maths.ox.ac.uk
#                      ----
#       Georgia S. Brennan          brennan@maths.ox.ac.uk
#                      ----
#       Andrew O'Heachteirn         oheachteirn@maths.ox.ac.uk
#                      ----
#       Christoffer Alexandersen    christoffer.alexandersen@maths.ox.ac.uk
#                      ----
#       Hadrien Oliveri             oliveri@maths.ox.ac.uk
#                      ----
#       Heather Harrington          harrington@maths.ox.ac.uk
#                      ----
#       Alain Goriely               goriely@maths.ox.ac.uk
#
# ============================================================================

from oxbridge.graphObjects.graphio import oxbridgePrYonGraphML as oxGraphML
from oxbridge.graphObjects.graphbase import oxbridgeNode  as Node
from oxbridge.graphObjects.graphbase import oxbridgeEdge  as Edge
from oxbridge.graphObjects.graphbase import oxbridgeGraph as Graph
from oxbridge.graphTools.autoGraphBuilders import autoGraphBase as autoGraph
from oxbridge.graphTools.autoGraphBuilders import parametricGraphBuilder as parametricGraph
import oxbridge.graphTools.parametric as pm

import math
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt




#=============================================================================
#            Examples using the oxbridge graph manipulation library
#=============================================================================

# ---------------
# Test the oxbridgePrYonGraphML read and write capabilities
def testgraphReadAndWrite(ingraph, outgraph):    
    # test input
    masterGraph = oxGraphML()
    masterGraph.setInputFile(ingraph)
    masterGraph.readAll()
    
    # test outpu
    masterGraph.setOutputFile(outgraph)
    masterGraph.writeAll()

# ---------------
# Test the scriptable creation of a simple 2-node graph and 
# write the graph to an PrYon formatted graphml file
def testTwoNodeGraph(outfile):
    # Make a test graph
    G = Graph() 
    
    # Make two nodes
    N1 = Node(1)
    N2 = Node(2)
    
    # Set the nodal (x,y,z) coordinates
    N1.setCoords(0.0,0.0,0.0)
    N2.setCoords(1.0,0.0,0.0)
    
    # Set the nodal data (name, region, freesurfer-name, hemisphere)
    N1.setData("lh.lobe_1", "cortical", "lobe_1", "left")
    N2.setData("rh.lobe_1", "cortical", "lobe_1", "right")
    
    # Add the nodes to the graph
    G.addNode(1, N1)
    G.addNode(2, N2)
    
    # now create an edge between them
    E1 = Edge(1,2,N1,N2)
    
    # set the edge properties
    E1.setEdgeWeights(1.0,2.0)
    E1.setOptionalValues(1.0,1.0,1.0)
    
    # add the edge to the graph
    G.addEdge(E1.getID(), E1)
    
    # Now create a oxbridgePryonGraphML object to write the graph
    # to disk in the PrYon format
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(G)
    graphout.writeAll()
    
# ---------------
def testThreeNodeAutoBuilder(outfile):

    # Make two nodes
    N1 = Node(1)
    N2 = Node(2)
    N3 = Node(3)
    
    # Set the nodal (x,y,z) coordinates
    N1.setCoords(0.0,0.0,0.0)
    N2.setCoords(1.0,0.0,0.0)
    N3.setCoords(0.5,0.5,0.0)
    
    # Set the nodal data (name, region, freesurfer-name, hemisphere)
    N1.setData("lh.alobe_1", "cortical", "alobe_1", "left")
    N2.setData("rh.blobe_1", "cortical", "blobe_1", "right")
    N3.setData("lh.clobe_1", "cortical", "clobe_1", "left")

    # now create edges between them.  The autoGraphBuilder knows
    # how to associate the actual nodes so all we need to specify
    # is the numeric IDs of the nodes for the edge
    E1 = Edge(1,2)
    E2 = Edge(1,3)
    E3 = Edge(2,3)
    
    # set the edge properties.  We can skip the "optional" properties
    # if we are using these files in PrYon.  
    E1.setEdgeWeights(1.0,2.0)
    E2.setEdgeWeights(2.0,3.0)
    E3.setEdgeWeights(3.0,4.0)

    myNodes = [N1,N2,N3]
    myEdges = [E1,E2,E3]
    
    # Assemble the graph object from the nodes and edges
    autoBuild = autoGraph(myNodes, myEdges)
    G = autoBuild.generate()
    
    # write the graph object to a file
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(G)
    graphout.writeAll()
    
    
# ---------------

# create a simple parametric curve (x(t),y(t),z(t)) = (t, t, t)
def testParametricLine(outfile):    
    tfunc = lambda t : t
    pcurvt = pm.parametricFunction(tfunc, lowerbd=0.0, upperbd=5.0)
    cline = pm.parametric3DCurve(pcurvt, pcurvt, pcurvt)
    
    pGraph = parametricGraph()
    pGraph.setDiscretization(1000)
    
    # add the parametric curve (with ID equal to 1)
    # and generate the resulting graph
    pGraph.addParametric3DCurve(cline, 1)
    Gline = pGraph.generate()
    
    # write the graph object to a file
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(Gline)
    graphout.writeAll()
    

# create a graph from two connected parametric curves 
# 1. (x(t),y(t),z(t)) = (t, t, t)
# 2. (X(t),Y(t),Z(t)) = (t, t, 0)
def testTwoParametricLines(outfile):    
    tfunc = lambda t : t
    zfunc = lambda t : 0.0
    
    # 
    pcurvt = pm.parametricFunction(tfunc)
    zcurvt = pm.parametricFunction(zfunc)
    
    cline = pm.parametric3DCurve(pcurvt, pcurvt, pcurvt)
    xyline = pm.parametric3DCurve(pcurvt, pcurvt, zcurvt)
    
    pGraph = parametricGraph()
    pGraph.setDiscretization(100)
    
    # add the parametric curve (with ID equal to 1)
    # and generate the resulting graph
    pGraph.addParametric3DCurve(cline, 1)
    pGraph.addParametric3DCurve(xyline, 2)
    
    
    Gline = pGraph.generate()
    
    # write the graph object to a file
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(Gline)
    graphout.writeAll()

# ----------
def testParametricCircle(outfile):    
    cosfunc = lambda t : math.cos(t)
    sinfunc = lambda t : math.sin(t)
    zrofunc = lambda t : 0.0
    
    cost = pm.parametricFunction(cosfunc, lowerbd=0.0, upperbd=2.0*math.pi)
    sint = pm.parametricFunction(sinfunc, lowerbd=0.0, upperbd=2.0*math.pi)
    zrot = pm.parametricFunction(zrofunc, lowerbd=0.0, upperbd=2.0*math.pi)
    
    circl = pm.parametric3DCurve(cost, sint, zrot)
     
    pGraph = parametricGraph(plist=[circl], N=100)
    Gline = pGraph.generate()
    
    # write the graph object to a file
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(Gline)
    graphout.writeAll()
    
# ----------
def testParametricCylinder(outfile):
    cosfunc = lambda t : math.cos(t)
    sinfunc = lambda t : math.sin(t)
    
    zplanes = []
    for z in range(0,100):
        # this might help if there is confusion around the lambda 
        # function created below
        # https://stackoverflow.com/questions/2295290/what-do-lambda-function-closures-capture
        zplanes.append(lambda t, z=z:0.02*z) 
    
    # x-y plane circle
    cost = pm.parametricFunction(cosfunc, lowerbd=0.0, upperbd=2.0*math.pi)
    sint = pm.parametricFunction(sinfunc, lowerbd=0.0, upperbd=2.0*math.pi)
    
    # elevated z functions
    zparametrics = [pm.parametricFunction(zp, lowerbd=0.0, upperbd=2.0*math.pi) for zp in zplanes]
        
    circs = [pm.parametric3DCurve(cost, sint, zpc) for zpc in zparametrics]
    
     
    pGraph = parametricGraph(plist=circs, N=100)
    GCyl = pGraph.generate()
    
    # write the graph object to a file
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(GCyl)
    graphout.writeAll()
# ---------------


# ---------------
def testUnitPlane(outfile):
    ident = pm.parametricFunction(lambda t : t)
    zerot = pm.parametricFunction(lambda t : 0.00)
    
    ylevels = []
    for y in range(0,100):
        ylevels.append(lambda t, y=y: y*(1/100))
    
    # y level set functions
    yparametrics = [pm.parametricFunction(yp) for yp in ylevels]
    
    # 3D parametric lines
    lines = [pm.parametric3DCurve(ident, ypc, zerot) for ypc in yparametrics]
    
    pGraph = parametricGraph(plist=lines, N=100)
    Glines = pGraph.generate()
    
    # write the graph object to a file
    graphout = oxGraphML()
    graphout.setOutputFile(outfile)
    graphout.setInnerGraph(Glines)
    graphout.writeAll()
    
    
# ---------------
def testNetworkX(infile):

    # compute a histogram showing the frequency of the (intermediate) constituents
    # in the shortest paths between the IDs of srcids and the IDs of 
    # trgids 
    # srcids: list of source ids
    # trgids: list of target ids
    # namedict: a dictionary whose keys are ids and values are the region names
    #           associated to that ID
    # Gx: Network x graph.  The IDs of srcid, trgid must be within the range 
    #   of the node IDs of this graph. (e.g. correspond to value node IDs)
    def pathHistogram(srcids, trgids, nodestolabels, Gx):
        
        # These will be the bin names
        binids = set({})
        fullresults = []

        for srcid in srcids:
            for trgid in trgids:
                path = nx.shortest_path(Gx, source=srcid, target=trgid)
                
                addtobins = None
                bDirect = False
                
                # the shortest path is a direct connection
                plen = len(path)
                if plen == 2:
                    addtobins = {'Direct Shortest Path'}
                    bDirect = True
                else:
                    # remove the first and last elemenet.  These are 
                    # the source and target and so we don't want them 
                    # to obscure the histogram results 
                    path.pop(plen-1)
                    path.pop(0)
                    addtobins = {nodestolabels[p] for p in path}
                
                binids.update(addtobins)
                
                if bDirect:
                    fullresults.append('Direct Shortest Path')
                else:
                    for p in path:
                        fullresults.append(nodestolabels[p])
        
        npbins = list(binids)
        npbins.sort()
        
        # we make a histogram by hand out of a dictionary which can 
        # be plotted using any auxilliary package.  
        myhist = {label: 0 for label in npbins}
        
        for f in fullresults:
            myhist[f] = myhist[f] + 1
        
        #retv = np.histogram(fullresults, bins=npbins)
        return myhist
    # -------------------------------------------------
    
    # get a list of the names of the neighbors of the node with id `nodeid'
    def getNeighborList(nodeids, nodestolabels, Gx):
        allneighbors = []
        
        for nodeid in nodeids:
            neighids = Gx.neighbors(nodeid)    
            neighbortags = [nodestolabels[nid] for nid in neighids]
            allneighbors.append(neighbortags)
        
        # remove duplicates that may be in the list
        allneighbors = list(dict.fromkeys(allneighbors))
        
        return allneighbors
    
    # returns the node IDs and weights of the top `nneigh' neighbors
    def getStronglyConnected(nodeid, G, nneigh=10):
        
        # get the degree of the node for safe bounds checking
        ndeg = G.degree[nodeid]
        
        if nneigh > ndeg:
            print(f"{nneigh} strongest connections requested by node {nodeid} has degree {ndeg}. using nneigh={ndeg}")
            nneigh = ndeg
        
        resv = sorted(G[nodeid].items(), key=lambda e: e[1]["weight"], reverse=True)[:nneigh]
        
        resd = {}
        
        # the `sorted' code above returns a list of tuples.  The first entry 
        # is a node ID of a neighbor and the second entry is a dictionary that
        # contains the edge information for the edge connecting `nodeid' to that 
        # neighbor.  We make an easily readable dictionary out of this by extracting
        # the neighbor node id and the value of the connecting edge 'weight' attribute
        for r in resv:
            resd[int(r[0])] = float(r[1]['weight'])
        
        return resd
    
    # == End internal functions ==
    # ============================

    masterGraph = oxGraphML()
    masterGraph.setInputFile(infile)
    masterGraph.readAll()
    
    # get the inner graph of the reader object
    G = masterGraph.getInnerGraph()
    
    # request a networkx version of the inner graph
    xG = G.getNetworkXGraph()    
    
    #------------------------------------------------
    
    associationCortexLeft = ['cortical.rostralmiddlefrontal.left',\
                                'cortical.caudalmiddlefrontal.left',\
                                'cortical.inferiorparietal.left',\
                                'cortical.superiortemporal.left']
    
    associationCortexRight = ['cortical.rostralmiddlefrontal.right',\
                                'cortical.caudalmiddlefrontal.right',\
                                'cortical.inferiorparietal.right',\
                                'cortical.superiortemporal.right']
    
    associationCortexAll = ['cortical.rostralmiddlefrontal.left',\
                                'cortical.rostralmiddlefrontal.right',\
                                'cortical.caudalmiddlefrontal.left',\
                                'cortical.caudalmiddlefrontal.right',\
                                'cortical.inferiorparietal.left',\
                                'cortical.inferiorparietal.right',\
                                'cortical.superiortemporal.left',\
                                'cortical.superiortemporal.right']

    #------------------------------------------------
    hippocampalLeft = 'subcortical.Left-Hippocampus.left'
    hippocampalRight = 'subcortical.Right-Hippocampus.right'
    hippocampalAll = ['subcortical.Right-Hippocampus.right', 'subcortical.Left-Hippocampus.left']
    #------------------------------------------------
    
    srcLeft = 'cortical.entorhinal.left'
    srcRight= 'cortical.entorhinal.right'
    
    #-----------------------------------------------
    
    # get the node names
    ndict = xG.nodes(data="name")
    namestonodes = {}
    nodestonames = {}
    
    # reverse the dictionary so that we can query by name
    # (multiple node IDs can be associated to each region name)
    for ntup in ndict:
        ndid = ntup[0]
        ndnm = ntup[1]

        nodestonames[ndid] = ndnm
        
        if ndnm in namestonodes:
            namestonodes[ndnm].append(ndid)
        else:
            namestonodes[ndnm] = [ndid]
    
    # Now we are going to make a histogram that will record the regional 
    # participants in the shortest paths between the source nodes and the 
    # various target nodes 

    # probe the connections of the hippocampus
    hlid = namestonodes[hippocampalLeft][0]
    hrid = namestonodes[hippocampalRight][0]
    
    # get the top 10 connections for each of these
    tophl = getStronglyConnected(hlid, xG, nneigh=15)
    tophr = getStronglyConnected(hrid, xG, nneigh=15)
    
    # now we extract the names of the nodes in the top most strongly connected
    # hippocampal neighbors
    hlregions = []
    hrregions = []
    
    for n in tophl:
        hlregions.append(nodestonames[n])
    
    for n in tophr:
        hrregions.append(nodestonames[n])
        

    # check to see if any of these regions are in the association cortex ROI
    shared = []
    for a in associationCortexAll:
        if a in hlregions or a in hrregions:
            shared.append(a)
    
    print(shared)


    # Lets do the left source nodes first
    #srcNodesLeft = namestonodes[srcLeft]
    #trgNodesLeft = namestonodes[hippocampalRight]
    #hist = pathHistogram(srcNodesLeft, trgNodesLeft, nodestonames, xG)
    
    # extract the leading order entries of the list
    def leadingorder(topconn):
        kys = list(topconn.keys())
        topval = topconn[kys[0]]
        
        leading = {}
        
        for k in kys:
            if (topval / topconn[k]) <= 10:
                leading[k] = topconn[k]
        
        return leading
    
    def printAssocHistCollection(histograms):
         # print ascii versions of the histogram
         for header in associationHistograms:
             print(f"= = = = = = = = = = {header} = = = = = = = = = =")
             print("")
             
             hist = histograms[header]
             
             # invert the histogram 
             inv_hist = {v: k for k, v in hist.items()}
             sortkeys = list(inv_hist.keys())
             sortkeys.sort()
             sortkeys.reverse()
             for val in sortkeys:
                 print(f"{inv_hist[val]}: {val}")
             print("================================================")
             print("")
    
    
    associationHistograms = {}
    for a in associationCortexLeft:
        associationHistograms[a] = {}
        
        nids = namestonodes[a]
        for n in nids:
            topconn = getStronglyConnected(n, xG, nneigh=15)
            
            leading = leadingorder(topconn)
            leadreg = []
            
            for r in leading:
                leadreg.append(nodestonames[r])
    
            # remove duplicates
            leadreg = list(dict.fromkeys(leadreg))
    
            for reg in leadreg:
                if reg in associationHistograms[a]:
                    associationHistograms[a][reg] += 1
                else:
                    associationHistograms[a][reg] = 1
    
    print("<><><><><><><><><> Histogram for Association Cortex Neighbor Connectivity <><><><><><><><><>")
    printAssocHistCollection(associationHistograms)
    # print ascii versions of the histogram
    #for header in associationHistograms:
    #    print(f"= = = = = = = = = = {header} = = = = = = = = = =")
    #    print("")
        
    #    hist = associationHistograms[header]
        
        # invert the histogram 
    #    inv_hist = {v: k for k, v in hist.items()}
    #    sortkeys = list(inv_hist.keys())
    #    sortkeys.sort()
    #    sortkeys.reverse()
        
    #    for val in sortkeys:
    #        print(f"{inv_hist[val]}: {val}")
        #for name in hist:
        #    print(f"{name}: " + '+'*hist[name])
        #    print("")
    #    print("================================================")
    #    print("")
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("")
    # Now we want to investigate the shortest paths between the 
    # entorhinal cortex and the nodes of the left association cortices
    leftAssociationShortestPathHistograms = {}
    srcnids = namestonodes[srcLeft]
    
    for a in associationCortexLeft:
        trgnids = namestonodes[a]
        leftAssociationShortestPathHistograms[a] = pathHistogram(srcnids, trgnids, nodestonames, xG)
    
  
    print("<><><><><><><><><> Histogram for EC-->Association Cortex Shortest Paths <><><><><><><><><>")
    printAssocHistCollection(leftAssociationShortestPathHistograms)
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    i=1
    
    # find the neighbors of the left and right hippocampus
    #hipleftid = namestonodes[hippocampalLeft]
    #hiprightid= namestonodes[hippocampalRight]
    
    #HPNLeft = getNeighborList(hipleftid, nodestonames, xG)
    #HPNRight = getNeighborList(hiprightid, nodestonames, xG)
    
    

def outputGraphLaplacianDiagonal(infile):
    # 0: length-free weights
    # 1: ballistic weights
    # 2: diffusive weights
    weightLengthPower = 2

    normalize = True

    masterGraph = oxGraphML()
    masterGraph.setInputFile(infile)
    masterGraph.readAll()
    
    G = masterGraph.getInnerGraph()

    # get the dictionary of neighborhoods for the graph.  This dictionary 
    # has keys equal to the node IDs and a value of the list which contains 
    # the node IDs of that node's neighbors
    ngbds = G.getAllNodeNeighborhoods()
    
    # sort the node IDs so that we can print an array that can 
    # be copied into C++
    nodeids = list(ngbds.keys())
    nodeids.sort()
    
    vals = []

    maxval = 0.0
    maxid = -1
    maxname = ""

    nodeclear = {}
    nodedvals = {}
    for n in nodeids:
        thisNgbd = ngbds[n]
        
        agg = 0.0
        for nbid in thisNgbd:
            # get the edge corresponding to n and its neighbor
            eid = (n, nbid)
            e = G.getEdge(eid)

            wij = e.getEdgeWeights()
            nij = wij[0]
            lij = wij[1]
            
            addwgt = nij/(lij**weightLengthPower)
            agg += addwgt
        
        vals.append(agg)
        
        if agg > maxval:
            maxval = agg
            maxid = n
            maxname = G.getNode(n).getName()


        nodenm = G.getNode(n).getGroupName()
        
        if agg in nodeclear:
            nodeclear[agg].append(nodenm)
        else:
            nodeclear[agg] = [nodenm]

        if nodenm in nodedvals:
            nodedvals[nodenm].append(agg)
        else:
            nodedvals[nodenm] = [agg]

    if normalize:
        maxv = np.amax(vals)
        vals = [v/maxv for v in vals]

    ranges = vals.copy()
    ranges.sort()
    
    
    # create a histogram of the values
    #plt.hist(vals, bins='auto')
    
    low = np.percentile(vals, 25)
    med = np.percentile(vals, 50)
    hgh = np.percentile(vals, 75)
    

    outstr=""
    print(*vals, sep="\n")
    outstr = f"initialVals = {str(vals)};"
    print(outstr)
   
    
   
    
# This function extracts a subgraph that contains a list of node IDs
# 
def extractSubgraph(G, nodelist=[], allEdges=False, quantile=0.75):
    pass

# This function takes as input a set of vertices and recursively builds a 
# superset of vertices.  The additional vertices in the superset  are determined 
# by taking an extension of the current set determined by adding those vertices 
# which are `strongly connected' to those in the base set.  The notion of 
# `strongly connected' is determined by taking those vertices whose edge weight 
# nij lies in the range of the input quantile (e.g. top 95% of nij values, etc).
#
# note: G must be an oxbridgeGraph object.  e.g.
#           masterGraph = oxGraphML()
#           masterGraph.setInputFile(infile)
#           masterGraph.readAll()
#           G = masterGraph.getInnerGraph()
#   
# nodelist must be a non-empty list of nodes.  This list represents a startpoint point
#   when newverts=[] and is used recursively thereafter to represent the new list.
def splayWellConnectedVertices(G, nodelist, newverts=[], vlimit=15, depthlimit=3, quantile=0.95):
    
    #-------------------------------------
    def splaynode(G, nid, quantile=0.75, matchHemi=True):
        allngbd = G.getNodeNeighbors(nid)
        ngbd = []
        
        if matchHemi == False:
            ngbd = allngbd
        else:
            nidhemi = (G.getNode(nid)).getHemisphere()
            for n in allngbd:
                if (G.getNode(n)).getHemisphere() == nidhemi:
                    ngbd.append(n)
            if len(ngbd) == 0:
                print(f"No nodes found in the same hemisphere of node {nid}")
                return None, None
        
        nijs = []
        lookup = {}
        foundneighbors = []
        
        # loop over is node's neighbors and get the associated edges
        for n in ngbd:
            foundedge = False
            edg = None
            
            if G.checkNodeExists(n) == False:
                print(f"node {n} not found")
                continue
                     
            if G.checkEdgeExists((nid,n)):
                edg = G.getEdge((nid,n))
                foundedge = True
            elif G.checkEdgeExists((n, nid)):
                edg = G.getEdge((n,nid))
                foundedge = True
        
            if foundedge:
                nij = (edg.getEdgeWeights())[0]
                nijs.append(nij)
                
                lookup[nij] = (n, edg)
            else:
                print(f"Could not find edge ({nid},{n}) or vice versa")
        
        
        npnij = np.array(nijs)
        qcut = np.quantile(npnij, quantile)
        
        # build the final neighbor ID list
        for k in list(lookup.keys()):
            if k > qcut:
                foundneighbors.append(lookup[k][0])
        
        return foundneighbors
    # -------------------------------------------

    if len(nodelist) == 0:
        print(f"Cannot splay the current vertex set, the list of nodes is empty.")
        return nodelist

    if (depthlimit == -1):
        print(f"Depth limit reached")
        return nodelist
    
    if len(nodelist + newverts) > vlimit:
        print(f"Node limit exceeded at depth {depthlimit}.  Returning the results of the last splay")
        return nodelist    
    
    # create an extended node set
    splayed = nodelist + newverts

    # if this is the first run, all of the vertices are `new'
    if len(newverts) == 0:
        newverts = nodelist
    
    extendby = []
    for v in newverts:
        newngb = splaynode(G, v, quantile=quantile, matchHemi=True)
        
        # gather novel vertices
        for s in newngb:
            if ((s in nodelist) == False) and ((s in newverts) == False) and ((s in extendby) == False):
                extendby.append(s)

    # call recursively
    return splayWellConnectedVertices(G, splayed, newverts=extendby, vlimit=vlimit, depthlimit=depthlimit-1, quantile=quantile)


def extractSubgraphFromVertices(G, vertexIDs):
    # the extracted graph
    eG = Graph()

    # first, we need to add all of the nodes to the graph
    for vid in vertexIDs:
        if G.checkNodeExists(vid) == False:
            print(f"Node with ID {vid} does not exist in the graph G")
            continue
        
        # add the node to the graph
        eG.addNode(vid, G.getNode(vid))
        
        
    # now we add the edges to the graph
    for vid in vertexIDs:
        for ng in vertexIDs:
            # skip the self node
            if ng == vid:
                continue
            
            if (G.checkEdgeExists( (vid, ng) ) ):
                eG.addEdge( (vid, ng), G.getEdge((vid, ng)) )
            else:
                continue
                # no edge (vid, ng) exists in the graph.  (ng, vid) might 
                # exist, but we will encounter this edge later.
    return eG

# This function extracts a well-connected subgraph that contains a given set of NodeIDs.  
#   infile: an oxmbm-formatted graphml file that contains the node IDs specified in the `nodeids' input list
#   primaryNodeids: a list of node ids (appearing in `infile') that act as the primary nodes for the subgraph extraction
#   quantile: keep the nodes connection strengths are in the indicated quantile
def extractWellConnectedSubgraph(infile, primaryNodeids=[], quantile=0.75):

    #-------------------------------------
    def getCutoffNeighborhood(G, nid, quantile=0.75, matchHemi=True):
        allngbd = G.getNodeNeighbors(nid)
        ngbd = []
        
        if matchHemi == False:
            ngbd = allngbd
        else:
            nidhemi = (G.getNode(nid)).getHemisphere()
            for n in allngbd:
                if (G.getNode(n)).getHemisphere() == nidhemi:
                    ngbd.append(n)
            if len(ngbd) == 0:
                print(f"No nodes found in the same hemisphere of node {nid}")
                return None, None
        
        nijs = []
        lookup = {}
        
        foundedges = []
        foundneighbors = []
        
        # loop over is node's neighbors and get the associated edges
        for n in ngbd:
            foundedge = False
            edg = None
            
            if G.checkNodeExists(n) == False:
                print(f"node {n} not found")
                continue
            
            
            if G.checkEdgeExists((pn,n)):
                edg = G.getEdge((pn,n))
                foundedge = True
            elif G.checkEdgeExists((n, pn)):
                edg = G.getEdge((n,pn))
                foundedge = True
        
            if foundedge:
                nij = (edg.getEdgeWeights())[0]
                nijs.append(nij)
                
                lookup[nij] = (n, edg)
            else:
                print(f"Could not find edge ({pn},{n}) or vice versa")
        
        
        npnij = np.array(nijs)
        qcut = np.quantile(npnij, quantile)
        
        # build the final neighbor ID list
        for k in list(lookup.keys()):
            if k > qcut:
                #final.append( (lookup[k])[0])
                foundneighbors.append(G.getNode( (lookup[k])[0] ))
                foundedges.append( (lookup[k])[1] )
        
        return foundneighbors, foundedges

    # -------------------------------------------
    # this routine completes the constructed subgraph at the indicated quantile
    # --> Fill me in: this function 
    def completeSubgraph(G, MG, quantile):
        return MG
    
    # -------------------------------------------


    if len(primaryNodeids) == 0:
        print("The list of primary node IDs must not be empty")
        return

    masterGraph = oxGraphML()
    masterGraph.setInputFile(infile)
    masterGraph.readAll()

    G = masterGraph.getInnerGraph()
    ngbds = G.getAllNodeNeighborhoods()


    allnodes = []
    alledges = []
    
    # create a new graph to start appending objects to
    MG = Graph()
    
    for pn in primaryNodeids:
        
        if G.checkNodeExists(pn) == True:
            
            # add the primary node
            if (MG.checkNodeExists(pn)) == False:
                MG.addNode(pn, G.getNode(pn))
            
            nodes, edges = getCutoffNeighborhood(G, pn, quantile=quantile)
            
            for n in nodes:
                nid = n.getID()
                if (MG.checkNodeExists(nid)) == False:
                    MG.addNode(nid, n)
            
            for e in edges:
                eid = e.getID()
                if (MG.checkEdgeExists(eid)) == False:
                    MG.addEdge(eid, e)

        else:
            print(f"Node with id {pn} does not exist in the {infile}")
    
    
    MGF = completeSubgraph(G, MG, quantile)
    
    return MGF


# entry point of the script 
if __name__ == "__main__": 
    # ---------------------------------------
    # Test graph input / output
    # ---------------------------------------
    #testgraphReadAndWrite("/home/zcxp/src/oxford/oxbdev/master-std500.graphml","/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")

    # ---------------------------------------
    # Test automated graph creation
    # ---------------------------------------
    #testTwoNodeGraph("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")
    #testThreeNodeAutoBuilder("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")    
    #testParametricLine("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")
    #testTwoParametricLines("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")
    #testParametricCircle("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")
    #testParametricCylinder("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")
    #testUnitPlane("/home/zcxp/src/oxford/oxbdev/oxbridge.graphml")
    #outputGraphLaplacianDiagonal("/home/zcxp/devops/oxb/master-std500.graphml")
    
    
    
    # --------------------------
    # Test the extraction of a well-connected subgraph from a master graph
    # --------------------------
    #extractWellConnectedSubgraph("/home/zcxp/devops/oxb/master-std33.graphml", primaryNodeids=[68, 71, 66, 70], quantile=0.5)
    #extractWellConnectedSubgraph("/home/zcxp/devops/oxb/master-std33.graphml", primaryNodeids=[68, 73, 81, 71, 66, 70, 60, 64, 72], quantile=0.75)
    #extractWellConnectedSubgraph("/home/zcxp/devops/oxb/master-std33.graphml", primaryNodeids=[68, 73, 81, 60, 72, 76, 71, 66, 70, 64, 75, 63, 65], quantile=0.75)
    masterGraph = oxGraphML()
    masterGraph.setInputFile("/home/zcxp/devops/oxb/master-std33.graphml")
    masterGraph.readAll()

    G = masterGraph.getInnerGraph()
    
    # initial vertex count (produced a quotient space that was too large)
    #splayedVertices = splayWellConnectedVertices(G, nodelist=[68, 71, 66, 70], vlimit=19, depthlimit=10, quantile=0.97)
    #splayedGraph = extractSubgraphFromVertices(G, splayedVertices)
    
    # smaller vertex counts (the resulting homotopy space size is here: https://oeis.org/A331554)
    vlimit=17
    splayedVertices = splayWellConnectedVertices(G, nodelist=[68, 71], vlimit=vlimit, depthlimit=10, quantile=0.97)
    splayedGraph = extractSubgraphFromVertices(G, splayedVertices)
    
    print(f"Found final vertex set {splayedVertices} with length {len(splayedVertices)}")
    
    eids = splayedGraph.getEdgeIDs()
    

    print(f"Edges needed are at least {(vlimit-1)*(vlimit-2)/2 + 1}")
    print(f"Complete graph has {vlimit*(vlimit-1)/2} vertices")
    print(f"Splayed graph has {len(eids)} edges")
    
    subGraph = oxGraphML()
    subGraph.setInnerGraph(splayedGraph)
    
    subGraph.setOutputFile("/home/zcxp/devops/oxb/adsubtype-graph-std33.graphml")
    subGraph.writeAll()
    
    # write the graph laplacian.  set lengthpower=X for different types
    #   X = 0: unweighted
    #   X = 1: ballistic weights
    #   X = 2: diffusive weights
    subGraph.writeGraphLaplacianCSV(csvpath="/home/zcxp/devops/oxb/", lengthpower=2)
    
    # ---------------------------------------
    # Test Network-X related functionality
    # ---------------------------------------
    #testNetworkX("/home/zcxp/src/oxford/oxbdev/master-std500.graphml")
    
    
    
    
    
