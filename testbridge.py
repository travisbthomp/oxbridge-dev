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
#import numpy as np
import networkx as nx



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
    
    
    for a in associationCortexLeft:
        nids = namestonodes[a]
        for n in nids:
            topconn = getStronglyConnected(n, xG, nneigh=15)
            topreg = []
            
            for r in topconn:
                topreg.append(nodestonames[r])
    
            print(topreg)
    
    
    
    # find the neighbors of the left and right hippocampus
    hipleftid = namestonodes[hippocampalLeft]
    hiprightid= namestonodes[hippocampalRight]
    
    HPNLeft = getNeighborList(hipleftid, nodestonames, xG)
    HPNRight = getNeighborList(hiprightid, nodestonames, xG)
    
    


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
    
    
    # ---------------------------------------
    # Test Network-X related functionality
    # ---------------------------------------
    testNetworkX("/home/zcxp/src/oxford/oxbdev/master-std500.graphml")
    
    
    
    
    
