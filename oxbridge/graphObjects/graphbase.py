#=============================================================================
# ** Oxford Mathematical Brain Modeling Group **
#   This source file defines the Node and Edge objects 
#   used by the graphml readers and writers. 
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

import networkx as nx

# Class representation of a node
class oxbridgeNode:
    
    def __init__(self, nid):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.id = int(nid)
        
        # --- Primary data
        self.region = ""
        self.fsname = ""
        self.hemisphere = ""
        self.name = ""
        
        # --- Secondary data
        self.fsGroup = ""
        self.fsGroupNumber = ""
        self.sepdelim = ""
        
        self.hemiGroup = ""
        self.hemiGroupNum = ""
        
        self.pryonName = ""

    def __setDerivedData(self, delim='_'):
        self.sepdelim=delim
        
        yat = self.name.find(delim)
        xat = self.fsname.find(delim)
        
        # ----- This grouping ID differentiates the hemispheres
        # i.e. rh.[name] from lh.[name]
        if yat != -1:
            self.hemiGroup = self.name[0:yat]
            self.hemiGroupNum = self.name[yat+1:]
        else:
            self.hemiGroup = self.name
            self.hemiGroupNum = ""
        
        # ---- This grouping ID does not differentiate the hemispheres
        # and is therefore bilateral.  E.g. both rh.[groupname] and lh.[groupname]
        # are produce the same labels.
        if xat != -1:
            self.fsGroup = self.fsname[0:xat]
            self.fsGroupNumber = self.fsname[xat+1:]
        else:
            self.fsGroup = self.fsname
            self.fsGroupNumber = ""
            
         
        self.pryonName = f"{self.region}.{self.fsGroup}.{self.hemisphere}"
        
    def setCoords(self, xin, yin, zin):
        self.x = xin
        self.y = yin
        self.z = zin
    
    
    def setData(self, name, region, fsname, hemisphere):
        self.name = name
        self.region = region
        self.fsname = fsname
        self.hemisphere = hemisphere
        
        self.__setDerivedData()
    
    
    def addNeighbor(self, edge, node):
        self.neighborhood[edge] = node
 
    def getID(self):
        return self.id
    
    def getRegion(self):
        return self.region
    
    def getFSName(self):
        return self.fsname
    
    def getHemisphere(self):
        return self.hemisphere
    
    def getName(self):
        return self.name
    
    def getGroupName(self):
        return self.hemiGroup
    
    def getNodeNumberInGroup(self):
        retv = 1
        
        if self.hemiGroupNum != "":
            retv = int(self.hemiGroupNum)
        
        return retv
    
    def getBilaterialGroupName(self):
        return self.fsGroup
    
    def getCoords(self):
        return (self.x, self.y, self.z)
    
    def getPryonName(self):
        return self.pryonName
    
# Class representation of an edge
class oxbridgeEdge:
    
    def __init__(self, src, trg, srcnode = None, trgnode = None):
        self.sid = src
        self.tid = trg
        self.nij = 1.0
        self.lij = 1.0
        
        # optional edge attributes
        self.lfiberstd = -1.0
        self.famean= -1.0
        self.fastd = -1.0
        
        self.srcnode = srcnode
        self.trgnode = trgnode
    
    def getID(self):
        return (self.sid, self.tid)
    
    def associateSource(self, sourcenode):
        self.srcnode = sourcenode
    
    def associateTarget(self, targetnode):
        self.trgnode = targetnode
        
    # set the edge connectivity number and length
    def setEdgeWeights(self, n, l):
        self.nij = float(n)
        self.lij = float(l)
    
    # set optional values for the edge.  These are
    # fiberstd: Fiber length standard deviation
    # famean: Fiber mean fractional anisotropy
    # fastd : Fiber fractional anisotropy standard deviation
    def setOptionalValues(self, fiberstd, famean, fastd):
        self.lfiberstd = float(fiberstd)
        self.famean = float(famean)
        self.fastd = float(fastd)
        
    # return the edge weights (note that lij is the 
    # mean fiber length and so is also returned by
    # getFiberLenthVals)
    def getEdgeWeights(self):
        return (self.nij, self.lij)
    
    # return the mean and standard deviation of 
    # the fractional isotropy for the edge
    def getFAVals(self):
        return (self.famean, self.fastd)
    
    # return the mean and standar deviation for 
    # fiber lengths and values
    def getFiberLengthVals(self):
        return (self.lij, self.lfiberstd)
    
class oxbridgeGraph:
    
    
    def __init__(self):
        # Dictionary: Node id --> node object
        self.nodes = {}
        
        # Dictionary: Edge src/target tuple (e.g. (1,2)) --> edge object
        self.edges = {}
        
        # Dictionary: Node id --> All associated edge object tuples
        self.nodeEdges = {}
        
        # Dictionary: Node id --> Anatomical region name
        self.nodeRegions = {}
        
        # Dictionary: Node id --> All associated neighbor node IDs
        self.nodeNeighbors = {}
        
        
        
    # add a node to the graph
    #   nid: integer node id 
    #   n:   a node object
    # [optional] bilateral = true identifies all nodes from the same region, 
    # on either side, with each other. E.g. it puts lh.entorhinal in the same 
    # group as rh.entorhinal
    def addNode(self, nid, n, bilateral=False):
        self.nodes[nid] = n
        
        rgn = n.getGroupName()
        
        if bilateral:
            rgn = n.getBilaterialGroupName()
        
        # ---- 
        # Add the node to its region list
        if rgn in self.nodeRegions:
            self.nodeRegions[rgn].append(nid)
        else:
            self.nodeRegions[rgn] = [nid]
        # ----
        

    def getNodeIDs(self):
        return list(self.nodes.keys())
    
    def getEdgeIDs(self):
        return list(self.edges.keys())


    def checkEdgeExists(self, eid):
        res = (eid in self.edges)
        return res
    
    def checkNodeExists(self, nid):
        res = (nid in self.nodes)
        return res
    
    
    def getNode(self, nid):
        nV = None
        
        if self.checkNodeExists(nid):
            nV = self.nodes[nid]
        
        return nV
    
    
    def getEdge(self, eid):
        eV = None
        
        if self.checkEdgeExists(eid):
            eV = self.edges[eid]
        else:
            # Its possible that the edge exists but in 
            # a reversed order.  Since our graph is 
            # undirected, lets check for that 
            erev = tuple(reversed(eid))
            if self.checkEdgeExists(erev):
                eV = self.edges[erev]
            
        return eV

    def getNodeNeighbors(self, nid):
        nB = None
        
        if self.checkNodeExists(nid):
            nB = self.nodeNeighbors[nid]
        
        return nB

    def getAllNodeNeighborhoods(self):
        return self.nodeNeighbors
    
        
    # add an edge 
    #   srctrg: a tuple with two integers indicating source and 
    #           target node IDs. Example (1,2)
    #   e: an edge object to add
    def addEdge(self, srctrg, e):
        
        if type(srctrg) == tuple:
            srcid = srctrg[0]
            trgid = srctrg[1]
            
            if srcid in self.nodes and trgid in self.nodes:
                self.edges[srctrg] = e
                
                #-----
                # add to the neighbor list
                if srcid in self.nodeNeighbors:
                    self.nodeNeighbors[srcid].append(trgid)
                else:
                    self.nodeNeighbors[srcid] = [trgid]
                    
                if trgid in self.nodeNeighbors:
                    self.nodeNeighbors[trgid].append(srcid)
                else:
                    self.nodeNeighbors[trgid] = [srcid]
                #----
                
                #----
                # append the edge (tuple) identifier to the node edge list
                for iid in srctrg:
                    if iid in self.nodeEdges:
                        self.nodeEdges[iid].append(srctrg)
                    else:
                        self.nodeEdges[iid] = [srctrg]
                #----
            else:
                print(f"Cannot find two nodes with source and target id {srcid}, {trgid}")
            
        else:
            print("** Edge not added")
            print("srctag should be a tuple of integer node IDs, e.g. (1,2)")
        
    
    # Return a networkx graph object of this graph
    # [Optional]
    #   Specify the weights as
    #   'N', 'Ballistic' or 'Diffusive'
    # These options determine the edge weights as nij/(lij)^K where
    # 'N': K=0
    # 'Ballistic': K=1
    # 'Diffusive': K=2
    def getNetworkXGraph(self, weights='Diffusive'):
        xG = None
        
        weightOptions = {'N': 0, 'Ballistic': 1, 'Diffusive': 2}
        
        if len(self.nodes) == 0 or len(self.edges) == 0:
            print("** Cannot build a NetworkX graph.  This object is missing nodes or edges.")
            print("** Please call `addNode' or `addEdge' to add the missing requirements")
            return xG
                
        if weights not in weightOptions:
            print(f"** Invalid option weights={weights}.  Please use one of: \'N\', \'Ballistic\', \'Diffusive\'")
            return xG
    
        K = weightOptions[weights]
    
        # start graph construction
        xG = nx.Graph()
        
        for nnum in self.nodes:
            n = self.nodes[nnum]
            nid = n.getID()
            ncoor = n.getCoords()
            ntag  = n.getPryonName()
            nodeInGroup = n.getNodeNumberInGroup()
            
            # add the node to the graph
            xG.add_node(nid, name=ntag, numInGroup=nodeInGroup, xcoor=ncoor[0], ycoor=ncoor[1], zcoor=ncoor[2])
            
            
        for en in self.edges:
            e = self.edges[en]
            
            # this returns the (source id, target id) pair
            eid = e.getID()
            ewgts = e.getEdgeWeights()
            nij = ewgts[0]
            lij = ewgts[1]
            
            # compute the edge weight
            ew = nij / (lij**K)
    
            xG.add_edge(eid[0], eid[1], weight=ew, n=nij, l=lij, typ=weights)
            
        return xG
    
    
    
    
    
    
    
    
    
    
    