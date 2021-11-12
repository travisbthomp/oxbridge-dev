#=============================================================================
# ** Oxford Mathematical Brain Modeling Group **
#   This source file defines graphBuilder tool.  This tool generates a graph 
#   from a collection of oxbridgeNodes and oxbridgeEdges
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

from oxbridge.graphObjects.graphbase import oxbridgeNode, oxbridgeEdge, oxbridgeGraph
#import oxbridge.graphTools.parametric as pm

import math

# This class builds a graph from a list of nodes and edges.  The nodes and 
# edges do not need to be associated beforehand.  Auto-associating the nodes 
# with the edges is the primary purpose of the graphBuilder object
class autoGraphBase:
    
    def __init__(self, nodelist=None, edgelist=None):
        self.allnodes = nodelist
        self.alledges = edgelist
        self.__nodedict = None
        self.generatedG = None
    
    def __buildNodeDictionary(self):
        self.__nodedict = {}
        
        for n in self.allnodes:
            nid = n.getID()
            self.__nodedict[nid] = n
        
    def setNodeList(self, nodelist):
        self.allnodes = nodelist
    
    def setEdgeList(self, edgelist):
        self.alledges = edgelist
    
    # assembles a graph from a list of nodes and edges
    def assembleGraph(self):
        
        G = None
        
        if self.allnodes == None:
            print("** Please associate a node list to the graph")
            print("** by calling setNodeList(myNodeList)")
            return G
        
        if self.alledges == None:
            print("** Please associate an edge list to the graph")
            print("** by calling setEdgeList(myEdgeList)")
            return G
        
        self.__buildNodeDictionary()
        
        G = oxbridgeGraph()
        
        # first, add all the nodes
        for n in self.allnodes:
            nid = n.getID()
            G.addNode(nid, n)
            
        for e in self.alledges:
            # get the edge ID
            eid = e.getID()
            
            srcnid = eid[0]
            trgnid = eid[1]
            
            # associate the nodes in case this
            # was not already done
            e.associateSource(self.__nodedict[srcnid])
            e.associateTarget(self.__nodedict[trgnid])
            
            # add the edge
            G.addEdge(eid, e)

            self.generatedG = G

        return G

    # Override this function to provide class-specific implementations
    def generate(self):
        return self.assembleGraph()

    # return the inner graph
    def getGraph(self):
        return self.generatedG



# =============================================================================
# =============================================================================
        
    
    
# This class represents an autograph which builds itself from 
# a family of parametric curves.  The order in which you specify 
# the parametric curves is important.  Given a set of M parametric 
# curves, {p1, p2, .. pM} this autograph creates the following geometry:
#
# 1. Creates N nodes on each curve
# 2. Creates edges between node J and node J+1 of curve pK
# 3. Creates edges between node J of curve pK and node J of curve p(K+1)
# 
# Note: Each curve is considered as its own "graph region" for PrYon labelling
#   purposes.  You can further specify that each point on each curve is 
#   also a separate graph region by using the function setLabelPointsAsRegions 
#   with (byPoint=True) to toggle this behavior on or off
class parametricGraphBuilder(autoGraphBase):
    
    # [Optional] 
    #   plist - a list of (non intersecting) parametric3DCuve objects
    #   N - the number of nodes per curve
    #   regionPrefix - a region name that will be prepended to the file
    #       region IDs of every node. (e.g. `lh.entorhinal' etc)
    #   hemisphere - set to "right" or "left"
    def __init__(self, plist=None, N=10, regionPrefix="", hemisphere="left"):
        self.pmlist = plist
        self.pmdict = {}
        self.regionPrefix = regionPrefix
        self.pointLabels = False
                
        self.N=N
        
        if hemisphere != "right" and hemisphere != "left":
            print("** the hemisphere option must be either `left' or `right'")
            hemisphere = "left"
        
        self.hemisphere = hemisphere
        self.hemi = ""
        
        if self.hemisphere == "left":
            self.hemi = "lh"
        else:
            self.hemi = "rh"
        
        
        if plist != None:
            if isinstance(plist, list) == False:
                print("** optional argument `plist' must specify a list of parametric3DCurve objects")
                print("** e.g. [p1, p2, p3, .. pM].  Resorting to default behavior, please add the")
                print("** parametric3DCurve objects manually or fix this oversight")
            else:
                idx = 1
                for p in plist:
                    self.pmdict[idx] = p
                    idx += 1
                    
        # we will set the internal lists ourselves so 
        # we call the base constructor without arguments
        super().__init__()        


    def addParametric3DCurve(self, p3D, nid):
        if type(nid) != int:
            print(f"** The identified {nid} for the parametric3DCurve is not of type int.")
        else:
            self.pmdict[nid] = p3D


    # remove the parametric curve with index `nid' from the internal 
    # parametric curve list
    def removeParametric3DCurve(self, nid):
        success = (nid in self.pmdict)
        
        if success:
            del self.pmdict[nid]
        
        return success

    # set the number of nodes to generate per curve
    def setDiscretization(self, N):
        self.N = N

    # If byPoint = True: 
    #   each region label has `pR' appended where R is 
    #   the number of the node.  This implies that all points on all curves 
    #   will denote `separate regions'
    #
    # If byPoint = False: (Default behavior of constructor)
    #   each curve will it a region and the point numbers will be in 
    #   that region.  e.g. "lh.[regionlabel]cK_pM"
    def setLabelPointsAsRegions(self, byPoint):
        self.pointLabels = byPoint

    # class specific version of the generate function
    # [optional]
    #   Set closeConstruction=True to connect the points from the last curve 
    #   to those of the first curve. (Default: False)
    def generate(self, closeConstruction=False):
        # we create the nodes for each parametric function
        kys = list(self.pmdict)
        nkys = len(kys)
        kys.sort()

        # the curve IDs don't have to be indexed by integers but we do 
        # want the node numbers to be integer indexed
        masterlist = {}
        
        # generate all points on all curves.  We need to ensure that they are 

        ttlnds = 0
        NN = 0
        for k in kys:
            crv = self.pmdict[k]
            nds, clsd = crv.evalNPoints(self.N)
            masterlist[k]=nds
            NN = len(nds)
            ttlnds += NN
        
        closedCurves = (NN == self.N-1)
        
        # this is a necessary, but not sufficient, condition to check
        # if all curves are either closed or all open
        if ttlnds != (self.N)*nkys and ttlnds != (self.N-1)*nkys:
            print("** All curves must be either closed curves (N-1 points) or")
            print("** open curves (N points).  This class does not currently")
            print("** support a mix of open and closed parametric curves")
            return
        
        # For simplicity, we assume that NN is either self.N or self.N-1 and 
        # use it for further indexing.  This is not strictly enforced by the 
        # necessary condition, above
        
        print("----------------------------------------------------------------")
        print(f"--- Generating nodes and edges from {nkys} parametric curves --- ")
        print("----------------------------------------------------------------")
        
        # node and edge lists to be populated
        super().setNodeList([])
        super().setEdgeList([])
        
        # the keys of pmdict are integers by construction - so we use
        # them for indexing the curves and the points
        globalNodeID = 1
        for k in kys:
    
            stradd = f"{self.regionPrefix}c{k}"            
            ndprefix = self.hemi + f".{stradd}"
            fsprefix = stradd
                        
            crvnds = masterlist[k] 
            
            # create nodes 
            for nd in crvnds:
                
                # get the number of this node relative to the given curve
                ndoncv = ((globalNodeID-1) % NN) + 1
                
                # finish out the labels depending on the naming options
                nodestamp = f"_{ndoncv}"
                if self.pointLabels:
                    nodestamp = f"p{ndoncv}"
                
                ndname = ndprefix + nodestamp
                fsname = fsprefix + nodestamp
                
                # extract the coordinates of the parametric curve nodes
                ndx = nd[0]
                ndy = nd[1]
                ndz = nd[2]
                
                # create the node object, add data
                N  = oxbridgeNode(globalNodeID)
                N.setCoords(ndx, ndy, ndz)
                N.setData(ndname, "cortical", fsname, self.hemisphere)
                
                self.allnodes.append(N)
                globalNodeID += 1
                
        
        # Now we create the edges corresponding to the nodes.  This 
        # is done quite simply.  We do the following:
        # 1. Connect one node to its (on-curve) neighbor
        # 2. Connect each node to its next-curve neighbor
        ndxs = range(len(self.allnodes))
            
        stepn = len(kys)
        steps = range(stepn)
            
        lines = {}
        for s in steps:
            # these are (zero-based) indices into the nodes list
            # which partition the full nodes list back into lines
            lines[s] = ndxs[s*NN:(s+1)*NN]
            
        # Now we create the edges from the lines
        for s in steps:
            # create the edges for this line
            lineNdxs = lines[s]
                
            for idx in range(NN):
                    
                # Schematic of nodes:
                # 
                #   un
                #   |
                #  tn -- nn
                    
                # this node
                tn = self.allnodes[lineNdxs[idx]]
                    
                # nodes next neighbor (on this curve)
                nn = None
                    
                # nodes 'upper' neighbor
                un = None
                    
                    
                if idx != NN-1:
                    nn = self.allnodes[lineNdxs[idx+1]]
                else: # close the curve if we need to do so
                    if closedCurves == True:
                        nn = self.allnodes[lineNdxs[0]]
                    
                if s !=stepn-1:
                    un = self.allnodes[lines[s+1][idx]]
                else: # close the entire construct by linking the nodes 
                    # in the last curve to those in the first curve
                    if closeConstruction == True and stepn>1:
                        un = self.allnodes[lines[0][idx]]
                    
                tnc = tn.getCoords()
                src = tn.getID()
                    
                # add the details for the node-->node 
                # connection along this curve
                if nn != None:
                    nnc = nn.getCoords()
                    trg = nn.getID()
                    
                    # distance between the nodes
                    L = math.sqrt( (tnc[0]-nnc[0])**2 + (tnc[1]-nnc[1])**2 + (tnc[2]-nnc[2])**2 )
                    # set up and add the edge
                    E = oxbridgeEdge(src, trg, tn, nn)
                    E.setEdgeWeights(1.0, float(L))
                    E.setOptionalValues(1.0, 1.0, 1.0)
                    self.alledges.append(E)
                    
                    
                if un != None:
                    unc = un.getCoords()
                    trg = un.getID()
        
                    # distance between the nodes
                    L = math.sqrt( (tnc[0]-unc[0])**2 + (tnc[1]-unc[1])**2 + (tnc[2]-unc[2])**2 )
                    # set up the edge and add it
                    E = oxbridgeEdge(src, trg, tn, nn)
                    E.setEdgeWeights(1.0, float(L))
                    E.setOptionalValues(1.0, 1.0, 1.0)
                    self.alledges.append(E)
                        
        # Call the generation procedure in the base class
        return super().generate()
        
        