#=============================================================================
# ** Oxford Mathematical Brain Modeling Group **
#   This source file defines the object that reads the 
#   OxMBM brain graphml files. 
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

import xml.etree.ElementTree as ET 
from oxbridge.graphObjects.graphbase import oxbridgeNode, oxbridgeEdge, oxbridgeGraph

class oxbridgeGraphWriterBase:
    
    def __init__(self):
        # add anything relevant here
        self.outfile = None
        
    def setOutputFile(self, flout):
        self.outfile = flout
        
        
    def writeHeader(self):
        print("** oxBridgeGraphWriterBase writeHeader() called.")
        print("** Your implementation should override this function.")
        
    def writeNodes(self):
        print("** oxBridgeGraphWriterBase writeGraphNodes() called.")
        print("** Your implementation should override this function.")
    
    def writeEdges(self):
        print("** oxBridgeGraphWriterBase writeGraphEdges() called.")
        print("** Your implementation should override this function.")
            
    def writeFooter(self):
        print("** oxBridgeGraphWriterBase writeFooter() called.")
        print("** Your implementation should override this function.")

    def writeAll(self):
        print("** oxBridgeGraphWriterBase writeFooter() called.")
        print("** Your implementation should override this function.")


class oxbridgeGraphReaderBase:
    
    def __init__(self):
        # add anything relevant here
        self.infile = None
    
    
    def readHeader(self):
        print("** oxbridgeGraphReaderBase readHeader() called.")
        print("** Your implementation should override this function.")

    def setInputFile(self, fln):
        self.infile = fln

    def readNodes(self):
        print("** oxbridgeGraphReaderBase readNodes() called.")
        print("** Your implementation should override this function.")

    def readEdges(self):
        print("** oxbridgeGraphReaderBase readEdges() called.")
        print("** Your implementation should override this function.")

    def readAll(self):
        print("** oxbridgeGraphReaderBase readAll() called.")
        print("** Your implementation should override this function.")


# This class reads a PrYon standardized-format graphml file into memory 
# which can be manipulated or output 
class oxbridgePrYonGraphML(oxbridgeGraphReaderBase, oxbridgeGraphWriterBase):

    def __init__(self):
        
        self.nodeAssociations = {"x-coord": "d0",  \
                                "y-coord" : "d1",  \
                                "z-coord" : "d2",  \
                                "node-id" : "d3",  \
                                "region"  : "d4",  \
                                "fsname"  : "d5",  \
                                "name"    : "d6",  \
                                "hemisphere": "d7"}
        
        self.edgeAssociations = {"nfibers": "d9",     \
                                 "famean" : "d10",    \
                                 "lfiberstd" : "d11", \
                                 "lfibermean": "d12", \
                                 "fastd": "d13"}
        
        
        self.nodesprocessed = False
        self.nodeswritten = False
    
        # the graph object we will read in 
        self.G = oxbridgeGraph()

        # the graphml (main) Element Tree 
        # and the graph (sub to graphml) Element tree used 
        # for the XML representation of the graph
        self.__graphml = None
        self.__graph   = None

    #------------------------------------------------------------------------          
    #------------------------------------------------------------------------          
    # Graph reader interface
    #------------------------------------------------------------------------          
    #------------------------------------------------------------------------          


    def __getNSKeyedStr(self, st):
        ns = 'http://graphml.graphdrawing.org/xmlns'
        ky = '{' + ns + '}' + st
        return ky


    def __processnode(self, node):
        
        # get the ID so we can create a new node
        iid = int(node.get('id'))
        n = oxbridgeNode(iid)
        
        ncoords = [0.0, 0.0, 0.0]
        nstrngs = ["", "", "", ""]
        
        # Make the keyed string for the 
        #   data elements
        datakys = self.__getNSKeyedStr('data')

        for dat in node.findall(datakys):
            #get the code in the file 
            # this will be something like d2, d3, etc
            code = dat.get('key')
            cval = dat.text
        
            if code == self.nodeAssociations["x-coord"]:
                ncoords[0] = float(cval)
            if code == self.nodeAssociations["y-coord"]:
                ncoords[1] = float(cval)
            if code == self.nodeAssociations["z-coord"]:
                ncoords[2] = float(cval)
            if code == self.nodeAssociations["name"]:
                nstrngs[0] = str(cval)
            if code == self.nodeAssociations["region"]:
                nstrngs[1] = str(cval)
            if code == self.nodeAssociations["fsname"]:
                nstrngs[2] = str(cval)
            if code == self.nodeAssociations["hemisphere"]:
                nstrngs[3] = str(cval)
      
        # set the data             
        n.setCoords(ncoords[0], ncoords[1], ncoords[2])
        n.setData(nstrngs[0], nstrngs[1], nstrngs[2], nstrngs[3])
        
        # add the node to the graph
        self.G.addNode(iid, n)
        
        
    def __processedge(self, edge):
        # get the source and target IDs
        srcid = int(edge.get('source'))
        trgid = int(edge.get('target'))

        srcn = self.G.getNode(srcid)
        trgn = self.G.getNode(trgid)

        if srcn == None or trgn == None:
            print(f"** Error: The nodes for edge ({srcid},{trgid}) do not exist in the graph.")
            return
       
        e = oxbridgeEdge(srcid, trgid, srcn, trgn)
       
        # Make the keyed string for the 
        #   data elements
        datakys= self.__getNSKeyedStr('data')
       
        edgeweights = [0.0, 0.0]
        edgeoptional = [-1.0, -1.0, -1.0]
        
        for dat in edge.findall(datakys):
           #get the code in the file 
           # this will be something like d2, d3, etc
           code = dat.get('key')
           cval = dat.text
           
           if code == self.edgeAssociations["nfibers"]:
               edgeweights[0] = float(cval)
           if code == self.edgeAssociations["lfibermean"]:
               edgeweights[1] = float(cval)
            
            
           if code == self.edgeAssociations["lfiberstd"]:
               edgeoptional[0] = float(cval)
           if code == self.edgeAssociations["famean"]:
               edgeoptional[1] = float(cval)
           if code == self.edgeAssociations["fastd"]:
               edgeoptional[2] = float(cval)

        # set weights, optional values
        e.setEdgeWeights(edgeweights[0], edgeweights[1])
        e.setOptionalValues(edgeoptional[0], edgeoptional[1], edgeoptional[2])
        
        # add the completed edge to the graph
        self.G.addEdge((srcid,trgid), e)

    
    def readNodes(self):
        tree = ET.parse(self.infile)
        root = tree.getroot()

        # get the graph main object
        gkey = self.__getNSKeyedStr('graph')
        graph = root.find(gkey)

        # the graph consists of node and
        #   edge children.  we want the 
        #   nodes here
        nkey = self.__getNSKeyedStr('node')
        allnodes = graph.findall(nkey)
        
        for n in allnodes:
            self.__processnode(n)
        
        self.nodesprocessed = True
            
    def readEdges(self):
        if self.nodesprocessed == False:
            print("** Please call readNodes() before readEdges() so that the appropriate links can be made at read time")
            return
        
        tree = ET.parse(self.infile)
        root = tree.getroot()

        # get the graph main object
        gkey = self.__getNSKeyedStr('graph')
        graph = root.find(gkey)
        
        nkey = self.__getNSKeyedStr('edge')
        alledges = graph.findall(nkey)

        for e in alledges:
            self.__processedge(e)


    # calls both read functions in the correct order
    def readAll(self):
        self.readNodes()
        self.readEdges()
    #------------------------------------------------------------------------
    #------------------------------------------------------------------------ 
    # Graph Object: Set / Retrieve
    #------------------------------------------------------------------------
    #------------------------------------------------------------------------ 

    

    # Retrieve the internal graph object.  
    # This may be desirable after a read operation.
    def getInnerGraph(self):
        return self.G
    
    # Set the internal graph object.  This may be 
    # desirable if an oxbridgeGraph has been created 
    # elsewhere and you want to write it out
    def setInnerGraph(self, Gin):
        self.G = Gin
        

    #------------------------------------------------------------------------
    #------------------------------------------------------------------------ 
    # GraphML Writing interface
    #------------------------------------------------------------------------
    #------------------------------------------------------------------------ 

    # fixes the `writes all to one line' standard behavior of the 
    #   ElementTree object. c.f.
    #   https://stackoverflow.com/questions/3095434/inserting-newlines-in-xml-file-generated-via-xml-etree-elementtree-in-python
    def __indent(self, elem, level=0):
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self.__indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i


    def writeHeader(self):
        self.__graphml = None
        self.__graph = None

        self.__graphml = ET.Element('graphml')
    
    
        # ----------------
        #  Standardized header
        # ----------------
        self.__graphml.set("xmlns","http://graphml.graphdrawing.org/xmlns")
        self.__graphml.set("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
        self.__graphml.set("xsi:schemaLocation","http://graphml.graphdrawing.org/xmlns"\
                +" http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd")

        # -- OxMBM Standardized key selection
        #<key attr.name="FA_std" attr.type="double" for="edge" id="d13" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","FA_std")
        akey.set("attr.type","double")
        akey.set("for","edge")
        akey.set("id", self.edgeAssociations["fastd"])  
  
        #<key attr.name="fiber_length_mean" attr.type="double" for="edge" id="d12" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","fiber_length_mean")
        akey.set("attr.type","double")
        akey.set("for","edge")
        akey.set("id", self.edgeAssociations["lfibermean"])
  
        #<key attr.name="fiber_length_std" attr.type="double" for="edge" id="d11" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","fiber_length_std")
        akey.set("attr.type","double")
        akey.set("for","edge")
        akey.set("id", self.edgeAssociations["lfiberstd"])

        #<key attr.name="FA_mean" attr.type="double" for="edge" id="d10" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","FA_mean")
        akey.set("attr.type","double")
        akey.set("for","edge")
        akey.set("id", self.edgeAssociations["famean"])
    
        #<key attr.name="number_of_fibers" attr.type="int" for="edge" id="d9" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","number_of_fibers")
        akey.set("attr.type","double")
        akey.set("for","edge")
        akey.set("id", self.edgeAssociations["nfibers"])
  
        #<key attr.name="dn_hemisphere" attr.type="string" for="node" id="d7" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_hemisphere")
        akey.set("attr.type","string")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["hemisphere"])

        #<key attr.name="dn_name" attr.type="string" for="node" id="d6" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_name")
        akey.set("attr.type","string")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["name"])

        #<key attr.name="dn_fsname" attr.type="string" for="node" id="d5" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_fsname")
        akey.set("attr.type","string")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["fsname"])

        #<key attr.name="dn_region" attr.type="string" for="node" id="d4" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_region")
        akey.set("attr.type","string")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["region"])

        #<key attr.name="dn_correspondence_id" attr.type="string" for="node" id="d3" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_correspondence_id")
        akey.set("attr.type","string")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["node-id"])

        #<key attr.name="dn_position_z" attr.type="double" for="node" id="d2" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_position_z")
        akey.set("attr.type","double")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["z-coord"])

        #<key attr.name="dn_position_y" attr.type="double" for="node" id="d1" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_position_y")
        akey.set("attr.type","double")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["y-coord"])

        #<key attr.name="dn_position_x" attr.type="double" for="node" id="d0" />
        akey = ET.SubElement(self.__graphml,'key')
        akey.set("attr.name","dn_position_x")
        akey.set("attr.type","double")
        akey.set("for","node")
        akey.set("id", self.nodeAssociations["x-coord"])    
        #---END---Standard Key Selection
        #------------------------------------------------------------
      
        #Create the main graph subnode of graphml
        self.__graph = ET.SubElement(self.__graphml,'graph')
        self.__graph.set("edgedefault","undirected")
        

    
    def writeNodes(self):
        if self.__graphml == None or self.__graph == None:
            print("**Error: Call writeHeader before calling writeNodes")
            return
        

        nodelist = self.G.getNodeIDs()
        
        for nd in nodelist:
            thisnode = self.G.getNode(nd)
            
            nid = thisnode.getID()
            ncoors = thisnode.getCoords()
            
            
            xmlnode = ET.SubElement(self.__graph,'node')
            xmlnode.set("id", str(nid))

            # extract the information from the node
            for ky in self.nodeAssociations.keys():
                
                fval = self.nodeAssociations[ky]
                nval = None
                
                if fval == self.nodeAssociations["x-coord"]:
                    nval = ncoors[0]
                elif fval == self.nodeAssociations["y-coord"]:
                    nval = ncoors[1]
                elif fval == self.nodeAssociations["z-coord"]:
                    nval = ncoors[2]
                elif fval == self.nodeAssociations["node-id"]:
                    nval = nid
                elif fval == self.nodeAssociations["name"]:
                    nval = thisnode.getName()
                elif fval == self.nodeAssociations["region"]:
                    nval = thisnode.getRegion()
                elif fval == self.nodeAssociations["fsname"]:
                    nval = thisnode.getFSName()
                elif fval == self.nodeAssociations["hemisphere"]:
                    nval = thisnode.getHemisphere()
                    
                akey = ET.SubElement(xmlnode, 'data')
                akey.set("key", fval)
                akey.text = str(nval)
        
        
        
    def writeEdges(self):
        if self.__graphml == None or self.__graph == None:
            print("**Error: Call writeHeader before calling writeEdges")
            return

        edgelist = self.G.getEdgeIDs()
        
        # Create all the edges
        for ed in edgelist: 
            thisedge = self.G.getEdge(ed)
            
            # get the various edge-related details
            srcid = ed[0]
            trgid = ed[1]
            eweights = thisedge.getEdgeWeights()
            favals = thisedge.getFAVals()
            flvals = thisedge.getFiberLengthVals()
            
            xmledge = ET.SubElement(self.__graph,'edge')
            xmledge.set("source", str(srcid))
            xmledge.set("target", str(trgid))
        
            # walk the edge associations list
            for ky in self.edgeAssociations.keys():
                fval = self.edgeAssociations[ky]
                nval = None
                
                if fval == self.edgeAssociations["nfibers"]:
                    nval = eweights[0]
                elif fval == self.edgeAssociations["lfibermean"]:
                    nval = eweights[1]
                elif fval == self.edgeAssociations["lfiberstd"]:
                    nval = flvals[1]
                elif fval == self.edgeAssociations["famean"]:
                    nval = favals[0]
                elif fval == self.edgeAssociations["fastd"]:
                    nval = favals[1]                

                akey = ET.SubElement(xmledge,'data')
                akey.set("key", fval)
                akey.text = str(nval)



    def writeFooter(self):
        if self.__graphml == None or self.__graph == None:
            print("**Error: Call writeHeader before calling writeFooter")
            return
           
        # we call the `indent' helper function before writing to go through 
        #   and give us proper newline characters after each key etc
        self.__indent(self.__graphml)

        xmldat = ET.tostring(self.__graphml, encoding="unicode")
        xmlof = open(self.outfile, "w")
        xmlof.write(xmldat)
        xmlof.close()    

        #reset the element tree objects for later use
        self.__graphml = None
        self.__graph = None
                   
        
    def writeAll(self):
        self.writeHeader()
        self.writeNodes()
        self.writeEdges()
        self.writeFooter()
        
        