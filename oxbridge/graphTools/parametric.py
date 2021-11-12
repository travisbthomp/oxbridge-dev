#=============================================================================
# ** Oxford Mathematical Brain Modeling Group **
#   This source file defines some surfaces that may be handy for generating 
#   graphs.  Some of these objects are also reflected in the capability of 
#   automatic graph builders in the autoGraphBuilders.py classes
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

import numpy as np
import math

# this class defines a parametric function of a single variable.
# e.g. x(t).  It is assumed that this is a real-valued function 
# of the form x : [a,b] --> R.  This function will be transformed to 
# of the form y : [0,1] --> R with y(0) = x(a), y(1) = x(b) etc. Thus 
# providing a standard domain across all parameteric functions   
class parametricFunction:
    
    # Inputs
    # fcn: a handle to a function which takes a floating point input
    # [optional]:
    #   lowerbd: lower bound of the parametric interval
    #   upperbd: upper bound of the parametric interval
    #
    def __init__(self, fcn, lowerbd=0, upperbd=1):
        
        if lowerbd == upperbd:
            print("** The parametric lower bound cannot coincide with the upper bound")
            print("** Using default values of lowerbd=0, upperbd=1")
            lowerbd=0
            upperbd=1
            
        self.lb = lowerbd
        self.ub = upperbd
        self.fcn = fcn
        
        # Maps T in [lowerbd, upperbd] to t in [0,1]
        self.phi = lambda T : (T-self.lb) / (self.ub - self.lb)
        
        # Maps t in [0,1] to T=phi(t) in [lowerbd, upperbd]
        self.phiinv = lambda t : t*self.ub + (1-t)*self.lb

    
    # evaluate the function between 0 and 1.  
    # [optional]
    #    set original=True to indicate that tin is a value in [lowerbd,upperbd]
    def eval(self, tin, original=False):
        
        # if original==False then we assume tin is a value in [0,1]
        if original == False and (tin < 0 or tin > 1):
            print("**Standard parametric evaluation: requires a value in [0,1]")
            print("**call eval(..) with convert=True to convert values")
            print(f"** in [{self.lb},{self.ub}] to [0,1]")
            return None

        T = tin
        
        # tin is a value in [0,1] so we convert it to the original 
        # bounds expected by the function passed in at construction
        if original == False:
            T = self.phiinv(tin)
        
        v = self.fcn(T)
        return v

    # evaluate the function at the values indicated in the numpy array npar
    # input: npar - a list or numpy array of the points at which to evaluate the 
    #   parametric function.  
    # [optional]
    #   set original=True to indicate that the points in the numpy array
    #   `npar' come from the interval [lowerbd,upperbd] instead of [0,1]
    def evalat(self, npar, original=False):
        res = []
        for tv in npar:
            res.append(self.eval(tv,original))
        return np.array(res)
        
    # generates a list of N equally spaced points a distance of 
    # dt = 1/N apart
    def evalNPoints(self,N):
        linsp = np.linspace(0, 1, num=N)
        retv = self.evalat(linsp)
        return retv



# a parametric curve embedded in 3D space
class parametric3DCurve:
    
    def __init__(self, xp, yp, zp):
        self.x = xp
        self.y = yp
        self.z = zp

    # input a list of tuples [(x1,y1,z1),...,(xN,yN,zN)]
    # this function checks to see if the first and last points coincide
    # and if fix=True, removes the last point from the list.
    # [optional]:
    #   set fix = False to keep all points in the list
    #
    # returns:
    # 1. The (potentially modified) input list
    # 2. True: if the curve is closed, False if not
    def __closed(self, lst, fix=True, eps=1e-12):
        pta = lst[0]
        ptb = lst[len(lst)-1]
        retl = lst
        
        # we check the distance between the first and last point
        L = math.sqrt( (pta[0]-ptb[0])**2 + (pta[1]-ptb[1])**2 + (pta[2]-ptb[2])**2 )
        bClosed = (L <= eps)
        
        # if the loop is closed and we want 
        # to fixed closed loops, delete the 
        # last entry in the list
        if bClosed and fix:
            del retl[-1]
        
        return retl, bClosed
    
    # turn lists into tuples
    def __tup(self, x, y, z):
        retlist = []
        lval = len(x)
        
        for idx in range(lval):
            retlist.append( (x[idx], y[idx], z[idx]) )
        
        return retlist


    # the curve analogue of the parametric evalat function.
    # [optional]
    #   set fixclosed = False to keep all points, even if 
    #   this is a closed curve
    def evalat(self, npar, original=False, fixclosed=True):
        xar = self.__xp.evalat(npar, original)
        yar = self.__yp.evalat(npar, original)
        zar = self.__zp.evalat(npar, original)
        
        retl = self.__tup(xar, yar, zar)
        retl, closedCurve = self.__closed(retl, fix=fixclosed)

        return retl, closedCurve

    # returns N evenly-spaced coordinates 
    # (xcoord,ycoord,zcoord) along the path.  
    # if the path is closed, so that 
    # (x(0), y(0), z(0)) = (x(1),y(1),z(1)), 
    # returns N-1 points and a true/false value
    # indicating if the curve is a closed curve
    def evalNPoints(self, N):
        xpts = self.x.evalNPoints(N)
        ypts = self.y.evalNPoints(N)
        zpts = self.z.evalNPoints(N)

        retlist = self.__tup(xpts, ypts, zpts)
        retlist, closedCurve = self.__closed(retlist)
        
        return retlist, closedCurve










